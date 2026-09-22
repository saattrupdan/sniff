"""Evidence-gated detector ringing diagnostics."""

import numpy as np
from calibration_helpers import identity_mass_axis

from sniff import ringing
from sniff.analyze import (
    _annotate_timebin_ringing,
    _is_noise_artifact,
    _ringing_artifact_diagnostics,
)


def _mass_after_delay(axis, mass, delay):
    return float(axis.tb_to_m(float(axis.m_to_tb(mass)) + delay))


def test_repeated_timebin_echo_requires_supported_parent_trace():
    axis = identity_mass_axis(a=10_000.0, b=0.0)
    cycles = np.arange(256, dtype=np.float64)
    trace = 100.0 + 30.0 * np.sin(cycles / 17.0) + 0.08 * cycles
    peaks = []
    traces = {}
    for parent_mass in (30.0, 35.0, 40.0, 45.0):
        child_mass = _mass_after_delay(axis, parent_mass, 309.0)
        peaks.extend(
            [
                {
                    "mz": parent_mass,
                    "height": 1_000.0,
                    "prominence": 900.0,
                    "candidates": [],
                },
                {
                    "mz": child_mass,
                    "height": 10.0,
                    "prominence": 9.0,
                    "candidates": [],
                },
            ]
        )
        traces[parent_mass] = trace
        traces[child_mass] = trace * 0.01

    supported_parent = 60.0
    supported_parent_trace = 60.1
    supported_child = _mass_after_delay(axis, supported_parent, 309.0)
    supported_child_trace = supported_child + 0.1
    unsupported_parent = 70.0
    unsupported_child = _mass_after_delay(axis, unsupported_parent, 309.0)
    peaks.extend(
        [
            {
                "mz": supported_parent_trace,
                "apex": supported_parent,
                "height": 1_000.0,
                "prominence": 900.0,
                "candidates": [],
            },
            {
                "mz": supported_child_trace,
                "apex": supported_child,
                "height": 300.0,
                "prominence": 250.0,
                "candidates": [],
            },
            {
                "mz": unsupported_parent,
                "height": 1_000.0,
                "prominence": 900.0,
                "candidates": [],
            },
            {
                "mz": unsupported_child,
                "height": 300.0,
                "prominence": 250.0,
                "candidates": [],
            },
        ]
    )
    traces[supported_parent_trace] = trace
    traces[supported_child_trace] = trace * 0.01
    traces[unsupported_parent] = trace
    traces[unsupported_child] = trace[::-1]

    _annotate_timebin_ringing(peaks, traces, axis)

    supported = next(peak for peak in peaks if peak.get("apex") == supported_child)
    unsupported = next(peak for peak in peaks if peak["mz"] == unsupported_child)
    assert supported["artifact_evidence"]["status"] == "likely"
    assert "ringing echo" in supported["likely_artifact"][0]
    assert supported["artifact_evidence"]["model"].endswith("ringing-v2")
    assert supported["artifact_evidence"]["response_model"]["parent_count"] >= 3
    assert supported["artifact_evidence"]["response_metrics"]["held_out_nrmse"] < 0.01
    assert supported["artifact_evidence"]["mode_seed_parent_count"] >= 3
    assert supported["artifact_evidence"]["parent_mz"] == supported_parent
    assert supported["artifact_evidence"]["parent_trace_mz"] == supported_parent_trace
    assert unsupported["artifact_evidence"]["status"] == "supporting"
    assert not unsupported.get("likely_artifact")


def test_response_model_rejects_inconsistent_candidate_gain():
    cycles = np.arange(160, dtype=np.float64)
    parent = 100.0 + 20.0 * np.sin(cycles / 13.0)
    references = []
    for index, gain in enumerate((0.0098, 0.0100, 0.0102, 0.0101)):
        metrics = ringing.response_metrics(parent, gain * parent)
        references.append(
            {
                "parent_id": index,
                "delay": 309.0 + 0.2 * index,
                "metrics": metrics,
            }
        )
    model = ringing.fit_response_model(references)
    candidate = ringing.response_metrics(parent, 0.03 * parent)
    score = ringing.score_response(candidate, model)

    assert model["reliable"]
    assert not score["supported"]
    assert any("prediction interval" in gate for gate in score["failed_gates"])


def test_response_model_counts_only_individually_valid_parents():
    cycles = np.arange(160, dtype=np.float64)
    parent = 100.0 + 20.0 * np.sin(cycles / 13.0)
    references = []
    for index in range(2):
        references.append(
            {
                "parent_id": index,
                "delay": 309.0,
                "metrics": ringing.response_metrics(parent, 0.01 * parent),
            }
        )
    bad = ringing.response_metrics(parent, np.roll(0.01 * parent, 40))
    references.append({"parent_id": 2, "delay": 309.0, "metrics": bad})

    model = ringing.fit_response_model(references)

    assert not model["reliable"]
    assert model["reason"] == "fewer than three independent response parents"


def test_connected_echo_chain_counts_as_one_response_family():
    metrics = {
        "usable": True,
        "held_out_nrmse": 0.01,
        "gain_relative_mad": 0.01,
        "gain": 0.01,
    }
    pairs = [
        {"parent_id": 1, "child_id": 2, "delay": 309.0, "metrics": metrics},
        {"parent_id": 2, "child_id": 3, "delay": 309.0, "metrics": metrics},
        {"parent_id": 3, "child_id": 4, "delay": 309.0, "metrics": metrics},
    ]
    families = ringing.response_families(pairs)
    for pair in pairs:
        pair["family_id"] = families[pair["parent_id"]]

    model = ringing.fit_response_model(pairs)

    assert len(set(families.values())) == 1
    assert not model["reliable"]
    assert "three independent" in model["reason"]


def test_ringing_summary_counts_version_two_evidence():
    summary = _ringing_artifact_diagnostics(
        [
            {
                "artifact_evidence": {
                    "model": "calibrated-timebin-ringing-v2",
                    "status": "likely",
                }
            },
            {
                "artifact_evidence": {
                    "model": "calibrated-timebin-ringing-v2",
                    "status": "supporting",
                }
            },
            {
                "artifact_evidence": {
                    "model": "calibrated-timebin-ringing-v1",
                    "status": "likely",
                }
            },
        ]
    )

    assert summary["model"] == "calibrated-timebin-ringing-v2"
    assert summary["likely"] == 1
    assert summary["supporting_only"] == 1


def test_mass_only_shoulder_hint_is_not_default_filtered():
    assert not _is_noise_artifact(["possible high-side shoulder of taller m/z 41.067"])


def test_timebin_pattern_does_not_override_formula_candidate():
    axis = identity_mass_axis(a=10_000.0, b=0.0)
    trace = np.linspace(1.0, 10.0, 256) ** 2
    peaks = []
    traces = {}
    for parent_mass in (30.0, 35.0, 40.0, 45.0):
        child_mass = _mass_after_delay(axis, parent_mass, 309.0)
        peaks.extend(
            [
                {
                    "mz": parent_mass,
                    "height": 1_000.0,
                    "prominence": 900.0,
                    "candidates": [],
                },
                {
                    "mz": child_mass,
                    "height": 10.0,
                    "prominence": 9.0,
                    "candidates": [],
                },
            ]
        )
        traces[parent_mass] = trace
        traces[child_mass] = trace * 0.01
    candidate_parent = 60.0
    candidate_child = _mass_after_delay(axis, candidate_parent, 309.0)
    peaks.extend(
        [
            {
                "mz": candidate_parent,
                "height": 1_000.0,
                "prominence": 900.0,
                "candidates": [],
            },
            {
                "mz": candidate_child,
                "height": 300.0,
                "prominence": 250.0,
                "candidates": [{"formula": "C4H10", "assignment_eligible": True}],
            },
        ]
    )
    traces[candidate_parent] = trace
    traces[candidate_child] = trace * 0.3

    _annotate_timebin_ringing(peaks, traces, axis)

    candidate = next(peak for peak in peaks if peak["mz"] == candidate_child)
    assert "artifact_evidence" not in candidate
    assert not candidate.get("likely_artifact")
