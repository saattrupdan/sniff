import numpy as np
import pytest

from sniff import formula_id, ion_roles


class Catalogue:
    def enrich_candidates(self, candidates):
        enriched = []
        for candidate in candidates:
            item = dict(candidate)
            if item["formula"] == "C2H6O":
                item["catalogue"] = [{"name": "ethanol"}]
            enriched.append(item)
        return enriched


def peak(mz, height=100.0):
    return {
        "mz": mz,
        "height": height,
        "formula_tolerance": {
            "ppm": 10.0,
            "proposal_ppm": 200.0,
            "score_sigma_ppm": 4.0,
            "mDa": mz * 10.0 / 1000.0,
        },
        "candidates": [],
    }


def test_dehydrated_ion_gets_separate_named_candidate():
    neutral_mass = formula_id.formula_mass({"C": 2, "H": 6, "O": 1})
    observed = neutral_mass + formula_id.PROTON - ion_roles.WATER_MASS
    item = peak(observed)

    ion_roles.annotate_ion_candidates(
        [item], Catalogue(), {"reagent": "H3O+", "status": "available"}
    )

    candidate = item["ion_candidates"][0]
    assert candidate["formula"] == "C2H6O"
    assert candidate["ion_notation"] == "[M+H-H2O]+"
    assert candidate["mass_eligible"] is True
    assert candidate["assignment_eligible"] is False
    assert candidate["catalogue"][0]["name"] == "ethanol"


def test_hydrated_search_uses_observed_ion_ppm_window():
    neutral_mass = formula_id.formula_mass({"C": 2, "H": 6, "O": 1})
    theoretical = neutral_mass + formula_id.PROTON + ion_roles.WATER_MASS
    item = peak(theoretical * (1.0 + 180e-6))

    ion_roles.annotate_ion_candidates(
        [item], Catalogue(), {"reagent": "H3O+", "status": "available"}
    )

    candidate = next(
        value
        for value in item["ion_candidates"]
        if value["ion_kind"] == "hydrated-protonated"
        and value["formula"] == "C2H6O"
    )
    assert candidate["mass_eligible"] is False
    assert candidate["delta_ppm"] == pytest.approx(180.0, abs=0.01)


def test_charge_transfer_requires_measured_o2_marker():
    neutral_mass = formula_id.formula_mass({"C": 2, "H": 6, "O": 1})
    without_marker = peak(neutral_mass - ion_roles.ELECTRON_MASS)
    with_marker = peak(neutral_mass - ion_roles.ELECTRON_MASS)

    ion_roles.annotate_ion_candidates(
        [without_marker, peak(31.989)],
        Catalogue(),
        {"reagent": "H3O+", "status": "available"},
    )
    ion_roles.annotate_ion_candidates(
        [with_marker],
        Catalogue(),
        {"reagent": "O2+", "status": "available"},
    )

    assert all(
        candidate["ion_kind"] != "charge-transfer"
        for candidate in without_marker.get("ion_candidates", [])
    )
    assert any(
        candidate["ion_kind"] == "charge-transfer"
        and candidate["formula"] == "C2H6O"
        and candidate["mass_eligible"] is True
        and candidate["delta_ppm"] == pytest.approx(0.0)
        for candidate in with_marker["ion_candidates"]
    )


def test_hydride_abstraction_exact_mass_includes_missing_electron():
    neutral_mass = formula_id.formula_mass({"C": 2, "H": 6, "O": 1})
    observed = neutral_mass - ion_roles.HYDROGEN_MASS - ion_roles.ELECTRON_MASS
    item = peak(observed)

    ion_roles.annotate_ion_candidates(
        [item],
        Catalogue(),
        {"reagent": "NO+", "status": "available"},
    )

    candidate = next(
        value
        for value in item["ion_candidates"]
        if value["ion_kind"] == "hydride-abstraction"
        and value["formula"] == "C2H6O"
    )
    assert candidate["mass_eligible"] is True
    assert candidate["delta_ppm"] == pytest.approx(0.0)


def test_two_covarying_satellites_support_doubly_charged_candidate():
    neutral_mass = formula_id.formula_mass({"C": 8, "H": 18, "O": 2})
    parent_mz = (neutral_mass + 2 * formula_id.PROTON) / 2
    masses = [
        parent_mz,
        parent_mz + formula_id.DM1 / 2,
        parent_mz + formula_id.DM1,
    ]
    peaks = [peak(masses[0], 100.0), peak(masses[1], 10.0), peak(masses[2], 2.0)]
    base = np.linspace(1.0, 20.0, 40)
    traces = {
        masses[0]: base,
        masses[1]: base * 0.1,
        masses[2]: base * 0.02,
    }

    ion_roles.annotate_ion_candidates(
        peaks,
        Catalogue(),
        {"reagent": "H3O+", "status": "available"},
        traces=traces,
    )

    charged = [
        candidate
        for candidate in peaks[0]["ion_candidates"]
        if candidate["ion_kind"] == "multiply-charged-protonated"
    ]
    assert charged
    assert charged[0]["charge"] == 2
    assert charged[0]["formula"] == "C8H18O2"
    assert len(charged[0]["charge_evidence"]["satellites"]) == 2


def test_one_satellite_cannot_create_charge_state_candidate():
    parent = peak(75.0, 100.0)
    satellite = peak(75.0 + formula_id.DM1 / 2, 10.0)
    base = np.linspace(1.0, 20.0, 40)

    ion_roles.annotate_ion_candidates(
        [parent, satellite],
        Catalogue(),
        {"reagent": "H3O+", "status": "available"},
        traces={parent["mz"]: base, satellite["mz"]: base * 0.1},
    )

    assert all(
        candidate["ion_kind"] != "multiply-charged-protonated"
        for candidate in parent.get("ion_candidates", [])
    )


def test_coverage_categories_partition_the_peak_count():
    peaks = [
        {"candidates": [{"formula": "CH2O", "name": "formaldehyde"}]},
        {"ion_candidates": [{"formula": "C2H6O", "catalogue": [{"name": "ethanol"}]}]},
        {
            "interpretation_candidates": [
                {"kind": "isotope", "candidate_formulas": ["C3H6O"]}
            ]
        },
        {"interpretation_candidates": [{"kind": "reagent"}]},
        {
            "interpretation_candidates": [
                {"kind": "authored", "formula": "C7H14"}
            ]
        },
        {"interpretation_candidates": [{"kind": "unresolved"}]},
    ]

    summary = ion_roles.coverage_summary(peaks)

    assert summary["with_formula_or_linked_candidate"] == 3
    assert summary["with_named_compound_candidate"] == 2
    assert summary["detected_peaks"] == len(peaks)
    assert summary["total_components"] == len(peaks)
    assert sum(summary["categories"].values()) == len(peaks)
    assert summary["candidate_coverage_percent"] == pytest.approx(50.0)
    assert summary["categories"]["authored_assignment"] == 1
    assert summary["signal_weighted_candidate_coverage_percent"] is None
    assert summary["signal_weight_source"].endswith("status is missing")


def test_signal_weighted_coverage_is_unavailable_for_withheld_component():
    peaks = [
        {
            "mz": 50.0,
            "role_signal": None,
            "role_trace_status": "unresolved",
            "interpretation_candidates": [{"kind": "unresolved-overlap"}],
        }
    ]

    summary = ion_roles.coverage_summary(peaks)

    assert summary["signal_weighted_candidate_coverage_percent"] is None
    assert summary["signal_weight_source"].startswith("unavailable")
    assert summary["categories"]["unknown"] == 1
    assert summary["categories"]["interpreted_noncompound"] == 0


@pytest.mark.parametrize("signal", [float("nan"), float("inf"), float("-inf")])
def test_nonfinite_signal_makes_weighted_coverage_unavailable(signal):
    summary = ion_roles.coverage_summary(
        [
            {
                "mz": 50.0,
                "role_signal": signal,
                "role_trace_status": "available",
                "candidates": [{"formula": "CH2O"}],
            }
        ]
    )

    assert summary["signal_weighted_candidate_coverage_percent"] is None
    assert summary["signal_weight_source"].startswith("unavailable")


def test_mixed_finite_and_nonfinite_component_withholds_weighted_coverage():
    peaks = [
        {
            "mz": 50.0,
            "role_signal": 10.0,
            "role_trace_status": "available",
            "candidates": [{"formula": "CH2O"}],
            "overlap": {"level": "unresolved", "neighbor": 50.01},
        },
        {
            "mz": 50.01,
            "role_signal": float("nan"),
            "role_trace_status": "available",
            "overlap": {"level": "unresolved", "neighbor": 50.0},
        },
    ]

    summary = ion_roles.coverage_summary(peaks)

    assert summary["total_components"] == 1
    assert summary["signal_weighted_candidate_coverage_percent"] is None
    assert summary["signal_weight_source"].startswith("unavailable")


def test_coverage_collapses_unresolved_neighbours_and_weights_once():
    peaks = [
        {
            "mz": 50.0,
            "height": 100.0,                "role_signal": 80.0,
                "role_trace_status": "available",
                "candidates": [{"formula": "C3H13"}],
            "overlap": {"level": "unresolved", "neighbor": 50.015},
        },
        {
            "mz": 50.015,
            "height": 90.0,                "role_signal": 70.0,
                "role_trace_status": "available",
                "overlap": {"level": "unresolved", "neighbor": 50.0},
        },
        {
            "mz": 75.0,
            "height": 20.0,                "role_signal": 20.0,
                "role_trace_status": "available",
                "interpretation_candidates": [{"kind": "unresolved"}],
        },
    ]

    summary = ion_roles.coverage_summary(peaks)

    assert summary["detected_peaks"] == 3
    assert summary["total_components"] == 2
    assert summary["with_formula_or_linked_candidate"] == 1
    assert summary["candidate_coverage_percent"] == pytest.approx(50.0)
    assert summary["signal_weighted_candidate_coverage_percent"] == pytest.approx(80.0)
