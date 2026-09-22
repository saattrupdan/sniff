from unittest import mock

from sniff import diagnostics


class DiagnosticCatalogue:
    def score_peak(self, mz, **_kwargs):
        if mz < 55.0:
            return [
                {
                    "formula": "C4H2",
                    "ion_mz": 50.0,
                    "delta_ppm": 480.0,
                    "delta_mDa": 24.0,
                }
            ]
        return []


def unresolved_peak(mz, *, neighbour=None, trace_status="available"):
    peak = {
        "mz": mz,
        "apex": mz,
        "height": 10.0,
        "candidates": [],
        "role_trace_status": trace_status,
        "formula_tolerance": {
            "ppm": 10.0,
            "proposal_ppm": 200.0,
        },
        "ion_role": {"kind": "unresolved"},
        "interpretation_candidates": [{"kind": "unresolved"}],
    }
    if neighbour is not None:
        peak["overlap"] = {"level": "unresolved", "neighbor": neighbour}
    return peak


def test_component_report_explains_every_unknown_component_once():
    peaks = [
        unresolved_peak(50.0),
        unresolved_peak(60.0, neighbour=60.01, trace_status="unresolved"),
        unresolved_peak(60.01, neighbour=60.0, trace_status="unresolved"),
        {
            "mz": 70.0,
            "height": 20.0,
            "candidates": [{"formula": "C4H6O"}],
        },
    ]

    with mock.patch.object(
        diagnostics.catalogue, "CompoundCatalogue", return_value=DiagnosticCatalogue()
    ):
        report = diagnostics.component_report(peaks)

    assert report["component_count"] == 2
    assert report["reason_counts"] == {
        "nearest_formula_outside_200_ppm_proposal_gate": 1,
        "no_plausible_formula_within_2000_ppm": 1,
    }
    overlap = next(item for item in report["components"] if len(item["member_mz"]) == 2)
    assert "inseparable_overlap" in overlap["additional_constraints"]
    assert overlap["nearest_formula_diagnostic"] is None
    assert report["components"][0]["route_checks"]["direct_protonated_formula"]


def test_baseline_comparison_preserves_the_prior_unknown_cohort():
    baseline = [
        unresolved_peak(50.0),
        {
            "mz": 70.0,
            "height": 20.0,
            "candidates": [{"formula": "C4H6O"}],
        },
    ]
    current = [
        {
            "mz": 50.0,
            "height": 10.0,
            "candidates": [{"formula": "C4H2", "delta_ppm": 2.0}],
        },
        unresolved_peak(70.0),
    ]

    with mock.patch.object(
        diagnostics.catalogue, "CompoundCatalogue", return_value=DiagnosticCatalogue()
    ):
        report = diagnostics.compare_baseline_unknowns(baseline, current)

    assert report["baseline_component_count"] == 1
    assert report["transition_counts"] == {"direct_protonated_formula": 1}
    assert report["gained_support"] == 1
    assert report["remaining_unknown"] == 0
    assert report["newly_unknown_component_ids"] == [2]
    assert report["components"][0]["baseline"]["component_id"] == 1


def test_component_report_explains_resolved_components_when_requested():
    peak = {
        "mz": 59.0491,
        "apex": 59.0492,
        "height": 20.0,
        "candidates": [
            {
                "formula": "C3H6O",
                "delta_ppm": 1.2,
                "mass_match": "broad-proposal",
                "assignment_eligible": False,
            }
        ],
    }

    with mock.patch.object(
        diagnostics.catalogue, "CompoundCatalogue", return_value=DiagnosticCatalogue()
    ):
        report = diagnostics.component_report([peak], include_resolved=True)

    component = report["components"][0]
    assert component["category"] == "direct_protonated_formula"
    assert component["primary_reason"] == "direct_formula_proposal_within_200_ppm"
    assert component["nearest_formula_diagnostic"] is None
    assert component["route_checks"]["direct_protonated_formula"] == "available"
    assert component["supporting_evidence"] == [peak["candidates"][0]]
