from pathlib import Path

from sniff import formula_id
from sniff.analyze import (
    _assign_suggested_identities,
    annotate_peaks,
    interpret_peak_roles,
)


def test_known_reagent_gets_non_analyte_interpretation():
    peaks = [{"mz": 29.998, "height": 100.0, "candidates": []}]

    interpret_peak_roles(peaks)

    interpretation = peaks[0]["interpretation_candidates"][0]
    assert interpretation["kind"] == "reagent"
    assert interpretation["label"] == "NO+"
    assert interpretation["exclude_from_analyte_assignment"] is True


def test_isobaric_water_cluster_does_not_suppress_valid_analyte_candidate():
    peak = {
        "mz": 73.056,
        "height": 100.0,
        "candidates": [{"formula": "C4H8O", "iso_pred": [0.045, 0.002]}],
    }

    interpret_peak_roles([peak])

    assert "interpretation_candidates" not in peak


def test_stronger_parent_supports_possible_isotope_interpretation():
    parent = {
        "mz": 55.0,
        "height": 100.0,
        "candidates": [{"formula": "C5H10", "iso_pred": [0.055, 0.002]}],
    }
    child = {
        "mz": 55.0 + formula_id.DM1 + 0.0119,
        "height": 10.0,
        "candidates": [],
    }

    unrelated = {
        "mz": 80.0,
        "height": 200.0,
        "candidates": [{"formula": "C2H8O3", "iso_pred": [0.02, 0.001]}],
    }

    interpret_peak_roles([parent, child, unrelated])

    interpretation = child["interpretation_candidates"][0]
    assert interpretation["kind"] == "isotope"
    assert interpretation["isotope_order"] == 1
    assert interpretation["related_mz"] == 55.0
    assert "spacing residual +11.9 mDa" in interpretation["evidence"]
    assert "observed ratio 0.1; predicted 0.055" in interpretation["evidence"]
    assert interpretation["candidate_formulas"] == ["C5H10"]


def test_isotope_spacing_outside_twelve_mda_is_not_claimed():
    parent = {
        "mz": 55.0,
        "height": 100.0,
        "candidates": [{"formula": "C5H10", "iso_pred": [0.055, 0.002]}],
    }
    child = {
        "mz": 55.0 + formula_id.DM1 + 0.0121,
        "height": 10.0,
        "candidates": [],
    }

    interpret_peak_roles([parent, child])

    assert child["interpretation_candidates"][0]["kind"] == "unresolved"


def test_valid_formula_candidate_is_not_reclassified_as_an_isotope():
    parent = {
        "mz": 45.0335,
        "height": 100.0,
        "candidates": [{"formula": "C2H4O", "iso_pred": [0.023, 0.002]}],
    }
    independent = {
        "mz": 47.0491,
        "height": 5.0,
        "candidates": [{"formula": "C2H6O", "iso_pred": [0.023, 0.002]}],
    }

    interpret_peak_roles([parent, independent])

    assert "interpretation_candidates" not in independent


def test_supported_fragment_link_is_an_interpretation_not_an_identity():
    fragment = {
        "mz": 43.018,
        "height": 20.0,
        "candidates": [],
        "fragmentation_links": [
            {
                "parent_mz": 73.065,
                "parent_formula": "C4H8O",
                "candidate_name": "test ketone",
                "expected_fragment_mz": 43.018,
                "level_correlation": 0.95,
                "change_correlation": 0.88,
            },
            {
                "parent_mz": 87.080,
                "parent_formula": "C5H10O",
                "candidate_name": "second ketone",
                "expected_fragment_mz": 43.018,
                "level_correlation": 0.90,
                "change_correlation": 0.80,
            },
        ],
    }

    interpret_peak_roles([fragment])

    interpretation = fragment["interpretation_candidates"][0]
    assert interpretation["kind"] == "fragment"
    assert interpretation["exclude_from_analyte_assignment"] is True
    assert "possible fragment" in interpretation["label"]
    assert any("not MS/MS proof" in item for item in interpretation["evidence"])
    assert interpretation["candidate_formulas"] == ["C4H8O", "C5H10O"]
    assert interpretation["compound_candidates"] == ["second ketone", "test ketone"]


def test_peak_without_formula_candidate_gets_explicit_unresolved_interpretation():
    peak = {
        "mz": 41.0541,
        "height": 100.0,
        "candidates": [],
        "overlap": {"neighbor": 41.0799},
    }

    interpret_peak_roles([peak])

    interpretation = peak["interpretation_candidates"][0]
    assert interpretation["kind"] == "unresolved"
    assert interpretation["label"] == "unresolved ion at m/z 41.0541"
    assert interpretation["evidence"] == [
        "no plausible protonated-neutral formula fits the current exact-mass tolerance",
        "the peak overlaps a neighbouring channel",
    ]


def test_interpretation_prevents_automatic_analyte_assignment():
    peaks = [
        {
            "mz": 56.0034,
            "height": 10.0,
            "candidates": [
                {
                    "formula": "C3H3O",
                    "name": "example",
                    "preferred_name": "example",
                    "probability": 1.0,
                }
            ],
            "id_confidence": 1.0,
            "interpretation_candidates": [
                {
                    "kind": "isotope",
                    "label": "possible M+1 isotope of m/z 55.0000",
                    "exclude_from_analyte_assignment": True,
                }
            ],
        }
    ]

    _assign_suggested_identities(peaks, assign_all_library=True)

    assert peaks[0]["suggested_label"] == "possible M+1 isotope of m/z 55.0000"
    assert "suggested_formula" not in peaks[0]


def test_alternative_ion_candidate_is_visible_but_not_auto_assignable():
    peak = {
        "mz": 47.0,
        "height": 20.0,
        "candidates": [],
        "ion_candidates": [
            {
                "formula": "C2H6O",
                "ion_notation": "[M+H+H2O]+",
                "pathway_reason": "measured H3O+ context",
                "delta_ppm": 2.0,
                "mass_match": "alternative-ion-within-run-tolerance",
                "limitation": "ion pathway is not established",
                "catalogue": [{"name": "ethanol"}],
                "assignment_eligible": False,
            }
        ],
    }

    interpret_peak_roles([peak])

    interpretation = peak["interpretation_candidates"][0]
    assert interpretation["kind"] == "alternative-ion"
    assert interpretation["candidate_formulas"] == ["C2H6O"]
    assert interpretation["compound_candidates"] == ["ethanol"]
    assert interpretation["exclude_from_analyte_assignment"] is True


def test_ptr_fixture_has_candidate_or_interpretation_for_every_peak():
    path = Path(__file__).parent / "fixtures" / "ptr_peak_masses.txt"
    masses = [
        float(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    assert len(masses) == 182

    _drift, peaks = annotate_peaks(
        [{"mz": mass, "height": 1.0, "rel_height": 1.0} for mass in masses],
        assign_all_library=True,
    )

    assert sum(bool(peak["candidates"]) for peak in peaks) >= 35
    assert all(
        peak.get("candidates") or peak.get("interpretation_candidates")
        for peak in peaks
    )
    assigned = [peak["suggested_formula"].upper() for peak in peaks if peak.get("suggested_formula")]
    assert len(assigned) == len(set(assigned))
