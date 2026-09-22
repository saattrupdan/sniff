"""Run-derived ppm tolerance and exact-mass scoring tests."""

import copy

import h5py
import numpy as np
from calibration_helpers import identity_mass_axis

from sniff import analyze, formula_id, ptrms


def test_ppm_window_scales_in_absolute_mass():
    axis = identity_mass_axis()

    low = ptrms.formula_assignment_tolerance(axis, 50.0)
    high = ptrms.formula_assignment_tolerance(axis, 200.0)

    assert low["ppm"] == high["ppm"] == 5.0
    assert low["mDa"] == 0.25
    assert high["mDa"] == 1.0


def test_formula_candidates_distinguish_assignment_and_proposal_boundaries():
    counts = formula_id.isotopes.parse_formula("C3H6O")
    theoretical = formula_id.formula_mass(counts) + formula_id.PROTON

    included = formula_id.score_peak(
        theoretical * (1.0 + 9e-6),
        1.0,
        elements=["C", "O"],
        tolerance_ppm=10.0,
    )
    provisional = formula_id.score_peak(
        theoretical * (1.0 + 11e-6),
        1.0,
        elements=["C", "O"],
        tolerance_ppm=10.0,
    )
    excluded = formula_id.score_peak(
        theoretical * (1.0 + 201e-6),
        1.0,
        elements=["C", "O"],
        tolerance_ppm=10.0,
    )

    assert included[0]["formula"] == "C3H6O"
    assert included[0]["delta_ppm"] == 9.0
    assert included[0]["assignment_eligible"] is True
    assert provisional[0]["formula"] == "C3H6O"
    assert provisional[0]["assignment_eligible"] is False
    assert provisional[0]["mass_match"] == "broad-proposal-only"
    assert excluded == []


def test_proposal_slots_survive_when_validated_candidate_limit_is_full():
    candidates = formula_id.score_peak(
        99.0,
        1.0,
        tolerance_ppm=10.0,
        proposal_tolerance_ppm=200.0,
        max_candidates=5,
    )

    eligible = [candidate for candidate in candidates if candidate["assignment_eligible"]]
    proposals = [
        candidate for candidate in candidates if not candidate["assignment_eligible"]
    ]
    assert len(eligible) == 5
    assert len(proposals) == 2


def test_formula_enumeration_supports_siloxanes_and_iodinated_references():
    assert formula_id.dbe({"Si": 1, "H": 4}) == 0
    assert formula_id.dbe({"C": 6, "H": 18, "O": 1, "Si": 2}) == 0
    assert formula_id.dbe({"C": 6, "H": 18, "O": 3, "Si": 3}) == 1
    formulas = (
        {"C": 2, "H": 6, "O": 1, "Si": 1},
        {"C": 6, "H": 5, "I": 1},
    )

    for counts in formulas:
        mass = formula_id.formula_mass(counts)
        matches = formula_id.enumerate_formulas(mass, 1e-6)
        assert formula_id.formula_str(counts) in {
            formula_id.formula_str(candidate) for candidate, _mass in matches
        }


def test_high_mass_search_uses_catalogue_instead_of_unbounded_enumeration(monkeypatch):
    def fail(*_args, **_kwargs):
        raise AssertionError("high-mass local enumeration should be skipped")

    monkeypatch.setattr(formula_id, "enumerate_formulas", fail)

    assert formula_id.score_peak(
        500.0,
        1.0,
        proposal_tolerance_ppm=200.0,
    ) == []


def test_held_out_mapping_residuals_set_tolerance_and_degrade_above_limit():
    masses = np.array([21.0221, 203.9430, 330.8480], dtype=np.float64)
    a, b = 10000.0, -200.0

    def model(residuals):
        shifted = masses * (1.0 + np.asarray(residuals) * 1e-6)
        mapping = np.column_stack((masses, a * np.sqrt(shifted) + b))
        with h5py.File("mapping", "w", driver="core", backing_store=False) as h5:
            h5.create_dataset("CALdata/Mapping", data=mapping)
            return ptrms._derive_formula_tolerance_model(h5, a, b)

    accepted = model([0.0, 0.0, 0.0])
    degraded = model([0.0, 20.0, 0.0])

    assert accepted["status"] == "accepted"
    assert 5.0 <= accepted["tolerance_ppm"] <= 10.0
    assert accepted["candidate_generation_allowed"] is True
    assert degraded["status"] == "degraded"
    assert degraded["candidate_generation_allowed"] is True
    assert degraded["automatic_assignment_allowed"] is False
    assert "exceeds the 10 ppm" in degraded["reason"]


def test_fallback_candidates_are_not_assigned_automatically():
    axis = identity_mass_axis()
    diagnostics = copy.deepcopy(axis.diagnostics)
    diagnostics.pop("formula_assignment_tolerance")
    legacy_axis = ptrms.MassAxisCalibration(
        axis.a,
        axis.b,
        scale=axis.scale,
        offset=axis.offset,
        diagnostics=diagnostics,
    )
    counts = formula_id.isotopes.parse_formula("C3H6O")
    theoretical = formula_id.formula_mass(counts) + formula_id.PROTON

    _drift, peaks = analyze.annotate_peaks(
        [{"mz": theoretical, "height": 100.0}],
        mass_axis=legacy_axis,
        assign_all_library=True,
    )

    assert peaks[0]["candidates"]
    assert peaks[0]["formula_tolerance"]["status"] == "fallback"
    assert peaks[0]["suggested_label"].startswith("unknown m/z")
    assert "suggested_formula" not in peaks[0]
