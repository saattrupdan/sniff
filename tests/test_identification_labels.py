#!/usr/bin/env python3
"""Compound naming and sample-specific inclusion carried into the review app.

Two review contracts live here:

* a compound is named once, so an auto-generated ``unknown m/z ...`` label must
  never sit next to an assigned formula (it says "unknown" and identifies the
  compound at the same time);
* ``peaks[].samples`` is carried through to the browser unchanged, while the
  delivered summary rows stay exactly as they were.
"""

from __future__ import annotations

import unittest
from unittest import mock

import h5py
import numpy as np
from calibration_helpers import identity_mass_axis

from sniff import analyze, formula_id, ptrms, viz


def _candidate(formula: str, name: str) -> dict:
    return {
        "formula": formula,
        "name": name,
        "ion_mz": 59.049,
        "delta_mDa": 1.0,
        "dbe": 1.0,
        "k": None,
        "k_estimated": True,
        "flags": [],
        "iso_pred": [0.01, 0.01],
        "iso_obs": None,
        "iso_used": False,
        "probability": 0.99,
    }


class IdentityLabelTest(unittest.TestCase):
    def test_formula_replaces_an_auto_generated_unknown_label(self):
        self.assertEqual(
            formula_id.identity_label("unknown m/z 73.029", "C4H8O"), "C4H8O"
        )

    def test_expert_wording_is_never_treated_as_a_placeholder(self):
        # "Unknown terpenes" is a statement about the sample, not a blank to fill
        self.assertEqual(
            formula_id.identity_label("Unknown terpenes", "C5H8"), "Unknown terpenes"
        )

    def test_an_unassigned_candidate_does_not_rename_the_peak(self):
        # the Identification card offers candidates; naming the peak from one would
        # turn an unreviewed guess into the compound's identity
        self.assertEqual(
            formula_id.identity_label("unknown m/z 59.049", ""), "unknown m/z 59.049"
        )

    def test_hand_written_label_wins(self):
        self.assertEqual(
            formula_id.identity_label("2-butanone", "C4H8O"), "2-butanone"
        )

    def test_unknown_label_without_a_formula_is_left_alone(self):
        self.assertEqual(
            formula_id.identity_label("unknown m/z 131.104", ""), "unknown m/z 131.104"
        )

    def test_missing_label_falls_back_to_the_formula(self):
        self.assertEqual(formula_id.identity_label("", "C6H6"), "C6H6")
        # nothing known at all stays blank in the CSV, as it always was
        self.assertEqual(formula_id.identity_label(None, None), "")


class ContextualPriorTest(unittest.TestCase):
    def test_a_compound_of_interest_modestly_boosts_its_formula(self):
        neutral_mass = 50.0
        candidates = [
            ({"C": 3, "H": 6, "O": 1}, neutral_mass),
            ({"C": 2, "H": 6}, neutral_mass),
        ]
        with (
            mock.patch.object(
                formula_id, "enumerate_formulas", return_value=candidates
            ),
            mock.patch.object(formula_id, "_prior", return_value=1.0),
            mock.patch.object(formula_id, "_known", return_value=None),
        ):
            scored = formula_id.score_peak(
                neutral_mass + formula_id.PROTON,
                1.0,
                compounds_of_interest=[
                    {"name": "acetone", "formula": "forged formula ignored"}
                ],
            )

        self.assertEqual(scored[0]["formula"], "C3H6O")
        self.assertEqual(scored[0]["interest_matches"], ["acetone"])
        self.assertEqual(scored[0]["probability"], 0.667)
        self.assertEqual(scored[1]["probability"], 0.333)

    def test_selected_isomer_becomes_the_preferred_editable_name(self):
        neutral_mass = 50.0
        with (
            mock.patch.object(
                formula_id,
                "enumerate_formulas",
                return_value=[({"C": 3, "H": 6, "O": 1}, neutral_mass)],
            ),
            mock.patch.object(formula_id, "_prior", return_value=1.0),
        ):
            candidate = formula_id.score_peak(
                neutral_mass + formula_id.PROTON,
                1.0,
                compounds_of_interest=["propanal"],
            )[0]

        self.assertEqual(candidate["name"], "acetone")
        self.assertEqual(candidate["preferred_name"], "propanal")
        self.assertIn("acetone", candidate["names"])
        self.assertIn("propanal", candidate["names"])


def _payload(peaks_cfg, candidates=None):
    """Build a review payload for `peaks_cfg` from a tiny synthetic file."""
    with h5py.File("in-memory", "w", driver="core", backing_store=False) as h5:
        h5.create_dataset("SPECdata/Intensities", data=np.zeros((2, 5)))
        h5.create_dataset("SPECdata/AverageSpec", data=np.ones(5))
        with (
            mock.patch.object(
                ptrms, "load_mass_axis", return_value=identity_mass_axis()
            ),
            mock.patch.object(
                ptrms,
                "load_transmission",
                return_value=(np.array([1.0, 1000.0]), np.array([1.0, 1.0])),
            ),
            mock.patch.object(ptrms, "spec_duration_s", return_value=1.0),
            mock.patch.object(ptrms, "extract_primary", return_value=None),
            mock.patch.object(ptrms, "water_cluster_ratio", return_value=None),
            mock.patch.object(ptrms, "build_discriminator", return_value=np.ones(2)),
            mock.patch.object(
                ptrms, "derive_molar_volume_info", return_value=(24.465, "test")
            ),
            mock.patch.object(ptrms, "derive_K", return_value=None),
            mock.patch.object(
                ptrms,
                "extract_traces",
                return_value=({2.0: (np.array([2.0, 4.0]), 2.0)}, (10.0, 1.0)),
            ),
            mock.patch.object(ptrms, "_cluster", return_value=[]),
            mock.patch.object(ptrms, "resolve_k", return_value={}),
            mock.patch.object(
                ptrms, "load_rate_constants", return_value={"compounds": []}
            ),
            mock.patch.object(
                viz.formula_id, "score_peak", return_value=candidates or []
            ),
        ):
            return viz.build_viz_data(
                h5,
                peaks_cfg=peaks_cfg,
                ranges_cfg=[
                    {"label": "sample_01", "start": 1, "end": 1, "class": "sample"},
                    {"label": "sample_02", "start": 2, "end": 2, "class": "sample"},
                ],
            )


class AutomaticAssignmentTest(unittest.TestCase):
    def test_library_candidates_fill_ambiguous_peaks_with_formulas(self):
        peaks = [
            {
                "mz": 59.049,
                "height": 20.0,
                "id_confidence": 0.45,
                "id_ambiguous": True,
                "candidates": [
                    {
                        "formula": "C3H6O",
                        "name": "acetone",
                        "preferred_name": "propanal",
                        "probability": 0.45,
                        "interest_matches": ["propanal"],
                    }
                ],
            }
        ]

        analyze._assign_suggested_identities(peaks, assign_all_library=True)

        self.assertEqual(peaks[0]["suggested_label"], "propanal")
        self.assertEqual(peaks[0]["suggested_formula"], "C3H6O")
        self.assertEqual(peaks[0]["suggested_candidate_rank"], 1)

    def test_global_matching_uses_distinct_fallback_compounds(self):
        peaks = [
            {
                "mz": 50.0,
                "height": 10.0,
                "candidates": [
                    {
                        "formula": "C3H6O",
                        "name": "acetone",
                        "probability": 0.7,
                    },
                    {
                        "formula": "C4H8",
                        "name": "butene",
                        "probability": 0.3,
                    },
                ],
            },
            {
                "mz": 51.0,
                "height": 20.0,
                "candidates": [
                    {
                        "formula": "C3H6O",
                        "name": "acetone",
                        "probability": 1.0,
                    }
                ],
            },
        ]

        analyze._assign_suggested_identities(peaks, assign_all_library=True)

        self.assertEqual(
            [(peak["suggested_label"], peak["suggested_formula"]) for peak in peaks],
            [("butene", "C4H8"), ("acetone", "C3H6O")],
        )
        self.assertEqual(peaks[0]["suggested_candidate_rank"], 2)

    def test_global_matching_maximises_quality_after_coverage(self):
        peaks = [
            {
                "mz": 50.0,
                "height": 20.0,
                "candidates": [
                    {"formula": "C3H6O", "name": "acetone", "probability": 0.99},
                    {"formula": "C4H8", "name": "butene", "probability": 0.01},
                ],
            },
            {
                "mz": 51.0,
                "height": 10.0,
                "candidates": [
                    {"formula": "C3H6O", "name": "acetone", "probability": 0.34},
                    {"formula": "C4H8", "name": "butene", "probability": 0.33},
                    {"formula": "C5H10", "name": "pentene", "probability": 0.33},
                ],
            },
        ]

        analyze._assign_suggested_identities(peaks, assign_all_library=True)

        self.assertEqual(peaks[0]["suggested_formula"], "C3H6O")
        self.assertEqual(peaks[1]["suggested_formula"], "C4H8")

    def test_duplicate_without_fallback_is_left_unassigned(self):
        candidate = {
            "formula": "C3H6O",
            "name": "acetone",
            "probability": 1.0,
        }
        peaks = [
            {"mz": 59.048, "height": 5.0, "candidates": [candidate]},
            {"mz": 59.049, "height": 20.0, "candidates": [candidate]},
        ]

        analyze._assign_suggested_identities(peaks, assign_all_library=True)

        self.assertEqual(peaks[1]["suggested_label"], "acetone")
        self.assertTrue(peaks[0]["suggested_label"].startswith("unknown m/z"))
        self.assertNotIn("suggested_formula", peaks[0])

    def test_named_duplicate_uses_the_next_formula_only_guess(self):
        candidates = [
            {
                "formula": "C3H6O",
                "name": "acetone",
                "probability": 0.95,
            },
            {"formula": "C3H8", "name": None, "probability": 0.05},
        ]
        peaks = [
            {
                "mz": 59.048,
                "height": 5.0,
                "id_confidence": 0.95,
                "candidates": candidates,
            },
            {
                "mz": 59.049,
                "height": 20.0,
                "id_confidence": 0.95,
                "candidates": candidates,
            },
        ]

        analyze._assign_suggested_identities(peaks, assign_all_library=True)

        self.assertEqual(peaks[1]["suggested_formula"], "C3H6O")
        self.assertEqual(peaks[0]["suggested_label"], "C3H8")
        self.assertEqual(peaks[0]["suggested_formula"], "C3H8")
        self.assertEqual(peaks[0]["suggested_candidate_rank"], 2)

    def test_formula_only_guesses_are_globally_unique(self):
        candidates = [
            {"formula": "C3H8", "name": None, "probability": 0.95},
            {"formula": "C2H6O", "name": None, "probability": 0.05},
        ]
        peaks = [
            {
                "mz": 45.0,
                "height": 5.0,
                "id_confidence": 0.95,
                "candidates": candidates,
            },
            {
                "mz": 45.001,
                "height": 20.0,
                "id_confidence": 0.95,
                "candidates": candidates,
            },
        ]

        analyze._assign_suggested_identities(peaks, assign_all_library=True)

        self.assertEqual(
            {peak["suggested_formula"] for peak in peaks}, {"C3H8", "C2H6O"}
        )
        self.assertTrue(all(peak["suggested_label"] for peak in peaks))

    def test_ambiguous_formula_only_peak_gets_an_app_review_default(self):
        peaks = [
            {
                "mz": 45.0,
                "height": 5.0,
                "id_confidence": 0.4,
                "id_ambiguous": True,
                "candidates": [
                    {"formula": "C3H8", "name": None, "probability": 0.4},
                    {"formula": "C2H6O", "name": None, "probability": 0.35},
                ],
            }
        ]

        analyze._assign_suggested_identities(peaks, assign_all_library=True)

        self.assertEqual(peaks[0]["suggested_label"], "C3H8")
        self.assertEqual(peaks[0]["suggested_formula"], "C3H8")
        self.assertEqual(peaks[0]["suggested_candidate_rank"], 1)

    def test_headless_default_keeps_a_low_share_library_match_unassigned(self):
        peaks = [
            {
                "mz": 59.049,
                "id_confidence": 0.45,
                "id_ambiguous": True,
                "candidates": [
                    {
                        "formula": "C3H6O",
                        "name": "acetone",
                        "probability": 0.45,
                    }
                ],
            }
        ]

        analyze._assign_suggested_identities(peaks)

        self.assertTrue(peaks[0]["suggested_label"].startswith("unknown m/z"))
        self.assertNotIn("suggested_formula", peaks[0])


class ReviewPayloadTest(unittest.TestCase):
    def test_assigned_formula_wins_over_unknown_label_in_payload(self):
        data = _payload(
            [{"mz": 2.0, "label": "unknown m/z 2.000", "formula": "C3H6O"}],
            candidates=[_candidate("C3H6O", "acetone")],
        )
        peak = data["peaks"][0]
        self.assertEqual(peak["label"], "C3H6O")
        self.assertNotIn("unknown", peak["label"].lower())
        # the name is still offered for review, it is just not the compound's name
        self.assertEqual(peak["candidates"][0]["name"], "acetone")

    def test_unknown_label_stands_alone_without_a_formula(self):
        data = _payload([{"mz": 2.0, "label": "unknown m/z 2.000"}])
        self.assertEqual(data["peaks"][0]["label"], "unknown m/z 2.000")
        self.assertEqual(data["peaks"][0]["formula"], "")

    def test_sample_selection_is_carried_into_the_payload(self):
        data = _payload(
            [
                {"mz": 2.0, "label": "acetone", "samples": ["sample_02"]},
                {"mz": 2.0, "label": "methanol"},
            ]
        )
        self.assertEqual(data["peaks"][0]["samples"], ["sample_02"])
        # no explicit list means every sample interval, as in older configs
        self.assertIsNone(data["peaks"][1]["samples"])


if __name__ == "__main__":
    unittest.main()


def test_sample_selection_is_inert_in_the_summary_analysis():
    """`samples` says which sample intervals a compound is part of.

    Until per-sample output exists it must stay metadata: a compound selected for
    at least one sample is summarised exactly as it was before the field existed.
    """
    from sniff import analyze

    plain = {"mz": 100.0, "label": "analyte", "window": 0.4, "use": True}
    tagged = dict(plain, samples=["sample_01"])

    assert analyze._peak_windows([plain]) == analyze._peak_windows([tagged])
    constants = ptrms.load_rate_constants()
    assert ptrms.resolve_k([plain], constants) == ptrms.resolve_k([tagged], constants)
