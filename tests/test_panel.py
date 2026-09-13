"""Tests for conservative cross-file peak-table adaptation."""

import unittest

from sniff import panel


class PanelAdaptationTest(unittest.TestCase):
    def test_identity_moves_to_target_mass_and_new_peaks_are_retained(self):
        template = [
            {
                "mz": 59.0491,
                "label": "acetone",
                "formula": "C3H6O",
                "k": 3.0,
                "window": {"left": 0.02, "right": 0.03},
                "samples": ["sample_01"],
            }
        ]
        detected = [{"mz": 59.052}, {"mz": 73.064}]

        peaks, diagnostics = panel.adapt_peak_table(template, detected)

        self.assertEqual([peak["mz"] for peak in peaks], [59.052, 73.064])
        self.assertEqual(peaks[0]["formula"], "C3H6O")
        self.assertEqual(peaks[0]["k"], 3.0)
        self.assertNotIn("window", peaks[0])
        self.assertNotIn("samples", peaks[0])
        self.assertEqual(diagnostics["n_matched"], 1)
        self.assertEqual(diagnostics["n_new"], 1)

    def test_transferred_label_does_not_keep_an_automatic_formula(self):
        peaks, _ = panel.adapt_peak_table(
            [{"mz": 59.049, "label": "reviewed unknown"}],
            [{"mz": 59.05, "label": "acetone", "formula": "C3H6O"}],
        )

        self.assertEqual(peaks[0]["label"], "reviewed unknown")
        self.assertNotIn("formula", peaks[0])

    def test_new_automatic_identity_cannot_duplicate_a_transferred_one(self):
        peaks, diagnostics = panel.adapt_peak_table(
            [{"mz": 59.049, "label": "acetone", "formula": "C3H6O"}],
            [
                {"mz": 59.05, "label": "acetone", "formula": "C3H6O"},
                {"mz": 59.07, "label": "acetone", "formula": "C3H6O"},
            ],
        )

        self.assertEqual(peaks[0]["label"], "acetone")
        self.assertEqual(peaks[0]["formula"], "C3H6O")
        self.assertEqual(peaks[1]["label"], "")
        self.assertNotIn("formula", peaks[1])
        self.assertEqual(diagnostics["n_automatic_duplicates_suppressed"], 1)

    def test_missing_target_is_not_snapped_to_an_unrelated_peak(self):
        peaks, diagnostics = panel.adapt_peak_table(
            [{"mz": 59.0, "formula": "C3H6O"}], [{"mz": 60.0}]
        )

        self.assertEqual(peaks, [{"mz": 60.0}])
        self.assertEqual(diagnostics["n_missing"], 1)

    def test_resolution_width_does_not_authorise_identity_transfer(self):
        peaks, diagnostics = panel.adapt_peak_table(
            [{"mz": 100.0, "formula": "C5H8O2"}],
            [{"mz": 100.03}],
            R_phys=1200.0,
        )

        self.assertNotIn("formula", peaks[0])
        self.assertEqual(diagnostics["n_missing"], 1)

    def test_ambiguous_match_does_not_transfer_identity(self):
        peaks, diagnostics = panel.adapt_peak_table(
            [{"mz": 100.0, "formula": "C5H8O2"}],
            [{"mz": 99.997}, {"mz": 100.003}],
            R_phys=2400.0,
        )

        self.assertTrue(all("formula" not in peak for peak in peaks))
        self.assertEqual(diagnostics["n_missing"], 1)
        self.assertIn("ambiguous", diagnostics["missing"][0]["reason"])


if __name__ == "__main__":
    unittest.main()
