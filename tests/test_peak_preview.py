"""Tests for authoritative hand-drawn peak previews."""

import unittest

import h5py
import numpy as np
from calibration_helpers import identity_mass_axis

from sniff import viz


class PeakPreviewTest(unittest.TestCase):
    def test_drawn_region_snaps_to_apex_and_returns_formula_candidates(self):
        a = 10000.0
        axis = identity_mass_axis(a=a, b=0.0)
        expected = 59.0491
        centre = a * np.sqrt(expected)
        average = np.ones(80000, dtype=np.float64)
        bins = np.arange(len(average))
        average += 1000.0 * np.exp(-0.5 * ((bins - centre) / 7.0) ** 2)

        with h5py.File("in-memory", "w", driver="core", backing_store=False) as h5:
            h5.create_dataset("SPECdata/AverageSpec", data=average)
            preview = viz.preview_peak(
                h5,
                59.0,
                59.1,
                mass_axis=axis,
            )

        self.assertAlmostEqual(preview["apex"], expected, delta=0.001)
        self.assertTrue(preview["candidates"])
        self.assertIn(
            "C3H6O",
            {candidate["formula"] for candidate in preview["candidates"]},
        )
        self.assertIn("isotope_model", preview["candidates"][0])


if __name__ == "__main__":
    unittest.main()
