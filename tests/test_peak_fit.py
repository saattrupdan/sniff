"""Tests for measured-shape fitting of overlapping peak groups."""

import unittest
import warnings

import numpy as np

from sniff import peak_fit


class PeakFitTest(unittest.TestCase):
    def test_empirical_profile_requires_clean_reference_peaks(self):
        spectrum = np.zeros(400)
        centres = [80.0, 180.0]
        sigmas = [3.0, 3.0]

        profile = peak_fit.estimate_empirical_profile(spectrum, centres, sigmas)

        self.assertFalse(profile["usable"])
        self.assertIn("three", profile["reason"])

    def test_empirical_profile_is_learned_from_three_peaks(self):
        x = np.arange(500)
        centres = [80.0, 220.0, 380.0]
        sigmas = [3.0, 3.0, 3.0]
        spectrum = np.full(500, 2.0)
        for centre in centres:
            profile_x = (x - centre) / 3.0
            shape = np.exp(-0.5 * profile_x**2)
            shape[profile_x > 0] *= np.exp(-0.15 * profile_x[profile_x > 0])
            spectrum += 100.0 * shape

        profile = peak_fit.estimate_empirical_profile(spectrum, centres, sigmas)

        self.assertTrue(profile["usable"])
        self.assertEqual(profile["n_reference_peaks"], 3)
        self.assertAlmostEqual(float(np.max(profile["y"])), 1.0)

    def test_group_fit_recovers_nonnegative_amplitudes_and_baseline(self):
        x = np.arange(220)
        centres = np.array([100.0, 107.0])
        sigmas = np.array([3.0, 3.0])
        profile = {
            "x": np.linspace(-4.5, 4.5, 181),
            "y": np.exp(-0.5 * np.linspace(-4.5, 4.5, 181) ** 2),
        }
        first = 80.0 * np.exp(-0.5 * ((x - centres[0]) / 3.0) ** 2)
        second = 35.0 * np.exp(-0.5 * ((x - centres[1]) / 3.0) ** 2)
        spectrum = 4.0 + first + second

        fit = peak_fit.fit_group_design(spectrum, centres, sigmas, profile)
        block = np.vstack((spectrum, 4.0 + 0.5 * first + 2.0 * second))
        amplitudes = peak_fit.apply_group_design(block, fit)

        self.assertEqual(fit["status"], "reliable")
        np.testing.assert_allclose(amplitudes[0], [80.0, 35.0], rtol=0.04)
        np.testing.assert_allclose(amplitudes[1], [40.0, 70.0], rtol=0.04)
        self.assertTrue((amplitudes >= 0).all())

    def test_extreme_spectrum_is_withheld_without_numerical_warnings(self):
        centres = np.array([100.0, 100.02, 100.04])
        sigmas = np.full(3, 3.0)
        profile_x = np.linspace(-4.5, 4.5, 181)
        profile = {"x": profile_x, "y": np.exp(-0.5 * profile_x**2)}
        spectrum = np.full(220, np.finfo(np.float64).max / 4.0)

        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            fit = peak_fit.fit_group_design(spectrum, centres, sigmas, profile)

        self.assertIn(fit["status"], {"unresolved", "failed"})

    def test_ill_conditioned_group_is_unresolved(self):
        x = np.arange(220)
        centres = np.array([100.0, 100.01])
        sigmas = np.array([3.0, 3.0])
        profile_x = np.linspace(-4.5, 4.5, 181)
        profile = {"x": profile_x, "y": np.exp(-0.5 * profile_x**2)}
        spectrum = 50.0 * np.exp(-0.5 * ((x - 100.0) / 3.0) ** 2)

        fit = peak_fit.fit_group_design(spectrum, centres, sigmas, profile)

        self.assertEqual(fit["status"], "unresolved")
        self.assertGreater(fit["condition"], 100.0)


if __name__ == "__main__":
    unittest.main()
