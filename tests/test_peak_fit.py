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

    def test_global_group_selects_regularisation_and_beats_reduced_model(self):
        rng = np.random.default_rng(4)
        x = np.arange(220)
        centres = np.array([100.0, 107.0])
        sigmas = np.array([3.0, 3.0])
        profile_x = np.linspace(-4.5, 4.5, 181)
        profile = {"x": profile_x, "y": np.exp(-0.5 * profile_x**2)}
        components = np.column_stack(
            [np.exp(-0.5 * ((x - centre) / 3.0) ** 2) for centre in centres]
        )
        cycles = np.arange(240, dtype=np.float64)
        traces = np.column_stack(
            (
                70.0 + 25.0 * np.sin(cycles / 19.0),
                45.0 + 18.0 * np.cos(cycles / 23.0),
            )
        )
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            spectra = 4.0 + traces @ components.T
        spectra += rng.normal(0.0, 0.15, size=spectra.shape)
        fit = peak_fit.fit_group_design(spectra.mean(axis=0), centres, sigmas, profile)

        global_fit = peak_fit.fit_global_group(
            spectra[:, fit["lo"] : fit["hi"]],
            fit,
            breaks=(120,),
            shape={
                "x": x[fit["lo"] : fit["hi"]],
                "centres": centres,
                "sigmas": sigmas,
                "profile_x": profile_x,
                "profile_y": profile["y"],
            },
        )

        self.assertEqual(global_fit["status"], "reliable")
        self.assertGreater(global_fit["selected_lambda"], 0.0)
        self.assertLess(global_fit["held_out_relative_rmse"], 0.05)
        self.assertGreater(global_fit["held_out_improvement"], 0.05)
        recovered = global_fit["amplitudes"]
        np.testing.assert_allclose(recovered, traces, rtol=0.04, atol=1.0)

    def test_global_group_polls_cancellation_during_validation(self):
        x = np.arange(220)
        centres = np.array([100.0, 107.0])
        sigmas = np.array([3.0, 3.0])
        profile_x = np.linspace(-4.5, 4.5, 181)
        profile = {"x": profile_x, "y": np.exp(-0.5 * profile_x**2)}
        spectrum = (
            3.0
            + 80.0 * np.exp(-0.5 * ((x - centres[0]) / 3.0) ** 2)
            + 40.0 * np.exp(-0.5 * ((x - centres[1]) / 3.0) ** 2)
        )
        fit = peak_fit.fit_group_design(spectrum, centres, sigmas, profile)
        calls = 0

        def stop():
            nonlocal calls
            calls += 1
            return calls > 3

        with self.assertRaises(peak_fit.PeakFitCancelled):
            peak_fit.fit_global_group(
                np.tile(spectrum[fit["lo"] : fit["hi"]], (80, 1)),
                fit,
                should_stop=stop,
            )

    def test_global_group_withholds_inactive_component(self):
        x = np.arange(220)
        centres = np.array([100.0, 107.0])
        sigmas = np.array([3.0, 3.0])
        profile_x = np.linspace(-4.5, 4.5, 181)
        profile = {"x": profile_x, "y": np.exp(-0.5 * profile_x**2)}
        first = np.exp(-0.5 * ((x - centres[0]) / 3.0) ** 2)
        spectra = np.tile(3.0 + 80.0 * first, (120, 1))
        fit = peak_fit.fit_group_design(spectra.mean(axis=0), centres, sigmas, profile)
        if fit["status"] != "reliable":
            self.skipTest("average-spectrum fit conservatively withheld the group")

        global_fit = peak_fit.fit_global_group(spectra[:, fit["lo"] : fit["hi"]], fit)

        self.assertEqual(global_fit["status"], "unresolved")
        self.assertTrue(np.isnan(global_fit["amplitudes"]).all())
        self.assertTrue(global_fit["failed_gates"])

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
