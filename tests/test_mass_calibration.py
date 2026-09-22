"""Regression tests for multi-point mass calibration loading."""

import unittest

import h5py
import numpy as np

from sniff import analyze, ptrms

# In-memory reproduction of the three-anchor calibration used by Data_10_26_33.
# Keeping these values here avoids depending on the external measurement file.
DATA_10_26_33_MAPPING = np.array(
    [
        [21.022100, 24299.236],
        [203.94299, 130446.04],
        [330.84799, 173231.97],
    ],
    dtype=np.float32,
)


class MassCalibrationTest(unittest.TestCase):
    @staticmethod
    def _file(mapping=None, spectrum=None):
        h5 = h5py.File("in-memory", "w", driver="core", backing_store=False)
        if mapping is not None:
            h5.create_dataset("CALdata/Mapping", data=mapping)
        if spectrum is not None:
            h5.create_dataset("CALdata/Spectrum", data=spectrum)
        return h5

    def test_three_mapping_anchors_fit_with_low_residual(self):
        with self._file(mapping=DATA_10_26_33_MAPPING) as h5:
            a, b = ptrms.load_mass_cal(h5)

        masses = DATA_10_26_33_MAPPING[:, 0]
        timebins = DATA_10_26_33_MAPPING[:, 1]
        inferred_masses = ((timebins - b) / a) ** 2
        residual_ppm = np.abs((inferred_masses - masses) / masses) * 1e6

        self.assertLessEqual(float(residual_ppm.max()), 10.0)
        self.assertGreater(a, 0.0)
        self.assertTrue(np.isfinite([a, b]).all())

    def test_two_mapping_anchors_keep_closed_form_calibration(self):
        mapping = np.array([[19.0, 500.0], [181.0, 1500.0]])
        expected_a = (1500.0 - 500.0) / (np.sqrt(181.0) - np.sqrt(19.0))
        expected_b = 500.0 - expected_a * np.sqrt(19.0)

        with self._file(mapping=mapping) as h5:
            actual_a, actual_b = ptrms.load_mass_cal(h5)

        self.assertEqual(actual_a, expected_a)
        self.assertEqual(actual_b, expected_b)

    def test_float32_two_mapping_anchors_match_legacy_expression_exactly(self):
        mapping = np.array([[19.0, 500.0], [181.0, 1500.0]], dtype=np.float32)
        expected_a = (mapping[1, 1] - mapping[0, 1]) / (
            np.sqrt(mapping[1, 0]) - np.sqrt(mapping[0, 0])
        )
        expected_b = mapping[0, 1] - expected_a * np.sqrt(mapping[0, 0])

        with self._file(mapping=mapping) as h5:
            actual_a, actual_b = ptrms.load_mass_cal(h5)

        self.assertEqual(actual_a, float(expected_a))
        self.assertEqual(actual_b, float(expected_b))

    def test_valid_mapping_is_preferred_over_spectrum_fallback(self):
        spectrum = np.array([[900.0, 1.0], [901.0, 1.0]])
        with self._file(mapping=DATA_10_26_33_MAPPING, spectrum=spectrum) as h5:
            a, b = ptrms.load_mass_cal(h5)

        mapping = DATA_10_26_33_MAPPING.astype(np.float64)
        design = np.column_stack((np.sqrt(mapping[:, 0]), np.ones(3)))
        expected_a, expected_b = np.linalg.lstsq(design, mapping[:, 1], rcond=None)[0]
        self.assertAlmostEqual(a, expected_a, places=10)
        self.assertAlmostEqual(b, expected_b, places=10)
        self.assertNotEqual((a, b), (900.5, 1.0))

    def test_malformed_mapping_falls_back_to_spectrum(self):
        mapping = np.ones((3, 3))
        spectrum = np.array([[10.0, 2.0], [12.0, 4.0]])

        with self._file(mapping=mapping, spectrum=spectrum) as h5:
            self.assertEqual(ptrms.load_mass_cal(h5), (11.0, 3.0))

    def test_malformed_two_mapping_anchors_fall_back_to_spectrum(self):
        mapping = np.array([[19.0, np.nan], [181.0, 1500.0]])
        spectrum = np.array([[10.0, 2.0], [12.0, 4.0]])

        with self._file(mapping=mapping, spectrum=spectrum) as h5:
            self.assertEqual(ptrms.load_mass_cal(h5), (11.0, 3.0))

    def test_non_monotonic_multi_mapping_falls_back_to_spectrum(self):
        mapping = np.array([[19.0, 500.0], [59.0, 700.0], [181.0, 650.0]])
        spectrum = np.array([[10.0, 2.0], [12.0, 4.0]])

        with self._file(mapping=mapping, spectrum=spectrum) as h5:
            self.assertEqual(ptrms.load_mass_cal(h5), (11.0, 3.0))

    def test_grossly_nonlinear_mapping_falls_back_to_spectrum(self):
        # The anchors are positive, strictly monotonic, and well-conditioned, but
        # do not describe the stated square-root time-of-flight model.
        mapping = np.array([[19.0, 500.0], [59.0, 1000.0], [181.0, 2000.0]])
        spectrum = np.array([[10.0, 2.0], [12.0, 4.0]])

        with self._file(mapping=mapping, spectrum=spectrum) as h5:
            self.assertEqual(ptrms.load_mass_cal(h5), (11.0, 3.0))

    def test_invalid_multi_mapping_row_falls_back_to_spectrum(self):
        mapping = np.array([[19.0, 500.0], [59.0, np.nan], [181.0, 1500.0]])
        spectrum = np.array([[10.0, 2.0], [12.0, 4.0]])

        with self._file(mapping=mapping, spectrum=spectrum) as h5:
            self.assertEqual(ptrms.load_mass_cal(h5), (11.0, 3.0))

    def test_non_positive_multi_mapping_row_falls_back_to_spectrum(self):
        mapping = np.array([[19.0, 500.0], [59.0, -1.0], [181.0, 1500.0]])
        spectrum = np.array([[10.0, 2.0], [12.0, 4.0]])

        with self._file(mapping=mapping, spectrum=spectrum) as h5:
            self.assertEqual(ptrms.load_mass_cal(h5), (11.0, 3.0))

    def test_ill_conditioned_multi_mapping_falls_back_to_spectrum(self):
        mapping = np.array([[100.0, 500.0], [100.00001, 501.0], [100.00002, 502.0]])
        spectrum = np.array([[10.0, 2.0], [12.0, 4.0]])

        with self._file(mapping=mapping, spectrum=spectrum) as h5:
            self.assertEqual(ptrms.load_mass_cal(h5), (11.0, 3.0))

    def test_degenerate_mapping_falls_back_to_spectrum(self):
        mapping = np.array([[19.0, 500.0], [19.0, 600.0], [19.0, 700.0]])
        spectrum = np.array([[10.0, 2.0], [12.0, 4.0]])

        with self._file(mapping=mapping, spectrum=spectrum) as h5:
            self.assertEqual(ptrms.load_mass_cal(h5), (11.0, 3.0))

    def test_unusable_calibration_raises_clear_error(self):
        mapping = np.array([[19.0, 500.0], [19.0, 600.0]])
        spectrum = np.array([[0.0, 2.0], [np.nan, 4.0]])

        with (
            self._file(mapping=mapping, spectrum=spectrum) as h5,
            self.assertRaisesRegex(ValueError, "no mass calibration"),
        ):
            ptrms.load_mass_cal(h5)


class InternalMassAxisCalibrationTest(unittest.TestCase):
    A = 1000.0
    B = 0.0
    SCALE = 1.0007
    OFFSET = -0.020
    NBIN = 15000

    @classmethod
    def _observed_mass(cls, corrected_mass):
        return (corrected_mass - cls.OFFSET) / cls.SCALE

    @classmethod
    def _spectrum(cls, peaks):
        spectrum = np.full(cls.NBIN, 4.0, dtype=np.float64)
        for mass, height in peaks:
            centre = cls.A * np.sqrt(mass) + cls.B
            lo = max(0, int(np.floor(centre)) - 5)
            hi = min(cls.NBIN, int(np.floor(centre)) + 7)
            bins = np.arange(lo, hi, dtype=np.float64)
            # A sampled parabola gives a known non-integer vertex, exercising the
            # sub-bin centre calculation without requiring SciPy fitting.
            shape = np.maximum(0.0, height * (1.0 - ((bins - centre) / 3.5) ** 2))
            spectrum[lo:hi] += shape
        return spectrum

    @classmethod
    def _file(cls, peaks=None, average=None):
        h5 = h5py.File("mass-axis", "w", driver="core", backing_store=False)
        h5.create_dataset("CALdata/Spectrum", data=np.array([[cls.A, cls.B]]))
        if average is None:
            average = cls._spectrum(peaks or [])
        h5.create_dataset("SPECdata/AverageSpec", data=average)
        if np.asarray(average).ndim == 1:
            h5.create_dataset(
                "SPECdata/Intensities",
                data=np.vstack([average, average]),
            )
        return h5

    @classmethod
    def _good_peaks(cls, intermediate=True):
        peaks = [
            (cls._observed_mass(37.028405), 1200.0),
            (cls._observed_mass(203.942993), 1000.0),
        ]
        if intermediate:
            peaks.append((cls._observed_mass(100.123), 800.0))
        return peaks

    def test_multi_point_mapping_is_authoritative_without_affine_shift(self):
        masses = np.array([21.0221, 203.9430, 330.8480], dtype=np.float64)
        mapping = np.column_stack((masses, self.A * np.sqrt(masses) + self.B))
        with self._file(peaks=[]) as h5:
            h5.create_dataset("CALdata/Mapping", data=mapping)
            axis = ptrms.load_mass_axis(h5)

        self.assertEqual(axis.scale, 1.0)
        self.assertEqual(axis.offset, 0.0)
        self.assertEqual(axis.diagnostics["authority"], "CALdata/Mapping")
        self.assertFalse(axis.diagnostics["mass_domain_correction_applied"])
        points = axis.diagnostics["mapping_calibration"]["points"]
        self.assertEqual(len(points), 3)
        self.assertLessEqual(
            max(abs(point["residual_ppm"]) for point in points), 10.0
        )

    def test_degraded_mapping_uses_stable_exact_mass_internal_references(self):
        masses = np.array([21.0221, 203.9430, 330.8480], dtype=np.float64)
        shifted = masses * (1.0 + np.array([0.0, 20.0, 0.0]) * 1e-6)
        mapping = np.column_stack(
            (masses, self.A * np.sqrt(shifted) + self.B)
        )
        with self._file(peaks=self._good_peaks()) as h5:
            h5.create_dataset("CALdata/Mapping", data=mapping)
            axis = ptrms.load_mass_axis(h5)

        self.assertEqual(
            axis.diagnostics["model"], ptrms.INTERNAL_MASS_CORRECTION_MODEL
        )
        self.assertIn("degraded CALdata/Mapping", axis.diagnostics["authority"])
        self.assertEqual(axis.diagnostics["mapping_validation"]["status"], "degraded")
        self.assertTrue(axis.diagnostics["mass_domain_correction_applied"])
        self.assertAlmostEqual(axis.scale, self.SCALE, places=4)
        self.assertAlmostEqual(axis.offset, self.OFFSET, places=3)

    def test_two_anchors_apply_shift_and_scale_to_the_whole_axis(self):
        with self._file(peaks=self._good_peaks()) as h5:
            axis = ptrms.load_mass_axis(h5)
            detected = analyze.detect_peaks(h5, mass_axis=axis)
            traces, _ = ptrms.extract_traces(h5, [100.123], R=1200.0, mass_axis=axis)

        self.assertTrue(axis.applied)
        self.assertNotEqual(axis.scale, 1.0)
        self.assertNotEqual(axis.offset, 0.0)
        self.assertAlmostEqual(axis.scale, self.SCALE, places=9)
        self.assertAlmostEqual(axis.offset, self.OFFSET, places=9)
        for anchor in axis.to_dict()["anchors"]:
            self.assertAlmostEqual(
                anchor["corrected_mz"], anchor["target_mz"], places=12
            )
            self.assertNotEqual(anchor["timebin"], round(anchor["timebin"]))
        intermediate_file_mass = self._observed_mass(100.123)
        self.assertAlmostEqual(
            float(axis.file_to_corrected(intermediate_file_mass)), 100.123, places=9
        )
        self.assertAlmostEqual(
            float(axis.tb_to_m(axis.m_to_tb(100.123))), 100.123, places=12
        )
        self.assertTrue(any(abs(peak["mz"] - 100.123) < 0.012 for peak in detected))
        self.assertGreater(float(traces[100.123][0].mean()), 100.0)
        self.assertLess(abs(traces[100.123][1] - 100.123), 0.015)

    def test_missing_anchors_raise_a_structured_calibration_error(self):
        with self._file(peaks=[]) as h5:
            with self.assertRaises(ptrms.MassCalibrationError) as caught:
                ptrms.load_mass_axis(h5)

        diagnostics = caught.exception.diagnostics
        self.assertFalse(diagnostics["applied"])
        self.assertIn("water_cluster anchor missing", diagnostics["fallback_reason"])
        self.assertIn("iodobenzene anchor missing", diagnostics["fallback_reason"])

    def test_weak_anchor_reports_thresholds_in_a_structured_error(self):
        peaks = self._good_peaks(intermediate=False)
        peaks[0] = (peaks[0][0], 5.0)
        with self._file(peaks=peaks) as h5:
            with self.assertRaises(ptrms.MassCalibrationError) as caught:
                ptrms.load_mass_axis(h5)

        diagnostics = caught.exception.diagnostics
        self.assertEqual(diagnostics["anchors"][0]["status"], "weak")
        self.assertIn("requires at least", diagnostics["fallback_reason"])

    def test_resolved_competing_anchor_is_ambiguous(self):
        peaks = self._good_peaks(intermediate=False)
        water = self._observed_mass(37.028405)
        peaks.extend([(water - 0.08, 900.0), (water + 0.08, 850.0)])
        # Remove the central water peak so two separate candidates compete.
        peaks = peaks[1:]
        with self._file(peaks=peaks) as h5:
            with self.assertRaises(ptrms.MassCalibrationError) as caught:
                ptrms.load_mass_axis(h5)

        diagnostics = caught.exception.diagnostics
        self.assertEqual(diagnostics["anchors"][0]["status"], "ambiguous")
        self.assertIn("multiple resolved maxima", diagnostics["fallback_reason"])

    def test_malformed_spectrum_has_precise_fallback(self):
        average = np.full((2, 20), np.nan)
        with self._file(average=average) as h5:
            with self.assertRaises(ptrms.MassCalibrationError) as caught:
                ptrms.load_mass_axis(h5)

        self.assertIn(
            "expected a one-dimensional array",
            caught.exception.diagnostics["fallback_reason"],
        )

    def test_nonfinite_spectrum_has_precise_fallback(self):
        with self._file(average=np.full(self.NBIN, np.nan)) as h5:
            with self.assertRaises(ptrms.MassCalibrationError) as caught:
                ptrms.load_mass_axis(h5)

        self.assertIn("no finite bins", caught.exception.diagnostics["fallback_reason"])

    def test_missing_raw_cycles_never_enable_calibration(self):
        with self._file(peaks=self._good_peaks()) as h5:
            del h5["SPECdata/Intensities"]
            with self.assertRaises(ptrms.MassCalibrationError) as caught:
                ptrms.load_mass_axis(h5)

        self.assertFalse(caught.exception.diagnostics["applied"])
        for anchor in caught.exception.diagnostics["anchors"]:
            self.assertFalse(anchor["persistence"]["available"])

    def test_malformed_raw_cycles_are_structured(self):
        with self._file(peaks=self._good_peaks()) as h5:
            del h5["SPECdata/Intensities"]
            h5.create_dataset("SPECdata/Intensities", data=np.ones(self.NBIN))
            with self.assertRaises(ptrms.MassCalibrationError) as caught:
                ptrms.load_mass_axis(h5)

        self.assertIn("malformed", caught.exception.diagnostics["fallback_reason"])

    def test_fewer_than_two_usable_raw_cycles_are_rejected(self):
        with self._file(peaks=self._good_peaks()) as h5:
            raw = h5["SPECdata/Intensities"]
            raw[1:, :] = 0.0
            with self.assertRaises(ptrms.MassCalibrationError) as caught:
                ptrms.load_mass_axis(h5)

        self.assertIn(
            "fewer than two usable", caught.exception.diagnostics["fallback_reason"]
        )

    def test_unrelated_nonfinite_raw_bins_do_not_reject_calibration(self):
        with self._file(peaks=self._good_peaks()) as h5:
            h5["SPECdata/Intensities"][0, 0] = np.nan
            axis = ptrms.load_mass_axis(h5)

        self.assertTrue(axis.applied)
        for anchor in axis.diagnostics["anchors"]:
            self.assertTrue(anchor["persistence"]["available"])

    def test_nonfinite_anchor_window_is_rejected_with_anchor_diagnostic(self):
        water = int(self.A * np.sqrt(self._observed_mass(37.028405)))
        with self._file(peaks=self._good_peaks()) as h5:
            h5["SPECdata/Intensities"][:, water - 30 : water + 31] = np.nan
            with self.assertRaises(ptrms.MassCalibrationError) as caught:
                ptrms.load_mass_axis(h5)

        diagnostics = caught.exception.diagnostics
        self.assertIn("water_cluster anchor", diagnostics["fallback_reason"])
        self.assertIn(
            "fewer than two usable raw cycles", diagnostics["fallback_reason"]
        )
        self.assertFalse(diagnostics["anchors"][0]["persistence"]["available"])

    def test_unusable_anchor_window_is_rejected_with_anchor_diagnostic(self):
        water = int(self.A * np.sqrt(self._observed_mass(37.028405)))
        with self._file(peaks=self._good_peaks()) as h5:
            h5["SPECdata/Intensities"][:, water - 30 : water + 31] = 0.0
            with self.assertRaises(ptrms.MassCalibrationError) as caught:
                ptrms.load_mass_axis(h5)

        diagnostics = caught.exception.diagnostics
        self.assertIn("water_cluster anchor", diagnostics["fallback_reason"])
        self.assertIn(
            "fewer than two usable raw cycles", diagnostics["fallback_reason"]
        )
        self.assertFalse(diagnostics["anchors"][0]["persistence"]["available"])

    def test_supplied_unapplied_axes_are_rejected_by_production_entry_points(self):
        axis = ptrms.MassAxisCalibration(
            self.A, self.B, diagnostics={"applied": False, "fallback_reason": "test"}
        )
        with self._file(peaks=self._good_peaks()) as h5:
            calls = (
                lambda: analyze.detect_peaks(h5, mass_axis=axis),
                lambda: ptrms.extract_primary(h5, mass_axis=axis),
                lambda: ptrms.water_cluster_ratio(h5, mass_axis=axis),
                lambda: ptrms.build_discriminator(h5, mass_axis=axis),
                lambda: ptrms.extract_traces(h5, [100.0], mass_axis=axis),
            )
            for call in calls:
                with self.subTest(entry_point=call):
                    with self.assertRaises(ptrms.MassCalibrationError):
                        call()

    def test_one_anchor_never_enables_partial_correction(self):
        water = (self._observed_mass(37.028405), 1200.0)
        with self._file(peaks=[water]) as h5:
            with self.assertRaises(ptrms.MassCalibrationError) as caught:
                ptrms.load_mass_axis(h5)

        diagnostics = caught.exception.diagnostics
        self.assertEqual(diagnostics["anchors"][0]["status"], "accepted")
        self.assertEqual(diagnostics["anchors"][1]["status"], "missing")
        self.assertIn("iodobenzene anchor missing", diagnostics["fallback_reason"])

    def test_supplied_axis_validation_rejects_contradictory_evidence(self):
        import copy

        with self._file(peaks=self._good_peaks()) as h5:
            axis = ptrms.load_mass_axis(h5)
        for field, value in (
            ("scale", axis.scale + 0.01),
            ("offset_da", axis.offset + 0.01),
        ):
            diagnostics = copy.deepcopy(axis.diagnostics)
            diagnostics[field] = value
            with self.subTest(field=field):
                with self.assertRaises(ptrms.MassCalibrationError):
                    ptrms.validate_mass_axis(
                        ptrms.MassAxisCalibration(
                            axis.a,
                            axis.b,
                            axis.scale,
                            axis.offset,
                            diagnostics,
                        )
                    )

        diagnostics = copy.deepcopy(axis.diagnostics)
        diagnostics["anchors"][0]["target_mz"] = 37.034
        with self.assertRaises(ptrms.MassCalibrationError):
            ptrms.validate_mass_axis(
                ptrms.MassAxisCalibration(
                    axis.a, axis.b, axis.scale, axis.offset, diagnostics
                )
            )

        diagnostics = copy.deepcopy(axis.diagnostics)
        persistence = diagnostics["anchors"][1]["persistence"]
        persistence["fraction"] = 1.0
        persistence["accepted_blocks"] = 7
        with self.assertRaises(ptrms.MassCalibrationError):
            ptrms.validate_mass_axis(
                ptrms.MassAxisCalibration(
                    axis.a, axis.b, axis.scale, axis.offset, diagnostics
                )
            )

    def test_mass_axis_progress_is_monotonic_and_cancellable(self):
        progress = []
        with self._file(peaks=self._good_peaks()) as h5:
            axis = ptrms.load_mass_axis(h5, progress=progress.append)
        self.assertTrue(axis.applied)
        self.assertEqual(progress[-1], 1.0)
        self.assertEqual(progress, sorted(progress))

        checks = []
        with self._file(peaks=self._good_peaks()) as h5:
            with self.assertRaises(ptrms.AnalysisCancelled):
                ptrms.load_mass_axis(
                    h5,
                    progress=lambda value: checks.append(value),
                    should_stop=lambda: len(checks) >= 2,
                )
        self.assertGreaterEqual(len(checks), 2)
        self.assertEqual(checks, sorted(checks))

    def test_implausible_affine_solution_is_rejected(self):
        peaks = [
            (37.028405 + 0.175, 1200.0),
            (203.942993 - 0.175, 1000.0),
        ]
        with self._file(peaks=peaks) as h5:
            with self.assertRaises(ptrms.MassCalibrationError) as caught:
                ptrms.load_mass_axis(h5)

        diagnostics = caught.exception.diagnostics
        self.assertIn(diagnostics["anchors"][1]["status"], {"ambiguous", "accepted"})
        self.assertFalse(diagnostics["applied"])


if __name__ == "__main__":
    unittest.main()
