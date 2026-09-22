"""Regression coverage for browser-review data preparation."""

import json
import unittest
from unittest import mock

import h5py
import numpy as np
from calibration_helpers import identity_mass_axis

from sniff import ptrms, viz


class VizDataTest(unittest.TestCase):
    def setUp(self):
        self._calibration = mock.patch.object(
            ptrms, "load_mass_axis", return_value=identity_mass_axis()
        )
        self._calibration.start()

    def tearDown(self):
        self._calibration.stop()

    def test_nonfinite_average_spectrum_bins_are_zero_filled(self):
        with h5py.File("in-memory", "w", driver="core", backing_store=False) as h5:
            h5.create_dataset("SPECdata/Intensities", data=np.zeros((2, 5)))
            h5.create_dataset(
                "SPECdata/AverageSpec",
                data=np.array([1.2, np.nan, np.inf, -np.inf, 4.8]),
            )
            with (
                mock.patch.object(ptrms, "load_mass_cal", return_value=(10.0, 1.0)),
                mock.patch.object(
                    ptrms,
                    "load_transmission",
                    return_value=(np.array([1.0]), np.array([1.0])),
                ),
                mock.patch.object(ptrms, "spec_duration_s", return_value=1.0),
                mock.patch.object(ptrms, "extract_primary", return_value=None),
                mock.patch.object(ptrms, "water_cluster_ratio", return_value=None),
                mock.patch.object(
                    ptrms, "build_discriminator", return_value=np.ones(2)
                ),
                mock.patch.object(
                    ptrms, "derive_molar_volume_info", return_value=(24.465, "test")
                ),
                mock.patch.object(ptrms, "derive_K", return_value=None),
                mock.patch.object(ptrms, "resolve_k", return_value={}),
                mock.patch.object(
                    ptrms, "load_rate_constants", return_value={"compounds": []}
                ),
            ):
                data = viz.build_viz_data(h5, peaks_cfg=[], ranges_cfg=[])

        self.assertEqual(data["spectrum"], [1, 0, 0, 0, 5])
        self.assertIn("mass_axis_calibration", data["meta"])
        self.assertTrue(data["meta"]["mass_axis_calibration"]["applied"])
        self.assertEqual(data["candidate_coverage"]["total_components"], 0)
        json.dumps(data, allow_nan=False)

    def test_peak_abundance_is_mean_integrated_raw_signal(self):
        with h5py.File("in-memory", "w", driver="core", backing_store=False) as h5:
            h5.create_dataset("SPECdata/Intensities", data=np.zeros((2, 5)))
            h5.create_dataset("SPECdata/AverageSpec", data=np.ones(5))
            with (
                mock.patch.object(ptrms, "load_mass_cal", return_value=(10.0, 1.0)),
                mock.patch.object(
                    ptrms,
                    "load_transmission",
                    return_value=(np.array([1.0]), np.array([1.0])),
                ),
                mock.patch.object(ptrms, "spec_duration_s", return_value=1.0),
                mock.patch.object(ptrms, "extract_primary", return_value=None),
                mock.patch.object(ptrms, "water_cluster_ratio", return_value=None),
                mock.patch.object(
                    ptrms, "build_discriminator", return_value=np.ones(2)
                ),
                mock.patch.object(
                    ptrms,
                    "derive_molar_volume_info",
                    return_value=(24.465, "test"),
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
                mock.patch.object(viz.formula_id, "score_peak", return_value=[]),
            ):
                data = viz.build_viz_data(
                    h5,
                    peaks_cfg=[{"mz": 2.0}],
                    ranges_cfg=[],
                )

        self.assertEqual(data["peaks"][0]["abundance"], 3.0)

    def _payload(self, peaks_cfg):
        """The review payload for ``peaks_cfg``, with the file layer mocked out."""
        with h5py.File("in-memory", "w", driver="core", backing_store=False) as h5:
            h5.create_dataset("SPECdata/Intensities", data=np.zeros((2, 5)))
            h5.create_dataset("SPECdata/AverageSpec", data=np.ones(5))
            with (
                mock.patch.object(ptrms, "load_mass_cal", return_value=(10.0, 1.0)),
                mock.patch.object(
                    ptrms,
                    "load_transmission",
                    return_value=(np.array([1.0]), np.array([1.0])),
                ),
                mock.patch.object(ptrms, "spec_duration_s", return_value=1.0),
                mock.patch.object(ptrms, "extract_primary", return_value=None),
                mock.patch.object(ptrms, "water_cluster_ratio", return_value=None),
                mock.patch.object(
                    ptrms, "build_discriminator", return_value=np.ones(2)
                ),
                mock.patch.object(
                    ptrms,
                    "derive_molar_volume_info",
                    return_value=(24.465, "test"),
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
                mock.patch.object(viz.formula_id, "score_peak", return_value=[]),
            ):
                return viz.build_viz_data(h5, peaks_cfg=peaks_cfg, ranges_cfg=[])

    def test_refined_candidates_fill_only_available_unique_fresh_review_blanks(self):
        config_peaks = [
            {"mz": 59.049, "label": "acetone", "formula": "C3H6O"},
            {"mz": 44.062, "label": ""},
            {"mz": 44.063, "label": ""},
        ]
        ethenamine = {
            "formula": "C2H5N",
            "name": "ethenamine",
            "probability": 0.9,
            "k": 2.1,
            "k_estimated": True,
            "flags": ["fragmentation reported"],
        }
        review_peaks = [
            {
                "_config_original": dict(config_peaks[0]),
                "mz": 59.049,
                "abundance": 200.0,
                "label": "acetone",
                "formula": "C3H6O",
                "candidates": [],
            },
            {
                "_config_original": dict(config_peaks[1]),
                "mz": 44.062,
                "abundance": 100.0,
                "label": "m44.062",
                "labelAuto": "m44.062",
                "formula": "",
                "candidates": [
                    {
                        "formula": "C3H6O",
                        "name": "acetone",
                        "probability": 0.1,
                        "k": 1.0,
                        "k_estimated": False,
                        "flags": [],
                    },
                    ethenamine,
                ],
            },
            {
                "_config_original": dict(config_peaks[2]),
                "mz": 44.063,
                "abundance": 10.0,
                "label": "m44.063",
                "labelAuto": "m44.063",
                "formula": "",
                "candidates": [{**ethenamine, "probability": 0.5}],
            },
        ]

        viz._apply_refined_identity_defaults(review_peaks, config_peaks)

        self.assertEqual(review_peaks[0]["label"], "acetone")
        self.assertEqual(review_peaks[1]["label"], "ethenamine")
        self.assertEqual(review_peaks[1]["formula"], "C2H5N")
        self.assertEqual(review_peaks[1]["k"], 2.1)
        self.assertTrue(review_peaks[1]["k_estimated"])
        self.assertEqual(review_peaks[1]["flags"], ["fragmentation reported"])
        self.assertNotIn("labelAuto", review_peaks[1])
        self.assertEqual(config_peaks[1]["label"], "ethenamine")
        self.assertEqual(config_peaks[1]["formula"], "C2H5N")
        self.assertEqual(review_peaks[2]["label"], "m44.063")
        self.assertIn("labelAuto", review_peaks[2])

    def test_provisional_candidate_cannot_win_refined_identity_assignment(self):
        config_peaks = [
            {"mz": 59.049, "label": ""},
            {"mz": 59.050, "label": ""},
        ]
        review_peaks = [
            {
                "_config_original": dict(config_peaks[0]),
                "mz": 59.049,
                "abundance": 100.0,
                "label": "m59.049",
                "labelAuto": "m59.049",
                "formula": "",
                "candidates": [
                    {
                        "formula": "C3H6O",
                        "name": "acetone",
                        "probability": 0.9,
                        "assignment_eligible": True,
                    },
                    {
                        "formula": "C2H2O2",
                        "probability": 0.8,
                        "assignment_eligible": True,
                    },
                ],
            },
            {
                "_config_original": dict(config_peaks[1]),
                "mz": 59.050,
                "abundance": 80.0,
                "label": "m59.050",
                "labelAuto": "m59.050",
                "formula": "",
                "candidates": [
                    {
                        "formula": "C3H6O",
                        "name": "acetone",
                        "probability": 0.99,
                        "assignment_eligible": False,
                        "mass_match": "broad-proposal",
                    }
                ],
            },
        ]

        viz._apply_refined_identity_defaults(review_peaks, config_peaks)

        self.assertEqual(review_peaks[0]["formula"], "C3H6O")
        self.assertEqual(review_peaks[0]["label"], "acetone")
        self.assertEqual(review_peaks[1]["formula"], "")
        self.assertEqual(review_peaks[1]["label"], "m59.050")
        self.assertIn("labelAuto", review_peaks[1])

    def test_a_name_the_tool_invented_is_never_saved_as_an_assignment(self):
        # An unnamed peak still needs something to draw on the spectrum, so a
        # mass-derived stand-in stands in for display. It must not reach the config:
        # a page autosave would otherwise turn 132 empty names into 132 names that
        # look assigned, and the CSV would print the mass twice.
        unnamed = self._payload([{"mz": 2.0}])["peaks"][0]
        named = self._payload(
            [{"mz": 2.0, "label": "ethanol", "formula": "C2H6O"}]
        )["peaks"][0]

        self.assertEqual(unnamed["label"], "m2.000")
        self.assertEqual(unnamed["labelAuto"], "m2.000")
        self.assertEqual(named["label"], "ethanol")
        self.assertNotIn("labelAuto", named)

    def test_irregular_pctimes_make_relative_and_absolute_axes(self):
        with h5py.File("in-memory", "w", driver="core", backing_store=False) as h5:
            h5.create_dataset("SPECdata/Intensities", data=np.zeros((3, 2)))
            h5.create_dataset(
                "SPECdata/PCTime",
                data=[[1_700_000_000], [1_700_000_002], [1_700_000_009]],
            )
            h5.attrs["Single Spec Duration (ms)"] = [1000.0]
            axes = ptrms.viz_x_axis_data(h5)

        self.assertEqual(axes["relative"], [0.0, 2.0, 9.0])
        self.assertEqual(axes["absolute"], [1700000000.0, 1700000002.0, 1700000009.0])
        self.assertTrue(axes["absolute_available"])

    def test_absolute_axis_applies_file_lab_timezone_offset(self):
        with h5py.File("in-memory", "w", driver="core", backing_store=False) as h5:
            h5.create_dataset("SPECdata/Intensities", data=np.zeros((2, 2)))
            h5.create_dataset("SPECdata/PCTime", data=[[1000.0], [1002.0]])
            h5.attrs["Single Spec Duration (ms)"] = [1000.0]
            h5.attrs["UTC_Offset"] = [3600.0]
            axes = ptrms.viz_x_axis_data(h5)

        self.assertEqual(axes["relative"], [0.0, 2.0])
        self.assertEqual(axes["absolute"], [4600.0, 4602.0])
        self.assertEqual(axes["absolute_offset_s"], 3600.0)

    def test_adjacent_sub_millisecond_pctimes_keep_distinct_axis_values(self):
        pctimes = [
            1_700_000_000.0004,
            1_700_000_000.0008,
            1_700_000_000.0012,
        ]
        with h5py.File("in-memory", "w", driver="core", backing_store=False) as h5:
            h5.create_dataset("SPECdata/Intensities", data=np.zeros((3, 2)))
            h5.create_dataset("SPECdata/PCTime", data=np.asarray(pctimes)[:, None])
            h5.attrs["Single Spec Duration (ms)"] = [1.0]
            axes = ptrms.viz_x_axis_data(h5)

        self.assertEqual(axes["absolute"], pctimes)
        self.assertTrue(np.all(np.diff(axes["absolute"]) > 0))
        json.loads(json.dumps(axes, allow_nan=False))

    def test_absolute_axis_accepts_year_zero_but_rejects_expanded_years(self):
        year_zero = -62167219200.0
        year_10000 = 253402300800.0
        cases = (
            ([year_zero, year_zero + 0.001], True),
            ([year_zero - 0.001, year_zero], False),
            ([year_10000 - 0.001, year_10000 - 0.0004], True),
            ([year_10000, year_10000 + 0.001], False),
        )
        for pctimes, available in cases:
            with self.subTest(pctimes=pctimes):
                with h5py.File(
                    "in-memory", "w", driver="core", backing_store=False
                ) as h5:
                    h5.create_dataset("SPECdata/Intensities", data=np.zeros((2, 2)))
                    h5.create_dataset(
                        "SPECdata/PCTime", data=np.asarray(pctimes)[:, None]
                    )
                    h5.attrs["Single Spec Duration (ms)"] = [1.0]
                    axes = ptrms.viz_x_axis_data(h5)

                self.assertEqual(axes["absolute_available"], available)
                if not available:
                    self.assertIsNone(axes["absolute"])

    def test_render_html_rejects_ambiguous_embedded_absolute_axis(self):
        data = {
            "meta": {
                "x_axis": {
                    "absolute": [1.0, 1.0],
                    "absolute_available": True,
                }
            }
        }
        with self.assertRaisesRegex(ValueError, "strictly increasing"):
            viz.render_html(data)

    def test_malformed_pctimes_disable_absolute_and_use_duration_fallback(self):
        with h5py.File("in-memory", "w", driver="core", backing_store=False) as h5:
            h5.create_dataset("SPECdata/Intensities", data=np.zeros((3, 2)))
            h5.create_dataset(
                "SPECdata/PCTime",
                data=[[1_700_000_000], [np.nan], [1_700_000_009]],
            )
            h5.attrs["Single Spec Duration (ms)"] = [2500.0]
            axes = ptrms.viz_x_axis_data(h5)

        self.assertIsNone(axes["absolute"])
        self.assertFalse(axes["absolute_available"])
        self.assertEqual(axes["relative"], [0.0, 2.5, 5.0])

    def test_nonpositive_or_nonfinite_duration_uses_safe_relative_domain(self):
        for duration_ms in (0.0, -1000.0, np.nan):
            with self.subTest(duration_ms=duration_ms):
                with h5py.File(
                    "in-memory", "w", driver="core", backing_store=False
                ) as h5:
                    h5.create_dataset("SPECdata/Intensities", data=np.zeros((3, 2)))
                    h5.attrs["Single Spec Duration (ms)"] = [duration_ms]
                    axes = ptrms.viz_x_axis_data(h5)

                self.assertEqual(axes["relative"], [0.0, 1.0, 2.0])
                self.assertTrue(np.all(np.isfinite(axes["relative"])))
                self.assertTrue(np.all(np.diff(axes["relative"]) > 0))

    def test_out_of_javascript_date_range_keeps_relative_but_disables_absolute(self):
        with h5py.File("in-memory", "w", driver="core", backing_store=False) as h5:
            h5.create_dataset("SPECdata/Intensities", data=np.zeros((3, 2)))
            h5.create_dataset(
                "SPECdata/PCTime",
                data=[[10_000_000_000_000], [10_000_000_000_002], [10_000_000_000_005]],
            )
            h5.attrs["Single Spec Duration (ms)"] = [1000.0]
            axes = ptrms.viz_x_axis_data(h5)

        self.assertIsNone(axes["absolute"])
        self.assertFalse(axes["absolute_available"])
        self.assertEqual(axes["relative"], [0.0, 2.0, 5.0])
        self.assertTrue(np.all(np.diff(axes["relative"]) > 0))


if __name__ == "__main__":
    unittest.main()
