"""Integration tests for measured-shape extraction and isotope-aware quantification."""

import unittest
from unittest import mock

import h5py
import numpy as np
from calibration_helpers import identity_mass_axis

from sniff import isotopes, ptrms


class EmpiricalExtractionTest(unittest.TestCase):
    def test_empirical_cluster_fit_runs_in_the_streaming_and_interval_paths(self):
        a = 10000.0
        bins = 105000
        cycles = 4
        axis = identity_mass_axis(a=a, b=0.0)
        masses = [40.0, 60.0, 80.0, 100.0, 100.025]
        first_amplitudes = np.array([100.0, 50.0, 80.0, 40.0])
        second_amplitudes = np.array([30.0, 80.0, 20.0, 60.0])
        x = np.arange(bins, dtype=np.float64)
        data = np.full((cycles, bins), 2.0, dtype=np.float32)

        def shape(mass):
            centre = a * np.sqrt(mass)
            sigma = ptrms._sigma_tb(mass, a, 2400.0, mass_axis=axis)
            normalised = (x - centre) / sigma
            profile = np.exp(-0.5 * normalised**2)
            profile[normalised > 0] *= np.exp(-0.12 * normalised[normalised > 0])
            return profile

        for mass in masses[:3]:
            data += (200.0 * shape(mass))[None, :]
        data += first_amplitudes[:, None] * shape(masses[3])[None, :]
        data += second_amplitudes[:, None] * shape(masses[4])[None, :]

        diagnostics = {}
        with h5py.File("in-memory", "w", driver="core", backing_store=False) as h5:
            h5.create_dataset("SPECdata/Intensities", data=data)
            h5.create_dataset("SPECdata/AverageSpec", data=data.mean(axis=0))
            traces, _ = ptrms.extract_traces(
                h5,
                masses,
                mass_axis=axis,
                peak_fit_model="empirical-v1",
                fit_diagnostics=diagnostics,
                per_range={"sample_01": (1, 2), "sample_02": (3, 4)},
                block=2,
            )
            profile_x = np.linspace(-4.5, 4.5, 181)
            unavailable_profile = {
                "usable": False,
                "reason": "no clean references",
                "n_reference_peaks": 0,
                "x": profile_x,
                "y": np.exp(-0.5 * profile_x**2),
            }
            fallback_diagnostics = {}
            with mock.patch.object(
                ptrms.peak_fit,
                "estimate_empirical_profile",
                return_value=unavailable_profile,
            ):
                fallback_traces, _ = ptrms.extract_traces(
                    h5,
                    masses,
                    mass_axis=axis,
                    peak_fit_model="empirical-v1",
                    fit_diagnostics=fallback_diagnostics,
                    block=2,
                )

        self.assertTrue(diagnostics["profile"]["usable"])
        self.assertEqual(diagnostics["clusters"][0]["method"], "empirical-v1")
        self.assertEqual(diagnostics["clusters"][0]["status"], "reliable")
        self.assertEqual(
            set(diagnostics["clusters"][0]["ranges"]),
            {"sample_01", "sample_02"},
        )
        first_trace = traces[masses[3]][0]
        second_trace = traces[masses[4]][0]
        self.assertAlmostEqual(first_trace[0] / first_trace[1], 2.0, delta=0.15)
        self.assertAlmostEqual(second_trace[0] / second_trace[1], 0.375, delta=0.05)
        fallback = fallback_diagnostics["clusters"][0]
        self.assertEqual(fallback["method"], "gaussian-fallback-v1")
        self.assertNotEqual(fallback["status"], "legacy")
        if fallback["status"] == "unresolved":
            self.assertTrue(np.isnan(fallback_traces[masses[3]][0]).all())

    def test_joint_temporal_fit_uses_one_whole_run_design(self):
        a = 10000.0
        bins = 105000
        cycles = 64
        axis = identity_mass_axis(a=a, b=0.0)
        masses = [40.0, 60.0, 80.0, 100.0, 100.025]
        x = np.arange(bins, dtype=np.float64)
        data = np.full((cycles, bins), 2.0, dtype=np.float32)

        def shape(mass):
            centre = a * np.sqrt(mass)
            sigma = ptrms._sigma_tb(mass, a, 2400.0, mass_axis=axis)
            normalised = (x - centre) / sigma
            return np.exp(-0.5 * normalised**2)

        for mass in masses[:3]:
            data += (200.0 * shape(mass))[None, :]
        cycle = np.arange(cycles, dtype=np.float64)
        first = 80.0 + 20.0 * np.sin(cycle / 8.0)
        second = 45.0 + 15.0 * np.cos(cycle / 11.0)
        data += first[:, None] * shape(masses[3])[None, :]
        data += second[:, None] * shape(masses[4])[None, :]
        diagnostics = {}
        progress = []

        with h5py.File("in-memory", "w", driver="core", backing_store=False) as h5:
            h5.create_dataset("SPECdata/Intensities", data=data)
            h5.create_dataset("SPECdata/AverageSpec", data=data.mean(axis=0))
            traces, _ = ptrms.extract_traces(
                h5,
                masses,
                mass_axis=axis,
                peak_fit_model="joint-temporal-v2",
                fit_diagnostics=diagnostics,
                per_range={"sample_01": (1, 32), "sample_02": (33, 64)},
                progress=progress.append,
                block=16,
            )

        report = diagnostics["clusters"][0]
        self.assertEqual(report["method"], "joint-temporal-v2")
        self.assertEqual(report["status"], "reliable")
        self.assertGreater(report["selected_lambda"], 0.0)
        self.assertNotIn("ranges", report)
        self.assertTrue(np.isfinite(traces[masses[3]][0]).all())
        self.assertTrue(np.isfinite(traces[masses[4]][0]).all())
        self.assertEqual(progress[-1], 1.0)
        self.assertEqual(progress, sorted(progress))
        self.assertTrue(any(0.0 < value < 1.0 for value in progress))


def test_gaussian_design_with_duplicate_centres_is_withheld():
    _lo, _hi, projection, _norm = ptrms._cluster_design(
        [59.049, 59.049],
        1000.0,
        0.0,
        nbin=10000,
        mass_axis=identity_mass_axis(a=1000.0, b=0.0),
    )

    assert np.isnan(projection).all()


def test_gaussian_design_with_resolved_centres_stays_finite():
    _lo, _hi, projection, _norm = ptrms._cluster_design(
        [59.049, 59.149],
        1000.0,
        0.0,
        nbin=10000,
        mass_axis=identity_mass_axis(a=1000.0, b=0.0),
    )

    assert np.isfinite(projection).all()


def test_evidence_traces_remove_shared_primary_ion_movement():
    traces = {59.0: np.array([10.0, 20.0, 40.0])}
    diagnostics = {}
    with mock.patch.object(
        ptrms,
        "load_transmission",
        return_value=(np.array([1.0, 200.0]), np.ones(2)),
    ):
        result = ptrms.normalise_evidence_traces(
            traces,
            {59.0: 59.0},
            object(),
            primary=np.array([1.0, 2.0, 4.0]),
            diagnostics=diagnostics,
        )

    np.testing.assert_allclose(result[59.0], [10.0, 10.0, 10.0])
    assert diagnostics["status"] == "available"


def test_evidence_traces_are_withheld_without_compatible_primary():
    traces = {59.0: np.array([10.0, 20.0, 40.0])}
    diagnostics = {}
    with mock.patch.object(
        ptrms,
        "load_transmission",
        return_value=(np.array([1.0, 200.0]), np.ones(2)),
    ):
        result = ptrms.normalise_evidence_traces(
            traces,
            {59.0: 59.0},
            object(),
            primary=np.array([1.0, 2.0]),
            diagnostics=diagnostics,
        )

    assert np.isnan(result[59.0]).all()
    assert diagnostics == {
        "model": "transmission-primary-normalised-evidence-v1",
        "status": "withheld",
        "signal_basis": None,
        "reason": "primary-ion trace is missing or incompatible",
    }


class IsotopeQuantificationTest(unittest.TestCase):
    def test_auxiliary_channels_do_not_create_rows(self):
        source = {"mz": 59.0, "formula": "C3H6O"}
        model = isotopes.formula_isotope_model(source["formula"])
        ratio = model["channels"][0]["ratio"]
        target_mass = source["mz"] + model["channels"][0]["shift"]
        target = {"mz": target_mass, "formula": "CH4O"}
        plan = isotopes.build_isotope_plan([source, target])
        traces = {
            source["mz"]: (np.full(4, 100.0), source["mz"]),
            target_mass: (
                np.full(4, 20.0 + 100.0 * ratio),
                target_mass,
            ),
        }

        with (
            mock.patch.object(
                ptrms,
                "load_transmission",
                return_value=(np.array([1.0, 200.0]), np.ones(2)),
            ),
            mock.patch.object(ptrms, "has_transmission", return_value=True),
        ):
            rows, params = ptrms.quantify(
                traces,
                object(),
                {"sample_01": (1, 4)},
                K=1.0,
                primary=np.ones(4),
                molar_volume=24.0,
                isotope_plan=plan,
            )

        self.assertEqual(len(rows), 2)
        target_row = next(row for row in rows if row["mass"] == target_mass)
        self.assertGreater(target_row["cor"]["Average"], 20.0)
        self.assertAlmostEqual(target_row["con"]["Average"], 20.0)
        self.assertTrue(params["isotopes"]["enabled"])
        self.assertEqual(
            params["isotopes"]["corrections"][1]["status"],
            "spillover-corrected",
        )

    def test_non_analyte_channel_keeps_signal_but_withholds_concentration(self):
        traces = {30.0: (np.full(4, 100.0), 30.0)}
        with (
            mock.patch.object(
                ptrms,
                "load_transmission",
                return_value=(np.array([1.0, 200.0]), np.ones(2)),
            ),
            mock.patch.object(ptrms, "has_transmission", return_value=True),
        ):
            rows, params = ptrms.quantify(
                traces,
                object(),
                {"sample_01": (1, 4)},
                K=1.0,
                primary=np.ones(4),
                molar_volume=24.0,
                non_analyte_masses={30.0},
            )

        self.assertEqual(rows[0]["raw"]["Average"], 100.0)
        self.assertEqual(rows[0]["cor"]["Average"], 100.0)
        self.assertTrue(np.isnan(rows[0]["con"]["Average"]))
        self.assertEqual(params["non_analyte_masses"], [30.0])


if __name__ == "__main__":
    unittest.main()
