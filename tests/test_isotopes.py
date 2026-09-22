"""Tests for formula-derived isotope channels and guarded correction."""

import unittest

import numpy as np

from sniff import isotopes


class IsotopeModelTest(unittest.TestCase):
    def test_carbon_m_plus_one_ratio_and_exact_shift(self):
        model = isotopes.formula_isotope_model("C")

        channel = model["channels"][0]
        self.assertAlmostEqual(channel["shift"], 1.0034, places=3)
        self.assertGreater(channel["ratio"], 0.010)
        self.assertLess(channel["ratio"], 0.012)
        self.assertLess(model["monoisotopic_fraction"], 1.0)

    def test_halogen_has_strong_m_plus_two_channel(self):
        model = isotopes.formula_isotope_model("C2H3Cl")

        self.assertGreater(model["channels"][1]["ratio"], 0.30)

    def test_silicon_and_iodine_have_versioned_isotope_models(self):
        silicon = isotopes.formula_isotope_model("C2H6OSi")
        iodine = isotopes.formula_isotope_model("C6H5I")

        self.assertGreater(silicon["channels"][0]["ratio"], 0.04)
        self.assertEqual(iodine["version"], isotopes.MODEL_VERSION)

    def test_formula_parser_rejects_unsupported_or_malformed_values(self):
        for formula in ("", "C2H5+", "C0H2", "XeH4"):
            with self.subTest(formula=formula):
                with self.assertRaises(ValueError):
                    isotopes.parse_formula(formula)

    def test_plan_adds_auxiliary_channels_without_analyte_rows(self):
        peaks = [{"mz": 59.0491, "formula": "C3H6O"}, {"mz": 80.0}]

        plan = isotopes.build_isotope_plan(peaks)

        self.assertEqual(len(plan["parents"]), 1)
        self.assertEqual(len(plan["parents"][0]["channels"]), 2)
        self.assertEqual(plan["parents"][0]["mz"], 59.0491)
        self.assertGreater(len(plan["extraction_masses"]), len(peaks))

    def test_spillover_is_removed_in_corrected_signal_space(self):
        source = {"mz": 59.0, "formula": "C3H6O"}
        source_model = isotopes.formula_isotope_model(source["formula"])
        ratio = source_model["channels"][0]["ratio"]
        target_mass = source["mz"] + source_model["channels"][0]["shift"]
        target = {"mz": target_mass, "formula": "CH4O"}
        plan = isotopes.build_isotope_plan([source, target])
        corrected = {
            source["mz"]: np.full(4, 100.0),
            target_mass: np.full(4, 20.0 + 100.0 * ratio),
        }

        net, diagnostics = isotopes.correct_parent_signals(corrected, plan)

        np.testing.assert_allclose(net[target_mass], 20.0)
        target_diagnostic = next(
            item for item in diagnostics if item["mz"] == target_mass
        )
        self.assertEqual(target_diagnostic["status"], "spillover-corrected")

    def test_abundance_scaling_happens_after_all_spillover_subtraction(self):
        source = {"mz": 59.0, "formula": "C20H20"}
        source_model = isotopes.formula_isotope_model(source["formula"])
        ratio = source_model["channels"][0]["ratio"]
        target_mass = source["mz"] + source_model["channels"][0]["shift"]
        target = {"mz": target_mass, "formula": "C2H6O"}
        target_model = isotopes.formula_isotope_model(target["formula"])
        plan = isotopes.build_isotope_plan([source, target])
        corrected = {
            source["mz"]: np.full(3, 100.0),
            target_mass: np.full(3, 20.0 + 100.0 * ratio),
        }

        net, _ = isotopes.correct_parent_signals(
            corrected, plan, abundance_basis="total"
        )

        np.testing.assert_allclose(
            net[target_mass],
            20.0 / target_model["monoisotopic_fraction"],
        )

    def test_abundance_scaling_preserves_only_unavailable_cycles(self):
        parent = {"mz": 59.0, "formula": "C2H6O"}
        model = isotopes.formula_isotope_model(parent["formula"])
        corrected = {parent["mz"]: np.array([10.0, np.nan, 30.0])}

        net, diagnostics = isotopes.correct_parent_signals(
            corrected,
            isotopes.build_isotope_plan([parent]),
            abundance_basis="total",
        )

        expected = np.array(
            [
                10.0 / model["monoisotopic_fraction"],
                np.nan,
                30.0 / model["monoisotopic_fraction"],
            ]
        )
        np.testing.assert_allclose(net[parent["mz"]], expected, equal_nan=True)
        self.assertTrue(diagnostics[0]["abundance_applied"])

    def test_isotope_observation_clustering_is_peak_order_independent(self):
        peaks = [
            {"mz": 100.0, "formula": "C6H12O"},
            {"mz": 100.02, "formula": "C5H10O2"},
        ]

        forward = isotopes.build_isotope_plan(peaks)
        reverse = isotopes.build_isotope_plan(list(reversed(peaks)))

        self.assertEqual(forward["extraction_masses"], reverse["extraction_masses"])
        forward_channels = sorted(
            (parent["mz"], channel["order"], channel["observation_mz"])
            for parent in forward["parents"]
            for channel in parent["channels"]
        )
        reverse_channels = sorted(
            (parent["mz"], channel["order"], channel["observation_mz"])
            for parent in reverse["parents"]
            for channel in parent["channels"]
        )
        self.assertEqual(forward_channels, reverse_channels)

    def test_weighted_nonnegative_solver_refits_after_negative_component(self):
        design = np.array([[1.0, 1.0], [1.0, 0.0]])
        observed = np.array([1.0, 2.0])

        estimate, active = isotopes._solve_weighted_nonnegative(design, observed)

        np.testing.assert_allclose(estimate, [1.5, 0.0])
        self.assertEqual(active, [0])
        self.assertTrue((estimate >= 0).all())

    def test_joint_envelope_withholds_parent_on_nonnegative_boundary(self):
        plan = {
            "version": isotopes.ENVELOPE_MODEL_VERSION,
            "parents": [
                {
                    "mz": 100.0,
                    "formula": "CH2O",
                    "monoisotopic_fraction": 1.0,
                    "channels": [
                        {
                            "observation_mz": 101.0,
                            "ratio": 1.0,
                            "order": 1,
                        }
                    ],
                },
                {
                    "mz": 101.0,
                    "formula": "CH4O",
                    "monoisotopic_fraction": 1.0,
                    "channels": [],
                },
            ],
        }
        corrected = {
            100.0: np.full(4, 2.0),
            101.0: np.full(4, 1.0),
        }

        net, diagnostics = isotopes.correct_parent_signals(corrected, plan)

        target = next(item for item in diagnostics if item["mz"] == 101.0)
        self.assertNotEqual(target["status"], "envelope-fitted")
        self.assertTrue(np.isnan(net[101.0]).all())

    def test_joint_envelope_fit_recovers_overlapping_parents(self):
        source = {"mz": 100.0, "formula": "C6H12O"}
        source_model = isotopes.formula_isotope_model(source["formula"])
        target = {
            "mz": source["mz"] + source_model["channels"][0]["shift"],
            "formula": "C5H10O2",
        }
        plan = isotopes.build_isotope_plan(
            [source, target],
            model=isotopes.ENVELOPE_MODEL_VERSION,
        )
        amplitudes = {
            source["mz"]: np.array([100.0, 120.0, 80.0, 110.0]),
            target["mz"]: np.array([20.0, 25.0, 18.0, 23.0]),
        }
        corrected = {
            mass: np.zeros(4, dtype=np.float64)
            for mass in plan["extraction_masses"]
        }
        for parent in plan["parents"]:
            signal = amplitudes[parent["mz"]]
            corrected[parent["mz"]] += signal
            for channel in parent["channels"]:
                corrected[channel["observation_mz"]] += signal * channel["ratio"]

        net, diagnostics = isotopes.correct_parent_signals(corrected, plan)

        np.testing.assert_allclose(net[source["mz"]], amplitudes[source["mz"]])
        np.testing.assert_allclose(net[target["mz"]], amplitudes[target["mz"]])
        self.assertTrue(all(item["status"] == "envelope-fitted" for item in diagnostics))
        self.assertTrue(all(item["model"] == "formula-envelope-v2" for item in diagnostics))
        self.assertTrue(all(item["median_relative_uncertainty"] > 0 for item in diagnostics))

    def test_joint_envelope_fit_withholds_rank_deficient_component(self):
        peaks = [
            {"mz": 100.0, "formula": "C6H12O"},
            {"mz": 100.0, "formula": "C6H12O"},
        ]
        plan = isotopes.build_isotope_plan(
            peaks,
            model=isotopes.ENVELOPE_MODEL_VERSION,
        )
        corrected = {
            mass: np.full(4, 100.0)
            for mass in plan["extraction_masses"]
        }

        net, diagnostics = isotopes.correct_parent_signals(corrected, plan)

        self.assertTrue(np.isnan(net[100.0]).all())
        self.assertTrue(
            all(
                item["status"] == "withheld-ill-conditioned-envelope"
                for item in diagnostics
            )
        )

    def test_negative_spillover_solution_is_withheld(self):
        source = {"mz": 59.0, "formula": "C20H20"}
        source_model = isotopes.formula_isotope_model(source["formula"])
        target_mass = source["mz"] + source_model["channels"][0]["shift"]
        target = {"mz": target_mass, "formula": "CH4O"}
        plan = isotopes.build_isotope_plan([source, target])
        corrected = {
            source["mz"]: np.full(3, 100.0),
            target_mass: np.full(3, 1.0),
        }

        net, diagnostics = isotopes.correct_parent_signals(corrected, plan)

        self.assertTrue(np.isnan(net[target_mass]).all())
        target_diagnostic = next(
            item for item in diagnostics if item["mz"] == target_mass
        )
        self.assertIn("withheld", target_diagnostic["status"])


if __name__ == "__main__":
    unittest.main()
