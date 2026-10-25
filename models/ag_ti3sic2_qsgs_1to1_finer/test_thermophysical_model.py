import sys
import unittest
from pathlib import Path

import numpy as np

PHASE_DIR = Path(__file__).resolve().parent / "phase_field_qsgs"
sys.path.insert(0, str(PHASE_DIR))

from stefan_benchmark import run_one_phase_melting  # noqa: E402
from thermophysical_model import AgEnthalpy, Ti3SiC2ThermalProperties  # noqa: E402


class AgEnthalpyTests(unittest.TestCase):
    def setUp(self):
        self.ag = AgEnthalpy()

    def test_sensible_and_liquid_enthalpy_round_trip(self):
        for temperature, fraction in ((300.0, 0.0), (800.0, 0.0), (1200.0, 0.0), (1300.0, 1.0)):
            with self.subTest(temperature=temperature):
                h = self.ag.specific_enthalpy_J_kg(temperature, fraction)
                state = self.ag.state_from_specific_enthalpy(h)
                self.assertAlmostEqual(float(state.temperature_K), temperature, places=8)
                self.assertAlmostEqual(float(state.liquid_fraction), fraction, places=12)

    def test_latent_heat_is_an_isothermal_enthalpy_plateau(self):
        tm = self.ag.melting_temperature_K
        h_s = self.ag.specific_enthalpy_J_kg(tm, 0.0)
        h_mid = self.ag.specific_enthalpy_J_kg(tm, 0.5)
        h_l = self.ag.specific_enthalpy_J_kg(tm, 1.0)
        self.assertAlmostEqual(h_l - h_s, self.ag.latent_heat_J_kg, places=8)
        for h, expected_fraction in ((h_s, 0.0), (h_mid, 0.5), (h_l, 1.0)):
            state = self.ag.state_from_specific_enthalpy(h)
            self.assertAlmostEqual(float(state.temperature_K), tm, places=10)
            self.assertAlmostEqual(float(state.liquid_fraction), expected_fraction, places=12)

    def test_vector_inverse_and_phase_fraction_are_monotone(self):
        enthalpy = np.array([
            self.ag.specific_enthalpy_J_kg(900.0, 0.0),
            self.ag.specific_enthalpy_J_kg(self.ag.melting_temperature_K, 0.25),
            self.ag.specific_enthalpy_J_kg(self.ag.melting_temperature_K, 0.75),
            self.ag.specific_enthalpy_J_kg(1300.0, 1.0),
        ])
        state = self.ag.state_from_specific_enthalpy(enthalpy)
        np.testing.assert_allclose(state.temperature_K, [900.0, self.ag.melting_temperature_K,
                                                         self.ag.melting_temperature_K, 1300.0], atol=1.0e-8)
        np.testing.assert_allclose(state.liquid_fraction, [0.0, 0.25, 0.75, 1.0], atol=1.0e-12)


class CompositePropertyTests(unittest.TestCase):
    def test_ti3sic2_mass_specific_heat_conversion_uses_formula_mass(self):
        scaffold = Ti3SiC2ThermalProperties()
        self.assertAlmostEqual(scaffold.molar_mass_kg_mol, 0.195708, places=9)
        self.assertAlmostEqual(scaffold.cp_J_kg_K(300.0), 110.0 / 0.195708, places=9)
        self.assertAlmostEqual(scaffold.cp_J_kg_K(1473.15), 155.0 / 0.195708, places=9)

    def test_ti3sic2_candidate_refuses_high_temperature_extrapolation(self):
        scaffold = Ti3SiC2ThermalProperties()
        with self.assertRaises(ValueError):
            scaffold.conductivity_W_m_K(1473.16)

    def test_stefan_front_and_energy_balance(self):
        ag = AgEnthalpy()
        for dx, cells in ((40.0e-6, 250), (20.0e-6, 500)):
            with self.subTest(dx=dx):
                result = run_one_phase_melting(dx, 0.01, 100.0, cells, ag)
                self.assertLess(abs(result.relative_interface_error), 0.02)
                self.assertLess(abs(result.energy_balance_residual_J_m2), 1.0e-5)


if __name__ == "__main__":
    unittest.main()
