from pathlib import Path
import sys
import unittest

import numpy as np


MODEL_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(MODEL_DIR / "phase_field_qsgs"))
from qsgs_enthalpy_solver import QSGSEnthalpyGrid, Ti3SiC2ThermalProperties  # noqa: E402


class QSGSEnthalpySolverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.grid_path = MODEL_DIR / "phase_field_qsgs" / "qsgs_phase_fractions_10A.npz"

    def test_uniform_enthalpy_round_trip_both_boundary_cases(self):
        for boundary in ("A", "B"):
            model = QSGSEnthalpyGrid(self.grid_path, boundary)
            shape = model.temperature_shape
            for temperature, liquid_fraction in ((300.0, 0.0), (900.0, 0.0), (1235.08, 0.37), (1300.0, 1.0), (1800.0, 1.0), (2000.0, 1.0)):
                t = np.full(shape, temperature)
                f = np.full(shape, liquid_fraction)
                h = model.enthalpy_density_from_state(t, f)
                state = model.state_from_enthalpy_density(h)
                active = model.ag_fraction > 1.0e-12
                self.assertTrue(np.allclose(state.temperature_K[active], temperature, atol=2.0e-9))
                self.assertTrue(np.allclose(state.ag_liquid_fraction[active], liquid_fraction, atol=2.0e-10))
                if boundary == "B":
                    self.assertTrue(np.allclose(state.temperature_K, temperature, atol=2.0e-9))

    def test_adiabatic_scaffold_only_cells_have_no_fake_temperature_rise(self):
        model = QSGSEnthalpyGrid(self.grid_path, "A")
        state = model.state_from_enthalpy_density(np.zeros(model.temperature_shape))
        scaffold_only = model.ag_fraction == 0.0
        self.assertTrue(np.all(state.temperature_K[scaffold_only] == 300.0))
        self.assertTrue(np.all(state.ag_liquid_fraction[scaffold_only] == 0.0))

    def test_diffusion_operator_conserves_energy_both_boundary_cases(self):
        rng = np.random.default_rng(8)
        for boundary in ("A", "B"):
            model = QSGSEnthalpyGrid(self.grid_path, boundary)
            temp = 450.0 + 200.0 * rng.random(model.temperature_shape)
            state = model.state_from_enthalpy_density(model.enthalpy_density_from_state(temp, np.zeros_like(temp)))
            diffusion = model.diffusion_power_density_W_m3(state.temperature_K, model.conductivity_W_m_K(state))
            integrated_rate = float(diffusion.sum() * model.cell_volume_m3)
            self.assertLess(abs(integrated_rate), 2.0e-12 * float(np.abs(diffusion).sum() * model.cell_volume_m3))

    def test_source_only_adds_masked_ag_energy_and_closes_balance(self):
        model = QSGSEnthalpyGrid(self.grid_path, "A")
        self.assertGreater(model.actual_source_power_W, 0.0)
        self.assertLess(model.actual_source_power_W, model.peak_flux_W_m2 * np.pi * (model.gaussian_radius_A * 1.0e-10) ** 2)
        dt = min(model.stable_timestep_s(), 1.0e-16)
        row = model.advance(dt, source_on=True)
        self.assertGreater(row["injected_energy_J"], 0.0)
        self.assertAlmostEqual(row["total_energy_J"], row["injected_energy_J"], delta=1.0e-11 * row["injected_energy_J"])

    def test_ti3sic2_properties_extend_to_2000k_and_reject_above(self):
        properties = Ti3SiC2ThermalProperties()
        self.assertAlmostEqual(
            properties.cp_J_kg_K(300.0), 560.3756617000839, delta=1.0e-9
        )
        self.assertAlmostEqual(
            properties.cp_J_kg_K(1473.15), 783.0773165577963, delta=1.0e-9
        )
        self.assertAlmostEqual(
            properties.cp_J_kg_K(2000.0), 798.0792813783802, delta=1.0e-9
        )
        self.assertAlmostEqual(
            properties.conductivity_W_m_K(2000.0), 29.6, delta=1.0e-12
        )
        self.assertAlmostEqual(
            properties.specific_enthalpy_J_kg(2000.0),
            1268886.2334097673,
            delta=1.0e-6,
        )
        with self.assertRaisesRegex(ValueError, "2000.00 K"):
            properties.cp_J_kg_K(2000.01)

    def test_case_b_enthalpy_guard_stops_above_2000k(self):
        model = QSGSEnthalpyGrid(self.grid_path, "B")
        liquid = np.ones(model.temperature_shape)
        with self.assertRaisesRegex(ValueError, "2000.00 K"):
            model.enthalpy_density_from_state(
                np.full(model.temperature_shape, 2000.01), liquid
            )
        enthalpy_at_limit = model.enthalpy_density_from_state(
            np.full(model.temperature_shape, 2000.0), liquid
        )
        with self.assertRaisesRegex(ValueError, "2000.00 K"):
            model.state_from_enthalpy_density(enthalpy_at_limit + 1.0e12)

    def test_explicit_timestep_guard(self):
        model = QSGSEnthalpyGrid(self.grid_path, "A")
        with self.assertRaises(ValueError):
            model.advance(model.stable_timestep_s() * 1.01, source_on=False)


if __name__ == "__main__":
    unittest.main()
