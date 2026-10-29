from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np


OUT_DIR = Path(__file__).resolve().parent
RUN_DIR = OUT_DIR.parent
PROJECT_DIR = OUT_DIR.parents[2]
PHASE_DIR = PROJECT_DIR / "phase_field_qsgs"
sys.path.insert(0, str(PHASE_DIR))

from qsgs_enthalpy_solver import QSGSEnthalpyGrid  # noqa: E402


GEOMETRY = RUN_DIR / "qsgs_morphology_ag70vol.npz"


class Ag70CaseBTests(unittest.TestCase):
    def model(self) -> QSGSEnthalpyGrid:
        return QSGSEnthalpyGrid(
            None,
            "B",
            peak_flux_eV_ps_A2=5.0,
            gaussian_radius_A=30.0,
            absorption_depth_A=20.0,
            native_geometry_npz=GEOMETRY,
            source_geometry_npz=GEOMETRY,
            mesh_refinement=1,
        )

    def energy(self, model: QSGSEnthalpyGrid) -> float:
        return float(model.enthalpy_density_J_m3.sum(dtype=np.float64) * model.cell_volume_m3)

    def test_ag70_geometry_and_source_are_bound_to_ag_only(self) -> None:
        model = self.model()
        self.assertEqual(model.temperature_shape, (100, 100, 50))
        self.assertAlmostEqual(float(model.ag_fraction.mean()), 0.699996, places=6)
        self.assertAlmostEqual(float(model.sc_fraction.mean()), 0.300004, places=6)
        self.assertEqual(model.source_geometry_npz, GEOMETRY)
        self.assertGreater(model.actual_source_power_W, 0.0)
        self.assertFalse(np.any((model.source_cell_power_W > 0.0) & (model.sc_fraction > 0.0)))
        self.assertTrue(np.all(model.state.temperature_K == 300.0))

    def test_uniform_two_phase_state_is_static_without_source(self) -> None:
        model = self.model()
        temperature = np.full(model.temperature_shape, 500.0)
        model.enthalpy_density_J_m3 = model.enthalpy_density_from_state(
            temperature, np.zeros_like(temperature)
        )
        model.state = model.state_from_enthalpy_density(model.enthalpy_density_J_m3)
        initial_energy = self.energy(model)
        model.run(0.001, source_on=False, sample_interval_ps=0.001)
        self.assertLess(float(np.max(np.abs(model.state.temperature_K - 500.0))), 1.0e-8)
        self.assertLess(abs(self.energy(model) - initial_energy), 1.0e-12 * initial_energy)

    def test_ag_to_scaffold_heat_transfer_conserves_system_energy(self) -> None:
        model = self.model()
        temperature = np.where(model.ag_fraction > 0.0, 600.0, 300.0)
        model.enthalpy_density_J_m3 = model.enthalpy_density_from_state(
            temperature, np.zeros_like(temperature)
        )
        model.state = model.state_from_enthalpy_density(model.enthalpy_density_J_m3)
        initial_energy = self.energy(model)
        ag_mean_before = float(np.sum(model.ag_fraction * model.state.temperature_K) / model.ag_fraction.sum())
        sc_max_before = float(model.state.temperature_K[model.sc_fraction > 0.0].max())
        dt = min(model.stable_timestep_s(), 0.01e-15)
        for _ in range(8):
            model.advance(dt, source_on=False, enforce_stability=True)
        final_energy = self.energy(model)
        ag_mean_after = float(np.sum(model.ag_fraction * model.state.temperature_K) / model.ag_fraction.sum())
        sc_max_after = float(model.state.temperature_K[model.sc_fraction > 0.0].max())
        self.assertLess(ag_mean_after, ag_mean_before)
        self.assertGreater(sc_max_after, sc_max_before)
        self.assertLess(abs(final_energy - initial_energy), 1.0e-12 * initial_energy)

    def test_explicit_dt_guard_and_tisc3_upper_property_limit(self) -> None:
        model = self.model()
        with self.assertRaises(ValueError):
            model.advance(model.stable_timestep_s() * 1.01, source_on=False)
        limit = float(model.sc.temperature_K[-1])
        at_limit = np.full(model.temperature_shape, limit)
        h_at_limit = model.enthalpy_density_from_state(at_limit, np.ones_like(at_limit))
        state = model.state_from_enthalpy_density(h_at_limit)
        self.assertLessEqual(float(state.temperature_K[model.sc_fraction > 0.0].max()), limit + 1.0e-9)
        above_limit = np.full(model.temperature_shape, limit + 0.01)
        h_above_limit = model.enthalpy_density_from_state(above_limit, np.ones_like(above_limit))
        with self.assertRaisesRegex(ValueError, "1473.15 K"):
            model.state_from_enthalpy_density(h_above_limit)


if __name__ == "__main__":
    unittest.main()
