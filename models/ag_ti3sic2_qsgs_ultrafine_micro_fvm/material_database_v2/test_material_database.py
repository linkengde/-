#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""文献物性模型与Numba实现的独立数值检查。"""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
PHASE = PROJECT / "phase_change_validation"
SOLVER = PROJECT / "solver_optimization"
for _path in (str(SOLVER), str(PHASE), str(PROJECT)):
    if _path in sys.path:
        sys.path.remove(_path)
    sys.path.insert(0, _path)

from material_database_v2.fast_materials import FastPankratzMaterials
from material_database_v2.material_model import PankratzMaterials


class MaterialDatabaseTests(unittest.TestCase):
    def setUp(self):
        self.material = PankratzMaterials(1.0)

    def test_primary_constants_and_phase_conductivities(self):
        m = self.material
        self.assertAlmostEqual(m.Tm_K, 1235.08, places=8)
        self.assertAlmostEqual(m.Lf_J_mol, 11296.8, places=8)
        self.assertAlmostEqual(m.Lf_J_kg, 104727.8067123, places=5)
        self.assertAlmostEqual(m.ag_cp_liquid_J_kg_K, 310.3046125, places=5)
        self.assertAlmostEqual(m.ag_conductivity_W_m_K(m.Ts_K), 358.0, delta=0.05)
        self.assertAlmostEqual(m.ag_conductivity_W_m_K(m.Tl_K), 175.0, delta=0.05)

    def test_enthalpy_inverse_across_solid_mushy_and_liquid(self):
        temperatures = np.array([300.0, 600.0, 1000.0, 1234.4, 1234.8,
                                 1235.08, 1235.4, 1235.8, 1300.0, 1800.0])
        enthalpy = self.material.ag_specific_enthalpy_J_kg(temperatures)
        recovered = self.material.ag_temperature_from_enthalpy_J_kg(enthalpy)
        np.testing.assert_allclose(recovered, temperatures, rtol=0.0, atol=2e-12)

    def test_enthalpy_derivative_matches_apparent_heat_capacity(self):
        for temperature in (500.0, 1000.0, 1234.9, 1235.08, 1235.3, 1300.0):
            eps = 1.0e-4
            derivative = (self.material.ag_specific_enthalpy_J_kg(temperature + eps)
                          - self.material.ag_specific_enthalpy_J_kg(temperature - eps)) / (2 * eps)
            capacity = self.material.ag_apparent_cp_J_kg_K(temperature)
            self.assertAlmostEqual(derivative, capacity, delta=max(2e-3, abs(capacity) * 2e-7))

    def test_latent_heat_is_added_once(self):
        m = self.material
        sensible_integral = 0.5 * (m.ag_cp_solid_mid + m.ag_cp_liquid_J_kg_K) * m.width_K
        measured = m.ag_specific_enthalpy_J_kg(m.Tl_K) - m.ag_specific_enthalpy_J_kg(m.Ts_K)
        self.assertAlmostEqual(measured - sensible_integral, m.Lf_J_kg, delta=1e-8)

    def test_mushy_width_sensitivity_preserves_full_latent_heat(self):
        for width in (0.5, 1.0, 2.0):
            m = PankratzMaterials(width)
            sensible = 0.5 * (m.ag_cp_solid_mid + m.ag_cp_liquid_J_kg_K) * width
            jump = m.ag_specific_enthalpy_J_kg(m.Tl_K) - m.ag_specific_enthalpy_J_kg(m.Ts_K)
            self.assertAlmostEqual(jump - sensible, m.Lf_J_kg, delta=2e-8)

    def test_material_ranges_reject_silent_extrapolation(self):
        with self.assertRaises(ValueError):
            self.material.ag_specific_enthalpy_J_kg(2000.01)
        with self.assertRaises(ValueError):
            self.material.ti_conductivity_W_m_K(1473.16)

    def test_numba_property_fields_match_python_reference(self):
        fast = FastPankratzMaterials(1.0)
        shape = (9, 8, 7)
        temperature = np.linspace(300.0, 1400.0, np.prod(shape)).reshape(shape).copy()
        ag = np.linspace(0.0, 1.0, np.prod(shape)).reshape(shape).copy()
        compiled = fast._evaluate(temperature, ag)
        reference = (
            self.material.enthalpy_density_J_m3(temperature, ag),
            self.material.volumetric_capacity_J_m3_K(temperature, ag),
            self.material.conductivity_W_m_K(temperature, ag),
        )
        for actual, expected in zip(compiled, reference):
            np.testing.assert_allclose(actual, expected, rtol=2e-14, atol=2e-12)

    def test_single_control_volume_energy_and_300k_hold(self):
        from optimized_solver import NumbaEnthalpyFVM
        from implicit_enthalpy_fvm import BoundarySet

        ag = np.ones((1, 1, 1), dtype=np.float64)
        spacing = (2.0e-6, 2.0e-6, 2.0e-6)
        sim = NumbaEnthalpyFVM(ag, spacing, materials=FastPankratzMaterials())
        quiet = sim.advance(1.0e-8, 0.0,
                            BoundarySet(dirichlet_K={"x_low": 300.0, "x_high": 300.0,
                                                     "y_low": 300.0, "y_high": 300.0,
                                                     "z_low": 300.0, "z_high": 300.0}))
        self.assertEqual(float(sim.temperature_K[0, 0, 0]), 300.0)
        self.assertEqual(float(quiet["energy_residual_J"]), 0.0)
        source = 1.0e15
        dt = 1.0e-7
        heated = sim.advance(dt, source)
        self.assertGreater(float(sim.temperature_K[0, 0, 0]), 300.0)
        self.assertAlmostEqual(float(sim.source_energy_J), source * dt * np.prod(spacing), delta=1e-24)
        self.assertLess(abs(float(heated["energy_residual_J"])) / max(abs(float(sim.source_energy_J)), 1e-30), 1e-8)

    def test_v2_checkpoint_restore_and_restart_are_identical(self):
        import tempfile
        from checkpointing import CheckpointStore
        from optimized_solver import NumbaEnthalpyFVM

        ag = np.full((4, 4, 4), 0.7, dtype=np.float64)
        spacing = (1.0e-6, 1.0e-6, 1.0e-6)
        source = np.full(ag.shape, 1.0e15, dtype=np.float64)
        config = {"materials": "Pankratz v2", "source": "uniform test"}
        sim1 = NumbaEnthalpyFVM(ag, spacing, materials=FastPankratzMaterials())
        sim1.advance(1.0e-8, source)
        with tempfile.TemporaryDirectory() as directory:
            store = CheckpointStore(directory, sim1, source, config)
            checkpoint = store.save("heating", 1.0e-8, {"test": True}, [], {"phase": "heating"})
            sim2 = NumbaEnthalpyFVM(ag, spacing, materials=FastPankratzMaterials())
            CheckpointStore.restore(checkpoint, sim2, source)
            np.testing.assert_array_equal(sim2.temperature_K, sim1.temperature_K)
            np.testing.assert_array_equal(sim2.enthalpy_density_J_m3, sim1.enthalpy_density_J_m3)
            self.assertEqual(sim2.time_s, sim1.time_s)
            sim1.advance(1.0e-8, source)
            sim2.advance(1.0e-8, source)
            np.testing.assert_allclose(sim2.temperature_K, sim1.temperature_K, rtol=0.0, atol=1e-12)


if __name__ == "__main__":
    unittest.main(verbosity=2)
