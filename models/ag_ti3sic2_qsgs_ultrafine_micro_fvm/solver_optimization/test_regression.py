#!/usr/bin/env python3
"""Small equivalence checks for kernels, nonlinear policies, and restart state."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from checkpointing import CheckpointStore
from diagnostics import phase_metrics
from implicit_enthalpy_fvm import BoundarySet, ImplicitEnthalpyFVM
from material_model import CandidateMaterials
from optimized_solver import FastCandidateMaterials, NumbaEnthalpyFVM
from profiled_advance import advance_profiled


class RegressionTests(unittest.TestCase):
    def setUp(self):
        from numba import set_num_threads
        set_num_threads(1)
        self.rng = np.random.default_rng(217)
        self.shape = (9, 8, 7)
        self.ag = np.ascontiguousarray(self.rng.uniform(0.0, 1.0, self.shape))
        self.temperature = np.ascontiguousarray(self.rng.uniform(295.0, 1300.0, self.shape))
        self.materials = CandidateMaterials()

    def test_fused_material_fields_match_reference(self):
        fast = FastCandidateMaterials()
        from numba_kernels import material_fields_3d
        optimized = material_fields_3d(self.temperature, self.ag, fast._kernel_parameters)
        reference = (
            self.materials.enthalpy_density_J_m3(self.temperature, self.ag),
            self.materials.volumetric_capacity_J_m3_K(self.temperature, self.ag),
            self.materials.conductivity_W_m_K(self.temperature, self.ag),
        )
        for actual, expected in zip(optimized, reference):
            relative = np.max(np.abs(actual - expected) / np.maximum(np.abs(expected), 1.0))
            self.assertLess(relative, 2.0e-14)

    def test_face_and_residual_operators_match_reference(self):
        base = ImplicitEnthalpyFVM(self.ag, (1.0e-7,) * 3, materials=CandidateMaterials())
        fast = NumbaEnthalpyFVM(self.ag, (1.0e-7,) * 3, materials=CandidateMaterials())
        base.temperature_K[...] = self.temperature
        fast.temperature_K[...] = self.temperature
        base.enthalpy_density_J_m3[...] = base.materials.enthalpy_density_J_m3(self.temperature, self.ag)
        fast.enthalpy_density_J_m3[...] = fast.materials.enthalpy_density_J_m3(self.temperature, self.ag)
        boundary = BoundarySet(dirichlet_K={"x_low": 300.0, "y_high": 300.0, "z_low": 300.0})
        source = np.zeros(self.shape)
        r0, g0, _, f0 = base._residual(self.temperature, base.enthalpy_density_J_m3, 1.0e-9, source, boundary)
        r1, g1, _, f1 = fast._residual(self.temperature, fast.enthalpy_density_J_m3, 1.0e-9, source, boundary)
        for a, b in zip(g0, g1):
            self.assertLess(np.max(np.abs(a - b) / np.maximum(np.abs(a), 1.0)), 3.0e-14)
        self.assertLess(np.max(np.abs(f0["conduction"] - f1["conduction"]))
                        / max(float(np.max(np.abs(f0["conduction"]))), 1.0), 3.0e-14)
        self.assertLess(np.max(np.abs(r0 - r1)) / max(float(np.max(np.abs(r0))), 1.0), 3.0e-14)

    def _new_solver(self, solver_class=NumbaEnthalpyFVM):
        material = FastCandidateMaterials() if solver_class is NumbaEnthalpyFVM else CandidateMaterials()
        return solver_class(self.ag, (1.0e-7,) * 3, materials=material,
                            maximum_nonlinear_iterations=48, maximum_linear_iterations=400)

    def test_profiled_reference_step_matches_original_advance(self):
        baseline = self._new_solver(ImplicitEnthalpyFVM)
        profiled = self._new_solver(ImplicitEnthalpyFVM)
        source = np.zeros(self.shape)
        source[3:6, 3:6, 2:5] = 4.0e18 * self.ag[3:6, 3:6, 2:5]
        base_record = baseline.advance(2.0e-10, source)
        profile_record, details = advance_profiled(profiled, 2.0e-10, source, linear_policy="fixed")
        np.testing.assert_allclose(profiled.temperature_K, baseline.temperature_K, rtol=2.0e-13, atol=1.0e-12)
        self.assertAlmostEqual(profiled.source_energy_J, baseline.source_energy_J, places=24)
        self.assertEqual(profile_record["pcg_iterations_total"], base_record["pcg_iterations_total"])
        self.assertEqual(len(details["newton_iterations"]), profile_record["nonlinear_iterations"])

    def test_numba_step_matches_reference_and_ew_closes_energy(self):
        source = np.zeros(self.shape)
        source[3:6, 3:6, 2:5] = 4.0e18 * self.ag[3:6, 3:6, 2:5]
        baseline = self._new_solver(ImplicitEnthalpyFVM)
        optimized = self._new_solver(NumbaEnthalpyFVM)
        r0, _ = advance_profiled(baseline, 2.0e-10, source, linear_policy="fixed")
        r1, prof = advance_profiled(optimized, 2.0e-10, source, linear_policy="fixed")
        np.testing.assert_allclose(optimized.temperature_K, baseline.temperature_K, rtol=2.0e-12, atol=2.0e-10)
        self.assertLess(abs(r1["energy_residual_J"]), 1.0e-18)
        adaptive = self._new_solver(NumbaEnthalpyFVM)
        r2, adaptive_log = advance_profiled(adaptive, 2.0e-10, source, linear_policy="ew")
        self.assertLess(abs(r2["energy_residual_J"]), 1.0e-15)
        self.assertLessEqual(r2["nonlinear_relative_residual"], adaptive.nonlinear_rtol * 1.01)
        self.assertEqual(len(adaptive_log["newton_iterations"]), r2["nonlinear_iterations"])

    def test_checkpoint_restart_matches_continuous_advance(self):
        source = np.zeros(self.shape)
        source[3:6, 3:6, 2:5] = 4.0e18 * self.ag[3:6, 3:6, 2:5]
        continuous = self._new_solver()
        first, profiles = advance_profiled(continuous, 2.0e-10, source)
        config = {"solver_parameters": {"nonlinear_rtol": 1.0e-8, "linear_rtol": 1.0e-8,
                                         "max_newton": 48, "max_pcg": 400,
                                         "interface_resistance_m2_K_W": 0.0}}
        with tempfile.TemporaryDirectory() as temp:
            store = CheckpointStore(temp, continuous, source, config)
            path = store.save("heating", 3.0e-10, first, profiles,
                              {"phase": "heating", "phase_elapsed_s": 2.0e-10})
            with np.load(path, allow_pickle=False) as saved:
                self.assertEqual(saved["temperature_K"].dtype, np.float64)
                self.assertEqual(saved["enthalpy_density_J_m3"].dtype, np.float64)
            resumed = self._new_solver()
            CheckpointStore.restore(path, resumed, source)
            advance_profiled(continuous, 3.0e-10, source)
            advance_profiled(resumed, 3.0e-10, source)
            np.testing.assert_array_equal(resumed.temperature_K, continuous.temperature_K)
            np.testing.assert_array_equal(resumed.enthalpy_density_J_m3, continuous.enthalpy_density_J_m3)
            self.assertEqual(resumed.source_energy_J, continuous.source_energy_J)
            self.assertEqual(resumed.time_s, continuous.time_s)

    def test_maximum_ag_liquid_fraction_excludes_ti_only_cells(self):
        ag = np.zeros((3, 3, 3), dtype=np.float64)
        ag[0, 0, 0] = 1.0
        sim = NumbaEnthalpyFVM(ag, (1.0e-7,) * 3, materials=FastCandidateMaterials())
        sim.temperature_K[2, 2, 2] = sim.materials.Tl_K + 1.0  # Ti-only cell
        self.assertEqual(phase_metrics(sim)["Ag_max_liquid_fraction"], 0.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
