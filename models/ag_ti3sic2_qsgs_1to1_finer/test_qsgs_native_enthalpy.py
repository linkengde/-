from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np


MODEL_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(MODEL_DIR / "phase_field_qsgs"))
from qsgs_enthalpy_solver import (  # noqa: E402
    QSGSEnthalpyGrid,
    ThermalState,
    longest_sampled_liquid_presence_ps,
)


class NativeQSGSEnthalpyTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.geometry_path = Path(self.tempdir.name) / "geometry.npz"

    def tearDown(self):
        self.tempdir.cleanup()

    def save_geometry(self, ag_mask):
        np.savez_compressed(
            self.geometry_path,
            skeleton=~np.asarray(ag_mask, dtype=bool),
            size_A=np.asarray(ag_mask.shape, dtype=float),
        )
        return self.geometry_path

    def test_native_grid_uses_binary_geometry_and_exact_source(self):
        ag = np.zeros((3, 2, 2), dtype=bool)
        ag[0, 0, 0] = True
        ag[1, 0, 0] = True
        self.save_geometry(ag)
        model = QSGSEnthalpyGrid(None, "A", native_geometry_npz=self.geometry_path)
        self.assertEqual(model.temperature_shape, ag.shape)
        np.testing.assert_array_equal(model.ag_fraction, ag.astype(float))
        self.assertIn("exact separable integral", model.source_integration_method)
        self.assertGreater(model.actual_source_power_W, 0.0)

    def test_refinement_subdivides_existing_voxels_without_changing_geometry_or_source_power(self):
        ag = np.zeros((3, 2, 2), dtype=bool)
        ag[0, 0, 0] = True
        ag[1, 0, 0] = True
        self.save_geometry(ag)
        coarse = QSGSEnthalpyGrid(None, "A", native_geometry_npz=self.geometry_path)
        fine = QSGSEnthalpyGrid(None, "A", native_geometry_npz=self.geometry_path, mesh_refinement=2)
        self.assertEqual(fine.temperature_shape, (6, 4, 4))
        self.assertTrue(np.allclose(fine.spacing_A, coarse.spacing_A / 2.0))
        self.assertAlmostEqual(fine.actual_source_power_W, coarse.actual_source_power_W, delta=1.0e-14)
        self.assertAlmostEqual(float(fine.ag_fraction.mean()), float(coarse.ag_fraction.mean()))

    def test_case_a_does_not_conduct_between_diagonal_ag_voxels(self):
        ag = np.zeros((3, 3, 3), dtype=bool)
        ag[0, 0, 0] = True
        ag[1, 1, 1] = True
        self.save_geometry(ag)
        model = QSGSEnthalpyGrid(None, "A", native_geometry_npz=self.geometry_path)
        temperature = np.full(ag.shape, 300.0)
        temperature[0, 0, 0] = 500.0
        conductivity = model.conductivity_W_m_K(ThermalState(temperature, np.zeros(ag.shape)))
        power = model.diffusion_power_density_W_m3(temperature, conductivity)
        self.assertEqual(np.count_nonzero(power), 0)

    def test_case_b_uses_ideal_scaffold_contact_and_conserves_pairwise_flux(self):
        ag = np.zeros((3, 3, 3), dtype=bool)
        ag[1, 1, 1] = True
        self.save_geometry(ag)
        model = QSGSEnthalpyGrid(None, "B", native_geometry_npz=self.geometry_path)
        temperature = np.full(ag.shape, 300.0)
        temperature[1, 1, 1] = 400.0
        conductivity = model.conductivity_W_m_K(ThermalState(temperature, np.zeros(ag.shape)))
        power = model.diffusion_power_density_W_m3(temperature, conductivity)
        self.assertLess(power[1, 1, 1], 0.0)
        self.assertGreater(np.count_nonzero(power > 0.0), 0)
        self.assertLess(abs(float(power.sum())), 1.0e-14 * float(np.abs(power).sum()))

    def test_sampled_liquid_presence_duration_resets_on_a_missing_sample(self):
        rows = [
            {"time_ps": 0.18, "connected_liquid_cluster_count": 1},
            {"time_ps": 0.19, "connected_liquid_cluster_count": 1},
            {"time_ps": 0.20, "connected_liquid_cluster_count": 0},
            {"time_ps": 0.21, "connected_liquid_cluster_count": 1},
            {"time_ps": 0.22, "connected_liquid_cluster_count": 1},
        ]
        self.assertAlmostEqual(longest_sampled_liquid_presence_ps(rows, 0.01), 0.01)


if __name__ == "__main__":
    unittest.main()
