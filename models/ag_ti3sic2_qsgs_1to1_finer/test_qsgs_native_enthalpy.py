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
from qsgs_liquid_tracking import LiquidClusterTracker  # noqa: E402


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

    def test_checkpoint_restart_preserves_fields_energy_and_tracker_identity(self):
        ag = np.ones((4, 4, 3), dtype=bool)
        self.save_geometry(ag)
        original = QSGSEnthalpyGrid(None, "A", native_geometry_npz=self.geometry_path)
        tracker = LiquidClusterTracker(original.spacing_A, original.size_A)
        tracker.update(0.0, original.state.ag_liquid_fraction, original.ag_fraction > 0.0)
        dt = original.stable_timestep_s() * 0.2
        original.advance(dt, source_on=True)
        tracker.update(original.time_s * 1.0e12, original.state.ag_liquid_fraction, original.ag_fraction > 0.0)
        checkpoint = Path(self.tempdir.name) / "restart.npz"
        original.save_checkpoint(checkpoint, tracker)

        restarted = QSGSEnthalpyGrid(None, "A", native_geometry_npz=self.geometry_path)
        restored_tracker = LiquidClusterTracker(restarted.spacing_A, restarted.size_A)
        restored = restarted.load_checkpoint(checkpoint, restored_tracker)
        self.assertAlmostEqual(restored["time_ps"], original.time_s * 1.0e12, places=15)
        np.testing.assert_array_equal(restarted.enthalpy_density_J_m3, original.enthalpy_density_J_m3)
        np.testing.assert_array_equal(restarted.state.temperature_K, original.state.temperature_K)
        np.testing.assert_array_equal(restarted.state.ag_liquid_fraction, original.state.ag_liquid_fraction)
        np.testing.assert_array_equal(restored_tracker.previous_track_labels, tracker.previous_track_labels)

        original.advance(dt, source_on=False)
        restarted.advance(dt, source_on=False)
        np.testing.assert_array_equal(restarted.enthalpy_density_J_m3, original.enthalpy_density_J_m3)
        np.testing.assert_array_equal(restarted.state.temperature_K, original.state.temperature_K)
        self.assertEqual(restarted.injected_energy_J, original.injected_energy_J)
        self.assertLess(abs(restored["energy_balance_error_J"]), 1.0e-25)

    def test_checkpoint_rejects_boundary_mismatch(self):
        self.save_geometry(np.ones((3, 3, 3), dtype=bool))
        model_a = QSGSEnthalpyGrid(None, "A", native_geometry_npz=self.geometry_path)
        checkpoint = Path(self.tempdir.name) / "case-a.npz"
        model_a.save_checkpoint(checkpoint)
        model_b = QSGSEnthalpyGrid(None, "B", native_geometry_npz=self.geometry_path)
        with self.assertRaisesRegex(ValueError, "identity does not match"):
            model_b.load_checkpoint(checkpoint)

    def test_latent_energy_ledger_matches_liquid_fraction_field(self):
        self.save_geometry(np.ones((3, 3, 3), dtype=bool))
        model = QSGSEnthalpyGrid(None, "A", native_geometry_npz=self.geometry_path)
        liquid = np.full(model.temperature_shape, 0.25)
        temperature = np.full(model.temperature_shape, model.melting_temperature_K)
        model.enthalpy_density_J_m3 = model.enthalpy_density_from_state(temperature, liquid)
        model.state = model.state_from_enthalpy_density(model.enthalpy_density_J_m3)
        expected = float(np.sum(model._latent_volumetric_J_m3 * liquid) * model.cell_volume_m3)
        self.assertAlmostEqual(model._ag_latent_energy_J(), expected, delta=1.0e-30)

    def test_run_split_at_checkpoint_matches_uninterrupted_run(self):
        self.save_geometry(np.ones((4, 4, 3), dtype=bool))
        uninterrupted = QSGSEnthalpyGrid(None, "A", native_geometry_npz=self.geometry_path)
        split = QSGSEnthalpyGrid(None, "A", native_geometry_npz=self.geometry_path)
        tracker_full = LiquidClusterTracker(uninterrupted.spacing_A, uninterrupted.size_A)
        tracker_split = LiquidClusterTracker(split.spacing_A, split.size_A)
        dt = uninterrupted.stable_timestep_s() * 0.2
        first_duration_ps = 5.0 * dt * 1.0e12
        total_duration_ps = 2.0 * first_duration_ps
        uninterrupted.run(
            total_duration_ps,
            source_on=True,
            sample_interval_ps=first_duration_ps,
            maximum_dt_s=dt,
            track_liquid_connectivity=True,
            tracker=tracker_full,
        )
        checkpoint_dir = Path(self.tempdir.name) / "split-run"
        split.run(
            first_duration_ps,
            source_on=True,
            sample_interval_ps=first_duration_ps,
            maximum_dt_s=dt,
            track_liquid_connectivity=True,
            tracker=tracker_split,
            checkpoint_times_ps=[first_duration_ps],
            checkpoint_dir=checkpoint_dir,
            checkpoint_prefix="heat",
        )

        restarted = QSGSEnthalpyGrid(None, "A", native_geometry_npz=self.geometry_path)
        restored_tracker = LiquidClusterTracker(restarted.spacing_A, restarted.size_A)
        restarted.load_checkpoint(next(checkpoint_dir.glob("heat_*.npz")), restored_tracker)
        restarted.run(
            first_duration_ps,
            source_on=True,
            sample_interval_ps=first_duration_ps,
            maximum_dt_s=dt,
            track_liquid_connectivity=True,
            tracker=restored_tracker,
        )
        np.testing.assert_allclose(
            restarted.enthalpy_density_J_m3,
            uninterrupted.enthalpy_density_J_m3,
            rtol=1.0e-13,
            atol=1.0e-20,
        )
        np.testing.assert_allclose(
            restarted.state.temperature_K,
            uninterrupted.state.temperature_K,
            rtol=0.0,
            atol=1.0e-10,
        )
        self.assertAlmostEqual(restarted.injected_energy_J, uninterrupted.injected_energy_J, delta=1.0e-30)
        self.assertAlmostEqual(restarted.time_s, uninterrupted.time_s, delta=1.0e-27)


if __name__ == "__main__":
    unittest.main()
