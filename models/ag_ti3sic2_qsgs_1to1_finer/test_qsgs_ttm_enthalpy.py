from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np


MODEL_DIR = Path(__file__).resolve().parent
PHASE_DIR = MODEL_DIR / "phase_field_qsgs"
sys.path.insert(0, str(PHASE_DIR))
from qsgs_ttm_enthalpy import EV_J, QSGSDoubleTemperatureGrid, field_snapshot_metrics  # noqa: E402
from qsgs_enthalpy_solver import PHASE_FRACTION_ROUNDOFF_TOL, QSGSEnthalpyGrid  # noqa: E402
from qsgs_face_contacts import exact_face_contact_fractions  # noqa: E402


class QSGSDoubleTemperatureTests(unittest.TestCase):
    def make_geometry(self, ag_fraction: np.ndarray, directory: str) -> Path:
        ag = np.asarray(ag_fraction, dtype=float)
        path = Path(directory) / "geometry.npz"
        np.savez_compressed(
            path,
            ag_fraction=ag,
            skeleton_fraction=1.0 - ag,
            size_A=np.asarray(ag.shape, dtype=float) * 10.0,
            spacing_A=np.full(3, 10.0),
        )
        return path

    def test_zero_dimensional_electron_lattice_exchange_conserves_energy(self):
        for boundary, fraction in (("A", 1.0), ("B", 0.5)):
            with tempfile.TemporaryDirectory() as directory:
                geometry = self.make_geometry(np.full((1, 1, 1), fraction), directory)
                model = QSGSDoubleTemperatureGrid(
                    geometry,
                    boundary,
                    peak_flux_eV_ps_A2=0.0,
                    initial_electron_temperature_K=5000.0,
                )
                initial_total = model.total_energy_J
                record = model.advance(0.001, heat_on=False)
                self.assertGreater(record["electron_to_lattice_energy_J"], 0.0)
                self.assertAlmostEqual(
                    model.total_energy_J,
                    initial_total,
                    delta=max(1.0e-30, abs(initial_total) * 1.0e-11),
                )
                self.assertGreater(model.lattice_energy_J, 0.0)
                self.assertLess(abs(record["energy_balance_residual_J"]), 1.0e-29)

    def test_one_dimensional_electron_diffusion_is_conservative(self):
        with tempfile.TemporaryDirectory() as directory:
            geometry = self.make_geometry(np.ones((5, 1, 1)), directory)
            initial_te = np.asarray([1200.0, 300.0, 300.0, 300.0, 300.0]).reshape(5, 1, 1)
            model = QSGSDoubleTemperatureGrid(
                geometry,
                "A",
                peak_flux_eV_ps_A2=0.0,
                initial_electron_temperature_K=initial_te,
            )
            initial_electron = model.electron_energy_J
            before = model.electron.temperature_K.copy()
            diffusion_rate, conductance_sum = model.electron._diffusion_rate(before)
            del diffusion_rate
            cap = model.electron.f_ag * model.electron.cell_volume_A3 * model.parameters.ce(before)
            explicit_limit_ps = float(np.min(cap[conductance_sum > 0.0] / conductance_sum[conductance_sum > 0.0]))
            dt_ps = explicit_limit_ps * 20.0
            record = model.advance(dt_ps, heat_on=False, coupling_on=False)
            after = model.electron.temperature_K
            self.assertLess(after[0, 0, 0], before[0, 0, 0])
            self.assertGreater(after[1, 0, 0], before[1, 0, 0])
            self.assertAlmostEqual(model.electron_energy_J, initial_electron, delta=1.0e-29)
            self.assertLess(abs(record["energy_balance_residual_J"]), 1.0e-29)

    def test_one_dimensional_lattice_diffusion_is_conservative(self):
        with tempfile.TemporaryDirectory() as directory:
            geometry = self.make_geometry(np.zeros((5, 1, 1)), directory)
            model = QSGSDoubleTemperatureGrid(
                geometry,
                "B",
                peak_flux_eV_ps_A2=0.0,
                initial_electron_temperature_K=300.0,
            )
            initial_temperature = np.asarray([400.0, 300.0, 300.0, 300.0, 300.0]).reshape(5, 1, 1)
            model.lattice.enthalpy_density_J_m3 = model.lattice.enthalpy_density_from_state(
                initial_temperature, np.zeros_like(initial_temperature)
            )
            model.lattice.state = model.lattice.state_from_enthalpy_density(
                model.lattice.enthalpy_density_J_m3
            )
            initial_enthalpy = model.lattice_energy_J
            step = 2.0e-1
            self.assertGreater(step, model.lattice_diffusion_timestep_ps * 20.0)
            record = model.advance(step, heat_on=False, coupling_on=False)
            final_temperature = model.lattice.state.temperature_K
            self.assertLess(final_temperature[0, 0, 0], initial_temperature[0, 0, 0])
            self.assertGreater(final_temperature[1, 0, 0], initial_temperature[1, 0, 0])
            self.assertAlmostEqual(model.lattice_energy_J, initial_enthalpy, delta=1.0e-28)
            self.assertLess(abs(record["lattice_diffusion_net_J"]), 1.0e-29)
            self.assertLess(abs(record["energy_balance_residual_J"]), 1.0e-29)

    def test_two_temperature_ag_lattice_conductivity_does_not_duplicate_total_k(self):
        geometry = PHASE_DIR / "qsgs_phase_fractions_10A.npz"
        model = QSGSDoubleTemperatureGrid(geometry, "A", peak_flux_eV_ps_A2=0.0)
        self.assertEqual(float(np.max(model._lattice_conductivity())), 0.0)
        one_temperature = QSGSEnthalpyGrid(geometry, "A", peak_flux_eV_ps_A2=0.0)
        self.assertAlmostEqual(
            float(one_temperature.conductivity_W_m_K().max()), 429.0, delta=0.1
        )
        mixed = QSGSDoubleTemperatureGrid(geometry, "B", peak_flux_eV_ps_A2=0.0)
        expected = mixed.lattice.sc_fraction * mixed.lattice.sc.conductivity_W_m_K(
            mixed.lattice.state.temperature_K
        )
        np.testing.assert_allclose(mixed._lattice_conductivity(), expected, rtol=0.0, atol=0.0)

    def test_active_cell_enthalpy_inverse_matches_full_grid_inverse(self):
        with tempfile.TemporaryDirectory() as directory:
            geometry = self.make_geometry(np.asarray([0.4, 0.6, 0.8]).reshape(3, 1, 1), directory)
            model = QSGSDoubleTemperatureGrid(geometry, "B", peak_flux_eV_ps_A2=0.0)
            temperature = np.asarray([500.0, 1235.08, 1300.0]).reshape(3, 1, 1)
            liquid = np.asarray([0.0, 0.4, 1.0]).reshape(3, 1, 1)
            enthalpy = model.lattice.enthalpy_density_from_state(temperature, liquid)
            expected = model.lattice.state_from_enthalpy_density(enthalpy).temperature_K.ravel()
            actual = model._active_lattice_temperature_from_enthalpy(
                enthalpy.ravel(), model.lattice.ag_fraction.ravel() > 0.0
            )
            np.testing.assert_allclose(actual, expected, rtol=0.0, atol=1.0e-10)

    def test_implicit_lattice_diffusion_conserves_energy_at_large_step(self):
        with tempfile.TemporaryDirectory() as directory:
            geometry = self.make_geometry(np.zeros((5, 1, 1)), directory)
            model = QSGSDoubleTemperatureGrid(
                geometry, "B", peak_flux_eV_ps_A2=0.0, initial_electron_temperature_K=300.0
            )
            initial_temperature = np.full((5, 1, 1), 300.0)
            initial_temperature[0, 0, 0] = 450.0
            model.lattice.enthalpy_density_J_m3 = model.lattice.enthalpy_density_from_state(
                initial_temperature, np.zeros_like(initial_temperature)
            )
            model.lattice.state = model.lattice.state_from_enthalpy_density(
                model.lattice.enthalpy_density_J_m3
            )
            initial_energy = model.lattice_energy_J
            record = model.advance(2.0e-1, heat_on=False, coupling_on=False)
            final_temperature = model.lattice.state.temperature_K
            self.assertLess(final_temperature[0, 0, 0], initial_temperature[0, 0, 0])
            self.assertGreater(final_temperature[1, 0, 0], initial_temperature[1, 0, 0])
            self.assertAlmostEqual(model.lattice_energy_J, initial_energy, delta=2.0e-27)
            self.assertLess(abs(record["lattice_diffusion_net_J"]), 2.0e-27)
            self.assertLess(abs(record["energy_balance_residual_J"]), 2.0e-27)

    def test_bottom_300k_sink_removes_energy_and_is_counted_in_ledger(self):
        with tempfile.TemporaryDirectory() as directory:
            ag_geometry = self.make_geometry(np.ones((1, 1, 3)), directory)
            electron_model = QSGSDoubleTemperatureGrid(
                ag_geometry,
                "A",
                peak_flux_eV_ps_A2=0.0,
                initial_electron_temperature_K=np.asarray([1000.0, 300.0, 300.0]).reshape(1, 1, 3),
                bottom_sink_temperature_K=300.0,
            )
            initial_total = electron_model.total_energy_J
            electron_record = electron_model.advance(0.01, heat_on=False, coupling_on=False)
            self.assertLess(electron_record["electron_diffusion_net_J"], 0.0)
            self.assertLess(electron_model.total_energy_J, initial_total)
            self.assertLess(abs(electron_record["energy_balance_residual_J"]), 1.0e-28)

            sc_directory = Path(directory) / "sc"
            sc_directory.mkdir()
            sc_geometry = self.make_geometry(np.zeros((1, 1, 3)), str(sc_directory))
            lattice_model = QSGSDoubleTemperatureGrid(
                sc_geometry,
                "B",
                peak_flux_eV_ps_A2=0.0,
                initial_electron_temperature_K=300.0,
                bottom_sink_temperature_K=300.0,
            )
            initial_temperature = np.asarray([500.0, 400.0, 300.0]).reshape(1, 1, 3)
            lattice_model.lattice.enthalpy_density_J_m3 = lattice_model.lattice.enthalpy_density_from_state(
                initial_temperature, np.zeros_like(initial_temperature)
            )
            lattice_model.lattice.state = lattice_model.lattice.state_from_enthalpy_density(
                lattice_model.lattice.enthalpy_density_J_m3
            )
            initial_total = lattice_model.total_energy_J
            lattice_record = lattice_model.advance(0.01, heat_on=False, coupling_on=False)
            self.assertLess(lattice_record["lattice_diffusion_net_J"], 0.0)
            self.assertLess(lattice_model.total_energy_J, initial_total)
            self.assertLess(abs(lattice_record["energy_balance_residual_J"]), 1.0e-28)

    def test_run_captures_exact_3d_fields_and_hotspot_metrics(self):
        with tempfile.TemporaryDirectory() as directory:
            geometry = self.make_geometry(np.ones((2, 2, 2)), directory)
            model = QSGSDoubleTemperatureGrid(geometry, "A", peak_flux_eV_ps_A2=0.0)
            model.run(
                0.001,
                0.0005,
                heat_on=False,
                sample_interval_ps=0.001,
                snapshot_times_ps=(0.0005, 0.001),
            )
            self.assertEqual(sorted(model.field_snapshots), [0.0005, 0.001])
            snapshot = model.field_snapshots[0.0005]
            self.assertEqual(snapshot["electron_temperature_K"].shape, (2, 2, 2))
            self.assertEqual(snapshot["electron_energy_density_eV_A3"].shape, (2, 2, 2))
            self.assertEqual(snapshot["lattice_enthalpy_density_J_m3"].shape, (2, 2, 2))
            self.assertEqual(
                snapshot["electron_to_lattice_cumulative_energy_density_J_m3"].shape,
                (2, 2, 2),
            )
            metrics = field_snapshot_metrics(
                0.0005,
                snapshot,
                model.lattice.ag_fraction,
                model.lattice.spacing_A,
                model.lattice.size_A,
            )
            self.assertEqual(metrics["time_ps"], 0.0005)
            self.assertEqual(metrics["ag_liquid_volume_fraction"], 0.0)
            self.assertGreaterEqual(metrics["electron_peak_K"], 300.0)
            self.assertEqual(metrics["electron_ag_weighted_q0.5_K"], 300.0)
            self.assertEqual(metrics["electron_ag_volume_fraction_above_300_plus_5K"], 0.0)
            self.assertEqual(metrics["electron_outside_source_ag_weighted_mean_K"], 300.0)

    def test_xy_tile_keeps_the_original_qsgs_geometry_at_the_center(self):
        with tempfile.TemporaryDirectory() as directory:
            geometry = self.make_geometry(np.asarray([[[1.0]], [[0.5]]]), directory)
            base = QSGSDoubleTemperatureGrid(geometry, "A", peak_flux_eV_ps_A2=0.0)
            tiled = QSGSDoubleTemperatureGrid(
                geometry, "A", peak_flux_eV_ps_A2=0.0, tile_xy=3
            )
            self.assertEqual(tiled.lattice.temperature_shape, (6, 3, 1))
            self.assertEqual(tiled.lattice.size_A.tolist(), [60.0, 30.0, 10.0])
            np.testing.assert_array_equal(tiled.lattice.ag_fraction[2:4, 1:2, :], base.lattice.ag_fraction)
            np.testing.assert_allclose(
                tiled.lattice.source_W_m3[2:4, 1:2, :], base.lattice.source_W_m3, rtol=0.0, atol=1.0e-12
            )

    def test_roundoff_only_phase_slivers_are_removed_from_thermal_domain(self):
        with tempfile.TemporaryDirectory() as directory:
            geometry = self.make_geometry(np.asarray([1.0e-16, 0.5]).reshape(2, 1, 1), directory)
            model = QSGSEnthalpyGrid(geometry, "A", peak_flux_eV_ps_A2=0.0)
            self.assertEqual(model.ag_fraction[0, 0, 0], 0.0)
            self.assertEqual(model.sc_fraction[0, 0, 0], 1.0)
            self.assertEqual(model.ag_fraction[1, 0, 0], 0.5)
            self.assertLessEqual(PHASE_FRACTION_ROUNDOFF_TOL, 1.0e-12)

    def test_exact_face_contacts_remove_spurious_electron_link_conservatively(self):
        with tempfile.TemporaryDirectory() as directory:
            geometry = self.make_geometry(np.full((2, 1, 1), 1.0 / 3.0), directory)
            source_ag = np.zeros((6, 1, 1), dtype=bool)
            source_ag[[1, 4], 0, 0] = True
            coarse_ag = np.full((2, 1, 1), 1.0 / 3.0)
            size_A = np.asarray([20.0, 10.0, 10.0])
            spacing_A = np.asarray([10.0, 10.0, 10.0])
            contacts = exact_face_contact_fractions(source_ag, coarse_ag, size_A, spacing_A)
            self.assertEqual(float(np.sum(contacts[0])), 0.0)
            face_file = Path(directory) / "contacts.npz"
            np.savez_compressed(
                face_file,
                contact_x=contacts[0],
                contact_y=contacts[1],
                contact_z=contacts[2],
                ag_fraction=coarse_ag,
                size_A=size_A,
                spacing_A=spacing_A,
            )
            temperatures = np.asarray([1200.0, 300.0]).reshape(2, 1, 1)
            corrected = QSGSDoubleTemperatureGrid(
                geometry,
                "A",
                peak_flux_eV_ps_A2=0.0,
                initial_electron_temperature_K=temperatures,
                face_contacts_npz=face_file,
            )
            legacy = QSGSDoubleTemperatureGrid(
                geometry,
                "A",
                peak_flux_eV_ps_A2=0.0,
                initial_electron_temperature_K=temperatures,
            )
            corrected_initial_energy = corrected.electron_energy_J
            corrected._implicit_electron_diffusion(1.0)
            np.testing.assert_array_equal(corrected.electron.temperature_K, temperatures)
            self.assertAlmostEqual(corrected.electron_energy_J, corrected_initial_energy, delta=1.0e-28)
            legacy._implicit_electron_diffusion(1.0)
            self.assertGreater(float(legacy.electron.temperature_K[1, 0, 0]), 300.0)

    def test_exact_face_contact_operator_conserves_electron_energy(self):
        with tempfile.TemporaryDirectory() as directory:
            geometry = self.make_geometry(np.ones((5, 1, 1)), directory)
            spacing_A = np.full(3, 10.0)
            size_A = np.asarray([50.0, 10.0, 10.0])
            face_file = Path(directory) / "contacts.npz"
            np.savez_compressed(
                face_file,
                contact_x=np.ones((5, 1, 1)),
                contact_y=np.zeros((5, 1, 1)),
                contact_z=np.zeros((5, 1, 1)),
                ag_fraction=np.ones((5, 1, 1)),
                size_A=size_A,
                spacing_A=spacing_A,
            )
            initial_temperature = np.asarray([1200.0, 300.0, 300.0, 300.0, 300.0]).reshape(5, 1, 1)
            model = QSGSDoubleTemperatureGrid(
                geometry,
                "A",
                peak_flux_eV_ps_A2=0.0,
                initial_electron_temperature_K=initial_temperature,
                face_contacts_npz=face_file,
            )
            initial_energy = model.electron_energy_J
            record = model.advance(1.0e-4, heat_on=False, coupling_on=False)
            self.assertLess(float(model.electron.temperature_K[0, 0, 0]), 1200.0)
            self.assertGreater(float(model.electron.temperature_K[1, 0, 0]), 300.0)
            self.assertAlmostEqual(model.electron_energy_J, initial_energy, delta=1.0e-28)
            self.assertLess(abs(record["energy_balance_residual_J"]), 1.0e-28)

    def test_sampled_run_aggregates_energy_over_the_reported_interval(self):
        with tempfile.TemporaryDirectory() as directory:
            geometry = self.make_geometry(np.ones((1, 1, 1)), directory)
            model = QSGSDoubleTemperatureGrid(
                geometry,
                "A",
                initial_electron_temperature_K=4000.0,
            )
            rows = model.run(0.001, 5.0e-4, heat_on=False, sample_interval_ps=0.001)
            self.assertEqual(len(rows), 1)
            self.assertAlmostEqual(rows[0]["dt_ps"], 0.001, places=14)
            self.assertGreaterEqual(rows[0]["electron_substeps"], 1.0)
            self.assertGreater(rows[0]["electron_to_lattice_energy_J"], 0.0)
            self.assertLess(abs(rows[0]["energy_balance_residual_J"]), 1.0e-29)

    def test_qsgs_source_integral_and_coupled_energy_ledger(self):
        geometry = PHASE_DIR / "qsgs_phase_fractions_10A.npz"
        for boundary in ("A", "B"):
            model = QSGSDoubleTemperatureGrid(geometry, boundary)
            self.assertAlmostEqual(
                model.actual_source_power_W,
                model.lattice.actual_source_power_W,
                delta=model.lattice.actual_source_power_W * 2.0e-15,
            )
            initial_total = model.total_energy_J
            record = model.advance(5.0e-4, heat_on=True)
            self.assertGreater(record["source_energy_J"], 0.0)
            self.assertGreater(record["electron_to_lattice_energy_J"], 0.0)
            self.assertGreater(model.electron.temperature_K.max(), 300.0)
            self.assertGreater(model.lattice.state.temperature_K.max(), 300.0)
            self.assertAlmostEqual(
                model.total_energy_J - initial_total,
                record["source_energy_J"],
                delta=max(1.0e-29, record["source_energy_J"] * 2.0e-10),
            )
            self.assertLess(abs(record["energy_balance_residual_J"]), 1.0e-28)
            self.assertAlmostEqual(
                float(model.coupling_energy_density_J_m3.sum() * model.lattice.cell_volume_m3),
                model.coupling_energy_J,
                delta=max(1.0e-31, model.coupling_energy_J * 1.0e-12),
            )
            self.assertAlmostEqual(
                model.actual_source_power_W,
                model.electron.source_power.sum() * EV_J / 1.0e-12,
                places=15,
            )


if __name__ == "__main__":
    unittest.main()
