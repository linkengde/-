from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np


MODEL_DIR = Path(__file__).resolve().parent
PHASE_DIR = MODEL_DIR / "phase_field_qsgs"
sys.path.insert(0, str(PHASE_DIR))
from qsgs_ttm_enthalpy import EV_J, QSGSDoubleTemperatureGrid  # noqa: E402


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
            record = model.advance(1.0e-4, heat_on=False, coupling_on=False)
            after = model.electron.temperature_K
            self.assertLess(after[0, 0, 0], before[0, 0, 0])
            self.assertGreater(after[1, 0, 0], before[1, 0, 0])
            self.assertAlmostEqual(model.electron_energy_J, initial_electron, delta=1.0e-29)
            self.assertLess(abs(record["energy_balance_residual_J"]), 1.0e-29)

    def test_one_dimensional_lattice_diffusion_is_conservative(self):
        with tempfile.TemporaryDirectory() as directory:
            geometry = self.make_geometry(np.ones((5, 1, 1)), directory)
            model = QSGSDoubleTemperatureGrid(
                geometry,
                "A",
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
            step = 0.5 * model.lattice_diffusion_timestep_ps
            record = model.advance(step, heat_on=False, coupling_on=False)
            final_temperature = model.lattice.state.temperature_K
            self.assertLess(final_temperature[0, 0, 0], initial_temperature[0, 0, 0])
            self.assertGreater(final_temperature[1, 0, 0], initial_temperature[1, 0, 0])
            self.assertAlmostEqual(model.lattice_energy_J, initial_enthalpy, delta=1.0e-28)
            self.assertLess(abs(record["lattice_diffusion_net_J"]), 1.0e-29)
            self.assertLess(abs(record["energy_balance_residual_J"]), 1.0e-29)

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
            self.assertGreater(rows[0]["electron_substeps"], 1.0)
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
                model.actual_source_power_W,
                model.electron.source_power.sum() * EV_J / 1.0e-12,
                places=15,
            )


if __name__ == "__main__":
    unittest.main()
