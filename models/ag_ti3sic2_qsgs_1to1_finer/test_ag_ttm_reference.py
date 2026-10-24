import unittest
from pathlib import Path

import numpy as np

from ag_ttm_reference import (
    AgParameters,
    AgTTMGrid,
    AtomGridMap,
    TTMRangeError,
    build_ag_grid_fields,
    cloud_in_cell_atom_map,
    lattice_temperature_by_cell,
    rescale_lattice_kinetic,
)


HERE = Path(__file__).resolve().parent
PARAMETERS = AgParameters.from_csv(HERE / "ttm_data" / "Ag_TTM_hybrid_figure_fit_300_20000K.csv")


class ParameterTableTests(unittest.TestCase):
    def test_candidate_table_range_and_energy_monotonicity(self):
        self.assertEqual(PARAMETERS.Tmin, 300.0)
        self.assertEqual(PARAMETERS.Tmax, 20000.0)
        self.assertTrue(np.all(np.diff(PARAMETERS.ue_eV_A3) > 0.0))
        self.assertGreater(float(PARAMETERS.ce(300.0)), 0.0)
        self.assertGreater(float(PARAMETERS.g(300.0)), 0.0)
        self.assertGreater(float(PARAMETERS.ke(300.0)), 0.0)

    def test_no_temperature_extrapolation(self):
        with self.assertRaises(TTMRangeError):
            PARAMETERS.ce(299.99)
        with self.assertRaises(TTMRangeError):
            PARAMETERS.ke(20000.01)

    def test_explicit_low_temperature_closure_is_opt_in_and_invertible(self):
        extended = AgParameters.from_csv(
            HERE / "ttm_data" / "Ag_TTM_hybrid_figure_fit_300_20000K.csv",
            allow_low_temperature_extension=True,
        )
        temperatures = np.array([1.0, 50.0, 150.0, 299.0, 300.0])
        self.assertEqual(extended.Tmin, 1.0)
        self.assertTrue(np.all(extended.ce(temperatures) > 0.0))
        self.assertTrue(np.all(np.diff(extended.ue(temperatures)) > 0.0))
        np.testing.assert_allclose(
            extended.temperature_from_ue(extended.ue(temperatures)),
            temperatures,
            atol=1.0e-8,
        )
        with self.assertRaises(TTMRangeError):
            extended.ce(0.99)

    def test_energy_temperature_round_trip(self):
        t = np.array([300.0, 500.0, 2000.0, 10000.0, 20000.0])
        np.testing.assert_allclose(PARAMETERS.temperature_from_ue(PARAMETERS.ue(t)), t, atol=1e-9)


class GridEnergyTests(unittest.TestCase):
    def test_source_energy_is_exact_and_can_be_switched_off(self):
        grid = AgTTMGrid(
            PARAMETERS,
            fraction_ag=np.ones((1, 1, 1)),
            spacing_A=(10.0, 10.0, 10.0),
            source_power_eV_ps=np.array([[[0.1]]]),
        )
        result = grid.advance(1.0e-4, lattice_temperature_K=300.0, coupling_on=False)
        self.assertAlmostEqual(result.source_to_electron_eV, 1.0e-5, places=13)
        self.assertAlmostEqual(result.conservation_residual_eV, 0.0, delta=1.0e-12)
        temp_after_heat = float(grid.temperature_K[0, 0, 0])
        self.assertGreater(temp_after_heat, 300.0)

        result_off = grid.advance(1.0e-4, lattice_temperature_K=300.0, heat_on=False, coupling_on=False)
        self.assertEqual(result_off.source_to_electron_eV, 0.0)
        self.assertEqual(result_off.substeps, 1)

    def test_electron_lattice_exchange_is_equal_and_opposite(self):
        grid = AgTTMGrid(
            PARAMETERS,
            fraction_ag=np.ones((1, 1, 1)),
            spacing_A=(10.0, 10.0, 10.0),
            initial_temperature_K=1000.0,
        )
        before = float(grid.electron_energy_eV.sum())
        result = grid.advance(1.0e-4, lattice_temperature_K=300.0, heat_on=False)
        after = float(grid.electron_energy_eV.sum())
        self.assertGreater(result.electron_to_lattice_eV, 0.0)
        self.assertAlmostEqual(after - before, -result.electron_to_lattice_eV, delta=1.0e-12)
        self.assertAlmostEqual(result.conservation_residual_eV, 0.0, delta=1.0e-12)

    def test_sub_300k_lattice_fluctuation_uses_explicit_low_t_model(self):
        extended = AgParameters.from_csv(
            HERE / "ttm_data" / "Ag_TTM_hybrid_figure_fit_300_20000K.csv",
            allow_low_temperature_extension=True,
        )
        grid = AgTTMGrid(extended, np.ones((1, 1, 1)), (10.0, 10.0, 10.0))
        before = float(grid.electron_energy_eV.sum())
        result = grid.advance(1.0e-4, lattice_temperature_K=250.0, heat_on=False)
        after = float(grid.electron_energy_eV.sum())
        self.assertGreater(result.electron_to_lattice_eV, 0.0)
        self.assertLess(float(grid.temperature_K[0, 0, 0]), 300.0)
        self.assertAlmostEqual(after - before, -result.electron_to_lattice_eV, delta=1.0e-12)
        self.assertAlmostEqual(result.conservation_residual_eV, 0.0, delta=1.0e-12)

    def test_diffusion_conserves_electron_energy_and_respects_zero_flux_topology(self):
        grid = AgTTMGrid(
            PARAMETERS,
            fraction_ag=np.ones((1, 1, 3)),
            spacing_A=(10.0, 10.0, 10.0),
        )
        grid.temperature_K[0, 0, :] = [1000.0, 300.0, 300.0]
        before = float(grid.electron_energy_eV.sum())
        result = grid.advance(1.0e-5, lattice_temperature_K=300.0, heat_on=False, coupling_on=False)
        after = float(grid.electron_energy_eV.sum())
        self.assertAlmostEqual(after, before, delta=1.0e-12)
        self.assertAlmostEqual(result.diffusion_net_eV, 0.0, delta=1.0e-12)
        self.assertGreater(grid.temperature_K[0, 0, 1], 300.0)
        self.assertGreater(grid.temperature_K[0, 0, 0], grid.temperature_K[0, 0, 1])

    def test_no_diffusion_through_empty_ag_cells(self):
        fraction = np.array([[[1.0, 0.0, 1.0]]])
        grid = AgTTMGrid(PARAMETERS, fraction, (10.0, 10.0, 10.0))
        grid.temperature_K[0, 0, 0] = 1500.0
        result = grid.advance(1.0e-5, lattice_temperature_K=300.0, heat_on=False, coupling_on=False)
        self.assertAlmostEqual(grid.temperature_K[0, 0, 2], 300.0, delta=1.0e-12)
        self.assertAlmostEqual(result.conservation_residual_eV, 0.0, delta=1.0e-12)

    def test_diffusion_is_periodic_xy(self):
        grid = AgTTMGrid(PARAMETERS, np.ones((2, 1, 1)), (10.0, 10.0, 10.0))
        grid.temperature_K[:, 0, 0] = [1500.0, 300.0]
        before = float(grid.electron_energy_eV.sum())
        grid.advance(1.0e-6, lattice_temperature_K=300.0, heat_on=False, coupling_on=False)
        after = float(grid.electron_energy_eV.sum())
        self.assertGreater(grid.temperature_K[1, 0, 0], 300.0)
        self.assertAlmostEqual(after, before, delta=1.0e-12)

    def test_fill_fraction_must_be_bounded(self):
        with self.assertRaises(ValueError):
            AgTTMGrid(PARAMETERS, np.array([[[1.001]]]), (10.0, 10.0, 10.0))


class LatticeTransferTests(unittest.TestCase):
    def test_cic_ag_fields_conserve_volume_and_source_power(self):
        xyz = np.array([[1.0, 1.0, 1.0], [19.0, 9.0, 19.0]])
        q = np.array([0.1, 0.2])
        mapping = build_ag_grid_fields(
            xyz,
            q,
            grid_shape=(2, 1, 2),
            spacing_A=(10.0, 10.0, 10.0),
            ag_atomic_volume_A3=17.0,
        )
        self.assertTrue(np.all(mapping.fraction_ag >= 0.0))
        self.assertTrue(np.all(mapping.fraction_ag <= 1.0))
        self.assertAlmostEqual(float(mapping.fraction_ag.sum() * 1000.0), 34.0, places=10)
        self.assertAlmostEqual(float(mapping.source_power_eV_ps.sum()), 5.1, places=10)
        self.assertTrue(np.allclose(mapping.atom_map.weights.sum(axis=1), 1.0))

    def test_atom_grid_map_wraps_xy_and_clamps_free_surface_z(self):
        mapping = cloud_in_cell_atom_map(
            np.array([[0.1, 9.9, 19.9]]), (2, 1, 2), (10.0, 10.0, 10.0)
        )
        self.assertTrue(np.all(mapping.cell_ids >= 0))
        self.assertTrue(np.all(mapping.cell_ids < 4))
        self.assertAlmostEqual(float(mapping.weights.sum()), 1.0, places=12)

    def test_cell_velocity_rescaling_matches_electron_energy_transfer(self):
        velocity = np.array(
            [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.1]]
        )
        masses = np.full(4, 107.8682)
        mapping = AtomGridMap.nearest_cell(np.array([0, 0, 1, 1]), (2, 1, 1))
        transfer = np.array([0.02, -0.00001])
        new_velocity, actual = rescale_lattice_kinetic(velocity, masses, mapping, transfer.reshape((2, 1, 1)))
        np.testing.assert_allclose(actual.ravel(), transfer, rtol=1.0e-10, atol=1.0e-12)
        self.assertTrue(np.all(np.isfinite(new_velocity)))

    def test_transfer_that_would_make_lattice_energy_negative_fails(self):
        velocity = np.zeros((1, 3))
        velocity[0, 0] = 0.01
        mapping = AtomGridMap.nearest_cell(np.array([0]), (1, 1, 1))
        with self.assertRaises(ValueError):
            rescale_lattice_kinetic(
                velocity, np.array([107.8682]), mapping, np.array([[[ -1.0 ]]])
            )

    def test_fractional_atom_to_grid_transfer_conserves_global_energy(self):
        velocity = np.array([[1.0, 0.2, 0.0], [-0.5, 0.0, 0.1]])
        masses = np.full(2, 107.8682)
        mapping = AtomGridMap(
            cell_ids=np.array([[0, 1], [0, 1]]),
            weights=np.array([[0.25, 0.75], [0.6, 0.4]]),
            grid_shape=(2, 1, 1),
        )
        transfer = np.array([0.01, -0.0001]).reshape((2, 1, 1))
        new_velocity, actual = rescale_lattice_kinetic(velocity, masses, mapping, transfer)
        self.assertAlmostEqual(float(actual.sum()), float(transfer.sum()), delta=1.0e-12)
        self.assertTrue(np.all(np.isfinite(new_velocity)))

    def test_fractional_map_gives_cell_lattice_temperatures(self):
        velocity = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        masses = np.full(2, 107.8682)
        mapping = AtomGridMap(
            cell_ids=np.array([[0, 1], [0, 1]]),
            weights=np.array([[0.5, 0.5], [0.25, 0.75]]),
            grid_shape=(2, 1, 1),
        )
        temperature = lattice_temperature_by_cell(velocity, masses, mapping)
        self.assertEqual(temperature.shape, (2, 1, 1))
        self.assertTrue(np.all(np.isfinite(temperature)))
        self.assertTrue(np.all(temperature > 0.0))


if __name__ == "__main__":
    unittest.main()
