"""Analytic constant-property checks for the implicit electron diffusion step."""

from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np


MODEL_DIR = Path(__file__).resolve().parent
PHASE_DIR = MODEL_DIR / "phase_field_qsgs"
sys.path.insert(0, str(PHASE_DIR))
from qsgs_ttm_enthalpy import QSGSDoubleTemperatureGrid  # noqa: E402
from qsgs_component_coarsegrain import build_component_graph, component_effective_conductance  # noqa: E402
from qsgs_face_contacts import exact_face_contact_fractions  # noqa: E402
from qsgs_transport_benchmark import (  # noqa: E402
    ROI_CASES,
    SOURCE as QSGS_SOURCE,
    effective_ag_conductance,
)
from qsgs_component_transient_benchmark import aggregate_component_temperature  # noqa: E402
from qsgs_component_ttm_benchmark import (  # noqa: E402
    ELECTRON_TABLE,
    _implicit_local_coupling,
    _lattice_h_vol_eV_A3,
)
from ag_ttm_reference import AgParameters  # noqa: E402
from thermophysical_model import AgEnthalpy  # noqa: E402


class ConstantElectronProperties:
    """Linear Ue(T) and constant ke for analytic diffusion modes."""

    def __init__(self, ce: float = 1.0, ke: float = 1.0) -> None:
        self.temperature_K = np.asarray([300.0, 20_000.0])
        self.ue_eV_A3 = ce * (self.temperature_K - 300.0)
        self.ce_value = float(ce)
        self.ke_value = float(ke)
        self.Tmin = 300.0
        self.Tmax = 20_000.0

    def ue(self, temperature_K):
        return self.ce_value * (np.asarray(temperature_K, dtype=float) - 300.0)

    def ke(self, temperature_K):
        return np.full_like(np.asarray(temperature_K, dtype=float), self.ke_value)

    def _checked_temperature(self, temperature_K):
        value = np.asarray(temperature_K, dtype=float)
        if np.any(value < self.Tmin) or np.any(value > self.Tmax):
            raise ValueError("benchmark temperature outside its test range")
        return value


def build_periodic_1d_model(cell_count: int, length_A: float, mode: int = 1):
    dx = length_A / cell_count
    with tempfile.TemporaryDirectory() as directory:
        geometry = Path(directory) / "pure_ag_geometry.npz"
        np.savez_compressed(
            geometry,
            ag_fraction=np.ones((cell_count, 1, 1)),
            skeleton_fraction=np.zeros((cell_count, 1, 1)),
            size_A=np.asarray([length_A, 10.0, 10.0]),
            spacing_A=np.asarray([dx, 10.0, 10.0]),
        )
        model = QSGSDoubleTemperatureGrid(
            geometry,
            "A",
            peak_flux_eV_ps_A2=0.0,
            initial_electron_temperature_K=300.0,
        )
    params = ConstantElectronProperties(ce=1.0, ke=1.0)
    model.parameters = params
    model.electron.parameters = params
    x = (np.arange(cell_count, dtype=float) + 0.5) * dx
    # Stay well inside the mocked table while retaining a clean Fourier mode.
    initial = 1_000.0 + 100.0 * np.sin(2.0 * np.pi * mode * x / length_A)
    model.electron.temperature_K[:, 0, 0] = initial
    return model, initial


def fourier_amplitude(values: np.ndarray, length_A: float, mode: int) -> float:
    values = np.asarray(values, dtype=float)
    phase = 2.0 * np.pi * mode * (np.arange(values.size) + 0.5) / values.size
    centered = values - float(values.mean())
    sine = 2.0 * float(np.dot(centered, np.sin(phase))) / values.size
    cosine = 2.0 * float(np.dot(centered, np.cos(phase))) / values.size
    return float(np.hypot(sine, cosine))


class ElectronDiffusionBenchmarkTests(unittest.TestCase):
    def test_component_initialization_is_volume_weighted_from_native_voxels(self):
        ag = np.asarray([[[True, True, False]], [[True, False, False]]])
        component_id = np.asarray([[[0, 0, -1]], [[0, -1, -1]]])
        temperature = np.asarray([[[300.0, 600.0, 0.0]], [[900.0, 0.0, 0.0]]])
        aggregated = aggregate_component_temperature(component_id, ag, temperature, 1)
        np.testing.assert_allclose(aggregated, [600.0], rtol=0.0, atol=1.0e-12)

    def test_component_ttm_local_coupling_conserves_electron_plus_lattice_energy(self):
        parameters = AgParameters.from_csv(ELECTRON_TABLE, allow_low_temperature_extension=False)
        lattice = AgEnthalpy()
        volume = np.asarray([120.0, 75.0])
        te0 = np.asarray([1500.0, 900.0])
        tl0 = np.asarray([300.0, 300.0])
        fl0 = np.zeros(2)
        total_before = float(
            np.dot(volume, parameters.ue(te0))
            + np.dot(volume, _lattice_h_vol_eV_A3(tl0, fl0, lattice))
        )
        te1, tl1, fl1, transfer = _implicit_local_coupling(
            volume, te0, tl0, fl0, 0.002, parameters, lattice
        )
        total_after = float(
            np.dot(volume, parameters.ue(te1))
            + np.dot(volume, _lattice_h_vol_eV_A3(tl1, fl1, lattice))
        )
        self.assertGreater(transfer, 0.0)
        self.assertTrue(np.all(te1 < te0))
        self.assertTrue(np.all(tl1 > tl0))
        self.assertTrue(np.all(fl1 == 0.0))
        self.assertAlmostEqual(total_after, total_before, delta=1.0e-10)

    def test_component_graph_preserves_uniform_conductivity_and_ag_capacity(self):
        ag = np.ones((4, 3, 2), dtype=bool)
        size = np.asarray([8.0, 9.0, 10.0])
        graph = build_component_graph(ag, size, coarse_voxels=(3, 3, 3))
        self.assertEqual(int(graph.voxel_count.sum()), int(ag.sum()))
        for axis in range(3):
            result = component_effective_conductance(graph, axis)
            self.assertAlmostEqual(result["effective_k_over_intrinsic_k"], 1.0, delta=1.0e-12)
            self.assertLess(result["hot_cold_flux_relative_balance"], 1.0e-12)

    def test_component_graph_does_not_bridge_disconnected_ag_islands(self):
        ag = np.zeros((3, 3, 3), dtype=bool)
        ag[0, 1, 1] = True
        ag[2, 1, 1] = True
        graph = build_component_graph(ag, np.asarray([9.0, 9.0, 9.0]), (3, 3, 3))
        self.assertEqual(graph.voxel_count.size, 2)
        self.assertEqual(graph.edge_conductance.size, 0)
        result = component_effective_conductance(graph, 0)
        self.assertEqual(result["effective_k_over_intrinsic_k"], 0.0)

    def test_qsgs_component_graph_matches_voxel_conductance_in_all_sampled_axes(self):
        with np.load(QSGS_SOURCE, allow_pickle=False) as data:
            full_ag = ~np.asarray(data["skeleton"], dtype=bool)
            full_size_A = np.asarray(data["size_A"], dtype=float)
        checked_connected = 0
        checked_blocked = 0
        for case in ROI_CASES:
            start, shape = case["start"], case["shape"]
            ag = full_ag[tuple(slice(s, s + n) for s, n in zip(start, shape))]
            size = full_size_A * np.asarray(shape) / np.asarray(full_ag.shape)
            spacing = size / np.asarray(shape)
            graph = build_component_graph(ag, size, coarse_voxels=(2, 1, 1))
            contacts = exact_face_contact_fractions(ag, ag.astype(float), size, spacing)
            for axis in range(3):
                reference = effective_ag_conductance(
                    ag, size, ag.astype(float), spacing, contacts, axis
                )
                coarse = component_effective_conductance(graph, axis)
                self.assertLess(coarse["hot_cold_flux_relative_balance"], 1.0e-8)
                if reference["effective_k_over_intrinsic_k"] > 1.0e-10:
                    relative_error = abs(
                        coarse["effective_k_over_intrinsic_k"]
                        - reference["effective_k_over_intrinsic_k"]
                    ) / reference["effective_k_over_intrinsic_k"]
                    self.assertLessEqual(relative_error, 0.05)
                    checked_connected += 1
                else:
                    self.assertEqual(coarse["effective_k_over_intrinsic_k"], 0.0)
                    checked_blocked += 1
        self.assertEqual(checked_connected, 21)
        self.assertEqual(checked_blocked, 6)

    def test_geometry_transport_reference_recovers_uniform_cube_conductivity(self):
        ag = np.ones((4, 3, 2), dtype=bool)
        size = np.asarray([8.0, 9.0, 10.0])
        spacing = size / np.asarray(ag.shape)
        contacts = exact_face_contact_fractions(ag, ag.astype(float), size, spacing)
        result = effective_ag_conductance(ag, size, ag.astype(float), spacing, contacts)
        self.assertAlmostEqual(result["effective_k_over_intrinsic_k"], 1.0, delta=1.0e-12)
        self.assertLess(result["hot_cold_flux_relative_balance"], 1.0e-12)

    def test_single_backward_euler_step_matches_straight_channel_fourier_mode(self):
        count, length, mode, dt = 32, 64.0, 4, 0.05
        model, initial = build_periodic_1d_model(count, length, mode)
        model.electron_face_contact_fractions = (
            np.ones_like(model.lattice.ag_fraction),
            np.zeros_like(model.lattice.ag_fraction),
            np.zeros_like(model.lattice.ag_fraction),
        )
        initial_mean = float(initial.mean())
        initial_amplitude = fourier_amplitude(initial, length, mode)
        dx = length / count
        discrete_eigenvalue = 4.0 * np.sin(np.pi * mode / count) ** 2 / dx**2
        expected_amplitude = initial_amplitude / (1.0 + dt * discrete_eigenvalue)
        energy_before = model.electron_energy_J

        diffusion_change = model._implicit_electron_diffusion(dt)
        final = model.electron.temperature_K[:, 0, 0]
        final_amplitude = fourier_amplitude(final, length, mode)

        self.assertAlmostEqual(final_amplitude, expected_amplitude, delta=2.0e-7)
        self.assertAlmostEqual(float(final.mean()), initial_mean, delta=2.0e-9)
        self.assertAlmostEqual(model.electron_energy_J, energy_before, delta=1.0e-28)
        self.assertLess(abs(diffusion_change), 1.0e-28)

    def test_spatial_refinement_follows_second_order_periodic_laplacian(self):
        length, mode, dt, nsteps = 64.0, 4, 0.002, 100
        solver_errors = []
        continuum_spatial_errors = []
        for count in (16, 32, 64):
            model, initial = build_periodic_1d_model(count, length, mode)
            initial_amplitude = fourier_amplitude(initial, length, mode)
            dx = length / count
            discrete_eigenvalue = 4.0 * np.sin(np.pi * mode / count) ** 2 / dx**2
            discrete_be = initial_amplitude * (1.0 + dt * discrete_eigenvalue) ** (-nsteps)
            for _ in range(nsteps):
                model._implicit_electron_diffusion(dt)
            final = model.electron.temperature_K[:, 0, 0]
            actual = fourier_amplitude(final, length, mode)
            solver_errors.append(abs(actual - discrete_be) / discrete_be)
            continuum_eigenvalue = (2.0 * np.pi * mode / length) ** 2
            semidiscrete_exact = initial_amplitude * np.exp(-discrete_eigenvalue * dt * nsteps)
            continuum_exact = initial_amplitude * np.exp(-continuum_eigenvalue * dt * nsteps)
            continuum_spatial_errors.append(abs(semidiscrete_exact - continuum_exact))
        # The numerical solution follows the exact semidiscrete BE mode, while
        # the semidiscrete Laplacian error against the continuum scales as dx².
        self.assertLess(max(solver_errors), 2.0e-7)
        self.assertGreater(continuum_spatial_errors[0] / continuum_spatial_errors[1], 3.5)
        self.assertGreater(continuum_spatial_errors[1] / continuum_spatial_errors[2], 3.5)

    def test_backward_euler_time_refinement_is_first_order(self):
        count, length, mode, final_time = 32, 64.0, 4, 0.5
        dx = length / count
        discrete_eigenvalue = 4.0 * np.sin(np.pi * mode / count) ** 2 / dx**2
        errors = []
        for dt in (0.05, 0.025, 0.0125):
            model, initial = build_periodic_1d_model(count, length, mode)
            amplitude0 = fourier_amplitude(initial, length, mode)
            nsteps = round(final_time / dt)
            for _ in range(nsteps):
                model._implicit_electron_diffusion(dt)
            final = model.electron.temperature_K[:, 0, 0]
            actual = fourier_amplitude(final, length, mode)
            # Compare to the exact-in-time decay of this grid's discrete mode,
            # isolating backward-Euler temporal error from spatial error.
            semidiscrete_exact = amplitude0 * np.exp(-discrete_eigenvalue * final_time)
            errors.append(abs(actual - semidiscrete_exact))
        self.assertGreater(errors[0], errors[1])
        self.assertGreater(errors[1], errors[2])
        self.assertGreater(errors[0] / errors[1], 1.5)
        self.assertGreater(errors[1] / errors[2], 1.5)


if __name__ == "__main__":
    unittest.main()
