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
from qsgs_face_contacts import exact_face_contact_fractions  # noqa: E402
from qsgs_transport_benchmark import effective_ag_conductance  # noqa: E402


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
