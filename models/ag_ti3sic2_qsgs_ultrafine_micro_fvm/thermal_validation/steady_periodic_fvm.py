#!/usr/bin/env python3
"""周期边界下的稳态有限体积导热与RVE均匀化。

相场为二值体素；面热阻由两侧半单元热阻和可选Ag/骨架界面热阻串联得到。
使用矩阵无显式存储的PCG，并以常系数周期拉普拉斯的FFT逆作为预条件器。
"""

from __future__ import annotations

from dataclasses import dataclass
import time

import numpy as np
from scipy import fft


@dataclass
class CGResult:
    solution: np.ndarray
    converged: bool
    iterations: int
    relative_residual: float
    residual_history: list[float]
    wall_seconds: float


class PeriodicConductivityFVM:
    """Cell-centred periodic corrector solver for a scalar heterogeneous k field."""

    def __init__(
        self,
        skeleton_mask: np.ndarray,
        spacing_m: tuple[float, float, float] | list[float] | np.ndarray,
        k_ag_W_m_K: float,
        k_skeleton_W_m_K: float,
        interface_resistance_m2_K_W: float = 0.0,
        fft_workers: int = 1,
    ) -> None:
        mask = np.asarray(skeleton_mask, dtype=bool)
        if mask.ndim != 3 or min(mask.shape) < 2:
            raise ValueError("skeleton_mask must be a 3D array with each dimension >= 2")
        spacing = np.asarray(spacing_m, dtype=np.float64)
        if spacing.shape != (3,) or np.any(~np.isfinite(spacing)) or np.any(spacing <= 0.0):
            raise ValueError("spacing_m must contain three finite positive lengths")
        if not np.isfinite(k_ag_W_m_K) or not np.isfinite(k_skeleton_W_m_K):
            raise ValueError("thermal conductivities must be finite")
        if min(k_ag_W_m_K, k_skeleton_W_m_K) <= 0.0:
            raise ValueError("thermal conductivities must be positive")
        if not np.isfinite(interface_resistance_m2_K_W) or interface_resistance_m2_K_W < 0.0:
            raise ValueError("interface resistance must be finite and nonnegative")
        self.mask = mask
        self.shape = mask.shape
        self.spacing_m = spacing
        self.k_ag = float(k_ag_W_m_K)
        self.k_skeleton = float(k_skeleton_W_m_K)
        self.interface_resistance = float(interface_resistance_m2_K_W)
        self.fft_workers = max(1, int(fft_workers))
        self.k_cell = np.where(mask, self.k_skeleton, self.k_ag).astype(np.float64, copy=False)
        self.cell_volume_m3 = float(np.prod(spacing))
        self.domain_size_m = spacing * np.asarray(mask.shape, dtype=np.float64)

        # C_axis is the conductance of each oriented cell face, in W/K.
        self.face_conductance_W_K: list[np.ndarray] = []
        self.interface_faces: list[np.ndarray] = []
        self.face_area_m2: list[float] = []
        for axis in range(3):
            neighbor_k = np.roll(self.k_cell, -1, axis=axis)
            interface = mask != np.roll(mask, -1, axis=axis)
            resistance = 0.5 * self.spacing_m[axis] * (1.0 / self.k_cell + 1.0 / neighbor_k)
            if self.interface_resistance:
                resistance += interface * self.interface_resistance
            face_area = self.cell_volume_m3 / self.spacing_m[axis]
            conductance = face_area / resistance
            self.face_conductance_W_K.append(np.asarray(conductance, dtype=np.float64))
            self.interface_faces.append(interface)
            self.face_area_m2.append(float(face_area))
            del neighbor_k, resistance, conductance

        self._fft_eigenvalues = self._build_fft_preconditioner_spectrum()

    def _build_fft_preconditioner_spectrum(self) -> np.ndarray:
        """Inverse spectrum for a constant-coefficient periodic graph Laplacian."""
        nx, ny, nz = self.shape
        means = [float(np.mean(face)) for face in self.face_conductance_W_K]
        fx = np.arange(nx, dtype=np.float64) / nx
        fy = np.arange(ny, dtype=np.float64) / ny
        fz = np.arange(nz // 2 + 1, dtype=np.float64) / nz
        lx = 2.0 * means[0] * (1.0 - np.cos(2.0 * np.pi * fx))
        ly = 2.0 * means[1] * (1.0 - np.cos(2.0 * np.pi * fy))
        lz = 2.0 * means[2] * (1.0 - np.cos(2.0 * np.pi * fz))
        eigenvalues = lx[:, None, None] + ly[None, :, None] + lz[None, None, :]
        eigenvalues[0, 0, 0] = 1.0  # zero-mean subspace; constant mode is projected out
        return eigenvalues

    def apply_operator(self, field: np.ndarray) -> np.ndarray:
        """Apply the conservative periodic FV stiffness matrix (W/K)."""
        x = np.asarray(field, dtype=np.float64)
        if x.shape != self.shape:
            raise ValueError("field shape must match the geometry")
        result = np.zeros_like(x)
        for axis, conductance in enumerate(self.face_conductance_W_K):
            forward_difference = np.roll(x, -1, axis=axis) - x
            forward_flux = conductance * forward_difference
            result += np.roll(forward_flux, 1, axis=axis)
            result -= forward_flux
        return result

    def _apply_preconditioner(self, residual: np.ndarray) -> np.ndarray:
        spectrum = fft.rfftn(residual, workers=self.fft_workers)
        spectrum /= self._fft_eigenvalues
        result = fft.irfftn(spectrum, s=self.shape, workers=self.fft_workers)
        result -= np.mean(result)
        return result

    def solve_cg(
        self,
        rhs: np.ndarray,
        relative_tolerance: float = 1.0e-9,
        maximum_iterations: int = 1000,
    ) -> CGResult:
        """Solve A x=b on the mean-zero subspace with FFT-preconditioned CG."""
        b = np.asarray(rhs, dtype=np.float64)
        if b.shape != self.shape:
            raise ValueError("rhs shape must match the geometry")
        if relative_tolerance <= 0.0 or maximum_iterations <= 0:
            raise ValueError("solver tolerance and iteration limit must be positive")
        b = b - np.mean(b)
        norm_b = float(np.linalg.norm(b.ravel()))
        x = np.zeros(self.shape, dtype=np.float64)
        if norm_b == 0.0:
            return CGResult(x, True, 0, 0.0, [0.0], 0.0)

        start = time.perf_counter()
        residual = b.copy()
        z = self._apply_preconditioner(residual)
        direction = z.copy()
        rz_old = float(np.vdot(residual, z).real)
        history = [1.0]
        converged = False
        relative_residual = 1.0
        iteration = 0

        for iteration in range(1, maximum_iterations + 1):
            applied = self.apply_operator(direction)
            denominator = float(np.vdot(direction, applied).real)
            if not np.isfinite(denominator) or denominator <= 0.0:
                raise RuntimeError("CG encountered a nonpositive curvature direction")
            alpha = rz_old / denominator
            x += alpha * direction
            residual -= alpha * applied
            residual -= np.mean(residual)
            relative_residual = float(np.linalg.norm(residual.ravel()) / norm_b)
            history.append(relative_residual)
            if relative_residual <= relative_tolerance:
                converged = True
                break
            z = self._apply_preconditioner(residual)
            rz_new = float(np.vdot(residual, z).real)
            if not np.isfinite(rz_new) or rz_new <= 0.0:
                raise RuntimeError("CG preconditioner produced a nonpositive residual product")
            beta = rz_new / rz_old
            direction *= beta
            direction += z
            rz_old = rz_new

        wall_seconds = time.perf_counter() - start
        if not converged:
            raise RuntimeError(
                f"CG did not converge in {maximum_iterations} iterations; "
                f"relative residual={relative_residual:.3e}"
            )
        return CGResult(x, True, iteration, relative_residual, history, wall_seconds)

    def solve_effective_tensor(
        self,
        relative_tolerance: float = 1.0e-9,
        maximum_iterations: int = 1000,
        snapshot_axis: int | str | None = None,
        snapshot_delta_temperature_K: float = 10.0,
        cold_temperature_K: float = 300.0,
    ) -> dict[str, object]:
        """Solve three periodic corrector problems and form the full K tensor."""
        if snapshot_axis not in (None, 0, 1, 2, "all"):
            raise ValueError("snapshot_axis must be None, 0, 1, 2, or 'all'")
        tensor = np.zeros((3, 3), dtype=np.float64)
        solves: list[CGResult] = []
        snapshots: dict[int, dict[str, object]] = {}
        contributions: list[dict[str, float]] = []
        interface_flux_statistics: list[dict[str, float | int]] = []

        for imposed_axis in range(3):
            axis = imposed_axis
            spacing = self.spacing_m[axis]
            face = self.face_conductance_W_K[axis]
            rhs = spacing * (face - np.roll(face, 1, axis=axis))
            result = self.solve_cg(rhs, relative_tolerance, maximum_iterations)
            solves.append(result)
            del rhs

            column_flux = []
            for flux_axis in range(3):
                chi_difference = np.roll(result.solution, -1, axis=flux_axis) - result.solution
                temperature_difference = chi_difference.copy()
                if flux_axis == imposed_axis:
                    temperature_difference += self.spacing_m[flux_axis]
                q_face = -(self.face_conductance_W_K[flux_axis] / self.face_area_m2[flux_axis]) * (
                    temperature_difference
                )
                average_flux = float(np.mean(q_face))
                tensor[flux_axis, imposed_axis] = -average_flux
                column_flux.append(q_face)

                phase_right = np.roll(self.mask, -1, axis=flux_axis)
                same_ag = (~self.mask) & (~phase_right)
                same_skeleton = self.mask & phase_right
                interface = self.interface_faces[flux_axis]
                scale = -1.0
                contributions.append({
                    "imposed_axis": imposed_axis,
                    "flux_axis": flux_axis,
                    "ag_W_m_K": scale * float(np.mean(q_face * same_ag)),
                    "skeleton_W_m_K": scale * float(np.mean(q_face * same_skeleton)),
                    "interface_W_m_K": scale * float(np.mean(q_face * interface)),
                })
                gradient_for_flux = snapshot_delta_temperature_K / self.domain_size_m[imposed_axis]
                interface_q = np.abs(q_face[interface] * gradient_for_flux)
                interface_flux_statistics.append({
                    "imposed_axis": imposed_axis,
                    "interface_face_normal_axis": flux_axis,
                    "face_count": int(interface_q.size),
                    "mean_abs_W_m2": float(interface_q.mean()) if interface_q.size else 0.0,
                    "p95_abs_W_m2": float(np.percentile(interface_q, 95)) if interface_q.size else 0.0,
                    "maximum_abs_W_m2": float(interface_q.max()) if interface_q.size else 0.0,
                })
                del interface_q
                del phase_right, same_ag, same_skeleton, chi_difference, temperature_difference, q_face

            if snapshot_axis == "all" or imposed_axis == snapshot_axis:
                interface = self.interface_faces[imposed_axis]
                normal_flux = column_flux[imposed_axis]
                gradient = snapshot_delta_temperature_K / self.domain_size_m[imposed_axis]
                abs_interface_flux = np.abs(normal_flux[interface] * gradient)
                if abs_interface_flux.size:
                    interface_stats = {
                        "interface_face_count": int(abs_interface_flux.size),
                        "mean_abs_W_m2": float(abs_interface_flux.mean()),
                        "p95_abs_W_m2": float(np.percentile(abs_interface_flux, 95)),
                        "maximum_abs_W_m2": float(abs_interface_flux.max()),
                    }
                else:
                    interface_stats = {
                        "interface_face_count": 0,
                        "mean_abs_W_m2": 0.0,
                        "p95_abs_W_m2": 0.0,
                        "maximum_abs_W_m2": 0.0,
                    }
                coordinate_shape = [1, 1, 1]
                coordinate_shape[imposed_axis] = self.shape[imposed_axis]
                coordinate = (
                    (np.arange(self.shape[imposed_axis], dtype=np.float32) + 0.5)
                    * self.spacing_m[imposed_axis]
                ).reshape(coordinate_shape)
                temperature = cold_temperature_K + gradient * (
                    coordinate + result.solution
                )
                one_snapshot: dict[str, object] = {
                    "temperature_K": np.asarray(temperature, dtype=np.float32),
                    "material_id": self.mask.astype(np.uint8),
                    "interface_stats": interface_stats,
                    "corrector_K": np.asarray(result.solution * gradient, dtype=np.float32),
                }
                for component_axis, component_face in enumerate(column_flux):
                    q_component = 0.5 * (
                        component_face + np.roll(component_face, 1, axis=component_axis)
                    ) * gradient
                    one_snapshot[f"heat_flux_{component_axis}_W_m2"] = np.asarray(
                        q_component, dtype=np.float32
                    )
                snapshots[imposed_axis] = one_snapshot
                del coordinate, temperature

        return {
            "effective_conductivity_tensor_W_m_K": tensor,
            "solver": [
                {
                    "imposed_axis": axis,
                    "iterations": solve.iterations,
                    "relative_residual": solve.relative_residual,
                    "wall_seconds": solve.wall_seconds,
                    "residual_history": solve.residual_history,
                }
                for axis, solve in enumerate(solves)
            ],
            "phase_flux_contributions": contributions,
            "interface_heat_flux_statistics": interface_flux_statistics,
            "snapshots_by_imposed_axis": snapshots,
        }
