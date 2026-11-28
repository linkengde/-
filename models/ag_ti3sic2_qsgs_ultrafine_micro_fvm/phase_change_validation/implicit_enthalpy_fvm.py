#!/usr/bin/env python3
"""单温焓法有限体积瞬态导热器（矩阵自由隐式Newton-PCG）。"""

from __future__ import annotations

from dataclasses import dataclass, field
import time

import numpy as np
from scipy import fft
from scipy.sparse.linalg import LinearOperator, cg

from material_model import CandidateMaterials


AXIS_NAMES = "xyz"


@dataclass(frozen=True)
class BoundarySet:
    """周期、绝热（默认）和定温边界。例：{x_low: 300.0, z_high: 300.0}。"""

    periodic_axes: tuple[int, ...] = ()
    dirichlet_K: dict[str, float] = field(default_factory=dict)

    def validate(self, ndim: int) -> None:
        if any(axis < 0 or axis >= ndim for axis in self.periodic_axes):
            raise ValueError("periodic axis is outside the field dimensions")
        for key, value in self.dirichlet_K.items():
            if key not in {f"{AXIS_NAMES[a]}_{side}" for a in range(ndim) for side in ("low", "high")}:
                raise ValueError(f"invalid Dirichlet face name: {key}")
            if not np.isfinite(value):
                raise ValueError("Dirichlet temperature must be finite")
            axis = AXIS_NAMES.index(key[0])
            if axis in self.periodic_axes:
                raise ValueError("an axis cannot be both periodic and Dirichlet")


class ImplicitEnthalpyFVM:
    """常密度、二相储热/导热、单温Ag相变，细胞中心正交FVM。"""

    def __init__(
        self,
        ag_fraction: np.ndarray,
        spacing_m: tuple[float, ...] | list[float],
        materials: CandidateMaterials | None = None,
        interface_resistance_m2_K_W: float = 0.0,
        initial_temperature_K: float = 300.0,
        nonlinear_relative_tolerance: float = 1.0e-8,
        linear_relative_tolerance: float = 1.0e-8,
        maximum_nonlinear_iterations: int = 20,
        maximum_linear_iterations: int = 200,
    ) -> None:
        self.ag_fraction = np.asarray(ag_fraction, dtype=np.float64)
        if self.ag_fraction.ndim < 1 or self.ag_fraction.ndim > 3:
            raise ValueError("the finite-volume field must have one to three dimensions")
        if np.any(~np.isfinite(self.ag_fraction)) or np.any((self.ag_fraction < 0.0) | (self.ag_fraction > 1.0)):
            raise ValueError("phase fractions must be finite and lie in [0,1]")
        self.shape = self.ag_fraction.shape
        self.ndim = self.ag_fraction.ndim
        self.spacing_m = np.asarray(spacing_m, dtype=np.float64)
        if self.spacing_m.shape != (self.ndim,) or np.any(self.spacing_m <= 0.0):
            raise ValueError("spacing_m must contain one positive length per field axis")
        if interface_resistance_m2_K_W < 0.0 or not np.isfinite(interface_resistance_m2_K_W):
            raise ValueError("interface resistance must be finite and nonnegative")
        self.interface_resistance_m2_K_W = float(interface_resistance_m2_K_W)
        self.materials = materials or CandidateMaterials()
        self.volume_m3 = float(np.prod(self.spacing_m))
        self.nonlinear_rtol = float(nonlinear_relative_tolerance)
        self.linear_rtol = float(linear_relative_tolerance)
        self.max_newton = int(maximum_nonlinear_iterations)
        self.max_cg = int(maximum_linear_iterations)
        if initial_temperature_K > self.materials.ti_Tmax_K and np.any(self.ag_fraction < 1.0):
            raise ValueError("initial temperature exceeds the Ti3SiC2 candidate-property limit")
        self.temperature_K = np.full(self.shape, float(initial_temperature_K), dtype=np.float64)
        self.enthalpy_density_J_m3 = self.materials.enthalpy_density_J_m3(self.temperature_K, self.ag_fraction)
        self.initial_stored_enthalpy_J = float(np.sum(self.enthalpy_density_J_m3) * self.volume_m3)
        self.time_s = 0.0
        self.source_energy_J = 0.0
        self.boundary_energy_out_J = 0.0
        self.step_count = 0

    def state(self, temperature_K: np.ndarray | None = None) -> dict[str, np.ndarray]:
        temperature = self.temperature_K if temperature_K is None else np.asarray(temperature_K, dtype=float)
        return {
            "temperature_K": temperature,
            "ag_liquid_fraction": self.materials.liquid_fraction(temperature),
            "enthalpy_density_J_m3": self.materials.enthalpy_density_J_m3(temperature, self.ag_fraction),
            "apparent_capacity_J_m3_K": self.materials.volumetric_capacity_J_m3_K(temperature, self.ag_fraction),
        }

    def temperature_from_enthalpy_density(
        self,
        enthalpy_density_J_m3: np.ndarray,
        initial_temperature_K: np.ndarray | None = None,
        tolerance: float = 2.0e-12,
        maximum_iterations: int = 24,
    ) -> np.ndarray:
        """对窄mushy-zone焓曲线做单调Newton反算，不裁剪液相率。"""
        target = np.asarray(enthalpy_density_J_m3, dtype=np.float64)
        if target.shape != self.shape or np.any(~np.isfinite(target)):
            raise ValueError("enthalpy field must match the grid and contain only finite values")
        if initial_temperature_K is None:
            t = np.full(self.shape, self.materials.reference_temperature_K, dtype=np.float64)
        else:
            t = np.asarray(initial_temperature_K, dtype=np.float64).copy()
        tmax = self.materials.ti_Tmax_K if np.any(self.ag_fraction < 1.0 - 1e-14) else self.materials.Tl_K + 5000.0
        for _ in range(maximum_iterations):
            h = self.materials.enthalpy_density_J_m3(t, self.ag_fraction)
            c = self.materials.volumetric_capacity_J_m3_K(t, self.ag_fraction)
            delta = (target - h) / c
            scale = np.maximum(np.maximum(np.abs(target), np.abs(h)), 1.0)
            if float(np.max(np.abs(target - h) / scale)) <= tolerance:
                return t
            # A monotone bracket is maintained by a conservative Newton step.
            trial = t + np.clip(delta, -100.0, 100.0)
            trial = np.clip(trial, 100.0, tmax)
            if float(np.max(np.abs(trial - t))) <= 1.0e-13:
                break
            t = trial
        h = self.materials.enthalpy_density_J_m3(t, self.ag_fraction)
        scale = np.maximum(np.maximum(np.abs(target), np.abs(h)), 1.0)
        if float(np.max(np.abs(target - h) / scale)) > 1.0e-9:
            raise RuntimeError("enthalpy-to-temperature inversion did not converge")
        return t

    def _boundary_face_conductance_density(
        self, conductivity_W_m_K: np.ndarray, boundary: BoundarySet
    ) -> dict[str, np.ndarray]:
        faces: dict[str, np.ndarray] = {}
        for key in boundary.dirichlet_K:
            axis = AXIS_NAMES.index(key[0])
            gb = np.zeros(self.shape, dtype=np.float64)
            selector = [slice(None)] * self.ndim
            selector[axis] = 0 if key.endswith("low") else -1
            area_over_volume = 1.0 / self.spacing_m[axis]
            gb[tuple(selector)] = (
                conductivity_W_m_K[tuple(selector)] * area_over_volume / (0.5 * self.spacing_m[axis])
            )
            faces[key] = gb
        return faces

    def _face_conductances_density(
        self, conductivity_W_m_K: np.ndarray, boundary: BoundarySet
    ) -> tuple[list[np.ndarray], dict[str, np.ndarray]]:
        faces: list[np.ndarray] = []
        sc_fraction = 1.0 - self.ag_fraction
        for axis in range(self.ndim):
            k_right = np.roll(conductivity_W_m_K, -1, axis=axis)
            ag_right = np.roll(self.ag_fraction, -1, axis=axis)
            unlike_phase_fraction = self.ag_fraction * (1.0 - ag_right) + sc_fraction * ag_right
            dx = self.spacing_m[axis]
            total_face_resistance = (
                0.5 * dx / conductivity_W_m_K + 0.5 * dx / k_right
                + self.interface_resistance_m2_K_W * unlike_phase_fraction
            )
            area_over_volume = 1.0 / dx
            g = area_over_volume / total_face_resistance
            if axis not in boundary.periodic_axes:
                last = [slice(None)] * self.ndim
                last[axis] = -1
                g[tuple(last)] = 0.0
            faces.append(g)
        return faces, self._boundary_face_conductance_density(conductivity_W_m_K, boundary)

    def _net_conduction_density(
        self,
        temperature_K: np.ndarray,
        conductance_faces: list[np.ndarray],
        boundary_faces: dict[str, np.ndarray],
        boundary: BoundarySet,
    ) -> np.ndarray:
        power_density = np.zeros(self.shape, dtype=np.float64)
        for axis, g in enumerate(conductance_faces):
            edge_flux_density = g * (np.roll(temperature_K, -1, axis=axis) - temperature_K)
            power_density += edge_flux_density - np.roll(edge_flux_density, 1, axis=axis)
        for key, gb in boundary_faces.items():
            power_density += gb * (float(boundary.dirichlet_K[key]) - temperature_K)
        return power_density

    def _apply_laplacian(
        self,
        field: np.ndarray,
        conductance_faces: list[np.ndarray],
        boundary_faces: dict[str, np.ndarray],
    ) -> np.ndarray:
        result = np.zeros(self.shape, dtype=np.float64)
        for axis, g in enumerate(conductance_faces):
            edge = g * (np.roll(field, -1, axis=axis) - field)
            result += np.roll(edge, 1, axis=axis) - edge
        for gb in boundary_faces.values():
            result += gb * field
        return result

    def _jacobi_diagonal(self, conductance_faces: list[np.ndarray], boundary_faces: dict[str, np.ndarray]) -> np.ndarray:
        diagonal = np.zeros(self.shape, dtype=np.float64)
        for axis, g in enumerate(conductance_faces):
            diagonal += g + np.roll(g, 1, axis=axis)
        for gb in boundary_faces.values():
            diagonal += gb
        return diagonal

    def _residual(
        self,
        temperature_K: np.ndarray,
        enthalpy_old_J_m3: np.ndarray,
        dt_s: float,
        source_W_m3: np.ndarray,
        boundary: BoundarySet,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, np.ndarray]]:
        k = self.materials.conductivity_W_m_K(temperature_K, self.ag_fraction)
        g, gb = self._face_conductances_density(k, boundary)
        h = self.materials.enthalpy_density_J_m3(temperature_K, self.ag_fraction)
        conduction = self._net_conduction_density(temperature_K, g, gb, boundary)
        residual = h - enthalpy_old_J_m3 - dt_s * (conduction + source_W_m3)
        return residual, g, gb, {"enthalpy": h, "conductivity": k, "conduction": conduction}

    @staticmethod
    def _relative_residual(residual: np.ndarray, *scales: np.ndarray) -> float:
        scale = max([float(np.max(np.abs(item))) for item in scales] + [1.0e5])
        return float(np.max(np.abs(residual)) / scale)

    def _boundary_out_rate_W(self, temperature_K: np.ndarray, k: np.ndarray, boundary: BoundarySet) -> float:
        boundary_faces = self._boundary_face_conductance_density(k, boundary)
        total_rate = 0.0
        for key, gb in boundary_faces.items():
            total_rate += float(np.sum(gb * (temperature_K - float(boundary.dirichlet_K[key]))) * self.volume_m3)
        return total_rate

    def advance(
        self,
        dt_s: float,
        source_W_m3: np.ndarray | float = 0.0,
        boundary: BoundarySet | None = None,
    ) -> dict[str, float]:
        if dt_s <= 0.0 or not np.isfinite(dt_s):
            raise ValueError("dt_s must be finite and positive")
        boundary = boundary or BoundarySet()
        boundary.validate(self.ndim)
        source = np.broadcast_to(np.asarray(source_W_m3, dtype=np.float64), self.shape)
        if np.any(~np.isfinite(source)):
            raise ValueError("source must be finite")
        h_old = self.enthalpy_density_J_m3
        t = self.temperature_K.copy()
        total_cg_iterations = 0
        max_cg_relative_residual = 0.0
        nonlinear_residual = np.inf
        start = time.perf_counter()

        for nonlinear_iteration in range(1, self.max_newton + 1):
            residual, gfaces, gbfaces, fields = self._residual(t, h_old, dt_s, source, boundary)
            nonlinear_residual = self._relative_residual(
                residual, h_old, dt_s * source, dt_s * fields["conduction"]
            )
            if nonlinear_residual <= self.nonlinear_rtol:
                break
            capacity = self.materials.volumetric_capacity_J_m3_K(t, self.ag_fraction)
            def matvec(flat: np.ndarray) -> np.ndarray:
                delta_t = flat.reshape(self.shape)
                out = capacity * delta_t + dt_s * self._apply_laplacian(delta_t, gfaces, gbfaces)
                return out.ravel()

            operator = LinearOperator((self.ag_fraction.size, self.ag_fraction.size), matvec=matvec, dtype=np.float64)
            # FFT常系数周期拉普拉斯预条件器：以当前平均体积热容和面导热系数
            # 构造质量项+扩散项谱，适用于绝热/定温外边界下的异质修正方程。
            # 预条件器是正定对称的；只用于加速迭代，不改变离散方程。
            mean_capacity = float(np.mean(capacity))
            spectrum = np.full(tuple(self.shape[:-1]) + (self.shape[-1] // 2 + 1,), mean_capacity,
                               dtype=np.float64)
            for axis, g in enumerate(gfaces):
                n_axis = self.shape[axis]
                if axis == self.ndim - 1:
                    frequencies = np.arange(n_axis // 2 + 1, dtype=np.float64) / n_axis
                else:
                    frequencies = np.arange(n_axis, dtype=np.float64) / n_axis
                shape = [1] * self.ndim
                shape[axis] = frequencies.size
                lam = 2.0 * float(np.mean(g)) * (1.0 - np.cos(2.0 * np.pi * frequencies))
                spectrum += dt_s * lam.reshape(shape)
            preconditioner = LinearOperator(
                (self.ag_fraction.size, self.ag_fraction.size),
                matvec=lambda flat: fft.irfftn(
                    fft.rfftn(flat.reshape(self.shape)) / spectrum, s=self.shape
                ).ravel(), dtype=np.float64,
            )
            cg_iterations = [0]

            def cg_callback(x: np.ndarray) -> None:
                cg_iterations[0] += 1

            delta, info = cg(
                operator, -residual.ravel(), M=preconditioner,
                rtol=self.linear_rtol, atol=0.0, maxiter=self.max_cg,
                callback=cg_callback,
            )
            cg_iterations = cg_iterations[0]
            total_cg_iterations += cg_iterations
            linear_res = float(np.linalg.norm(matvec(delta) + residual.ravel()) / max(np.linalg.norm(residual), 1.0e-30))
            max_cg_relative_residual = max(max_cg_relative_residual, linear_res)
            if info != 0:
                raise RuntimeError(f"PCG did not converge (info={info}, iterations={cg_iterations}, rel={linear_res:.3e})")

            accepted = False
            damping = 1.0
            # 相变焓曲线在窄糊状区很陡，首个Newton步可能跨越数十K；
            # 允许充分回溯以守住单调残差下降，而不是强行裁剪液相率。
            for _ in range(32):
                candidate = t + damping * delta.reshape(self.shape)
                candidate_max = self.materials.ti_Tmax_K if np.any(self.ag_fraction < 1.0 - 1e-14) else 10000.0
                if float(candidate.min()) < 100.0 or float(candidate.max()) > candidate_max + 1e-8:
                    damping *= 0.5
                    continue
                new_residual, _, _, new_fields = self._residual(candidate, h_old, dt_s, source, boundary)
                candidate_norm = self._relative_residual(
                    new_residual, h_old, dt_s * source, dt_s * new_fields["conduction"]
                )
                if candidate_norm < nonlinear_residual or candidate_norm <= self.nonlinear_rtol:
                    t = candidate
                    nonlinear_residual = candidate_norm
                    accepted = True
                    break
                damping *= 0.5
            if not accepted:
                raise RuntimeError(
                    f"Newton line search failed at iteration {nonlinear_iteration}; "
                    f"relative residual={nonlinear_residual:.3e}"
                )
            if nonlinear_residual <= self.nonlinear_rtol:
                break
        else:
            raise RuntimeError(
                f"Newton did not converge in {self.max_newton} iterations; relative residual={nonlinear_residual:.3e}"
            )

        final_residual, _, _, final_fields = self._residual(t, h_old, dt_s, source, boundary)
        nonlinear_residual = self._relative_residual(
            final_residual, h_old, dt_s * source, dt_s * final_fields["conduction"]
        )
        if nonlinear_residual > max(self.nonlinear_rtol * 1.01, 1.0e-12):
            raise RuntimeError(f"post-step nonlinear residual is too high: {nonlinear_residual:.3e}")
        if float(np.min(self.materials.liquid_fraction(t))) < -1e-12 or float(np.max(self.materials.liquid_fraction(t))) > 1.0 + 1e-12:
            raise RuntimeError("nonphysical Ag liquid fraction; no clipping is applied")
        if np.any(self.ag_fraction < 1.0 - 1e-14) and float(t.max()) > self.materials.ti_Tmax_K + 1e-8:
            raise RuntimeError("Ti3SiC2 temperature exceeded the audited candidate range")

        source_rate_W = float(np.sum(source) * self.volume_m3)
        out_rate_W = self._boundary_out_rate_W(t, final_fields["conductivity"], boundary)
        source_increment = dt_s * source_rate_W
        out_increment = dt_s * out_rate_W
        storage_change = float(np.sum(final_fields["enthalpy"] - h_old) * self.volume_m3)
        self.temperature_K = t
        self.enthalpy_density_J_m3 = final_fields["enthalpy"]
        self.time_s += dt_s
        self.source_energy_J += source_increment
        self.boundary_energy_out_J += out_increment
        self.step_count += 1
        storage_total = float(np.sum(self.enthalpy_density_J_m3) * self.volume_m3)
        energy_residual = (storage_total - self.initial_stored_enthalpy_J
                           - self.source_energy_J + self.boundary_energy_out_J)
        return {
            "time_s": self.time_s,
            "dt_s": dt_s,
            "nonlinear_iterations": nonlinear_iteration,
            "pcg_iterations_total": total_cg_iterations,
            "pcg_relative_residual_max": max_cg_relative_residual,
            "nonlinear_relative_residual": nonlinear_residual,
            "max_temperature_K": float(t.max()),
            "mean_ag_temperature_K": float(np.sum(self.ag_fraction * t) / max(np.sum(self.ag_fraction), 1.0e-30)),
            "max_ag_liquid_fraction": float(self.materials.liquid_fraction(t).max()),
            "ag_liquid_volume_fraction": float(np.sum(self.ag_fraction * self.materials.liquid_fraction(t)) / max(np.sum(self.ag_fraction), 1.0e-30)),
            "source_energy_step_J": source_increment,
            "boundary_energy_out_step_J": out_increment,
            "storage_change_step_J": storage_change,
            "total_stored_enthalpy_J": storage_total,
            "cumulative_source_energy_J": self.source_energy_J,
            "cumulative_boundary_out_J": self.boundary_energy_out_J,
            "energy_residual_J": energy_residual,
            "step_wall_seconds": time.perf_counter() - start,
        }
