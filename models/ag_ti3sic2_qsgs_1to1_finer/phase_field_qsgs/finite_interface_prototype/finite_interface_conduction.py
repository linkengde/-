"""Conservative two-phase finite-volume heat conduction with interface resistance.

This module is an isolated conduction-kernel prototype. It does not modify or
replace ``qsgs_enthalpy_solver.py``. Each voxel has one phase-local temperature;
unlike-phase neighboring voxels exchange heat through an explicit interfacial
conductance, which permits a temperature jump across that interface. Material
properties are cell-centered fields supplied by the caller and are held fixed
for each solve/step. Latent heat and temperature-dependent property updates are
deliberately not implemented in this small validation kernel.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


BoundaryKey = tuple[int, str]


@dataclass(frozen=True)
class StepAccounting:
    """Energy ledger for one explicit conduction step."""

    energy_before_J: float
    energy_after_J: float
    source_energy_J: float
    boundary_energy_J: float
    residual_J: float
    max_temperature_K: float


class FiniteInterfaceFV:
    """Cell-centered finite-volume conduction on a binary/labelled phase mask.

    Parameters
    ----------
    phase_id:
        Three-dimensional integer/bool phase labels. Faces between different
        labels receive the requested thermal boundary conductance.
    spacing_m:
        Orthogonal cell pitches `(dx, dy, dz)` in metres.
    conductivity_W_m_K, volumetric_capacity_J_m3_K:
        Positive scalar or fields broadcastable to `phase_id.shape`.
    interface_conductance_W_m2_K:
        Nonnegative finite value, zero (adiabatic interface), or positive
        infinity (perfect-contact limit).
    dirichlet_boundaries:
        Optional fixed-temperature outer faces keyed by `(axis, "low"|"high")`.
        Unspecified external faces are insulated.
    periodic_axes:
        Axes with periodic wrap faces. Dirichlet and periodic boundaries cannot
        both be specified on one axis.
    """

    def __init__(
        self,
        phase_id: np.ndarray,
        spacing_m: tuple[float, float, float],
        conductivity_W_m_K: np.ndarray | float,
        volumetric_capacity_J_m3_K: np.ndarray | float,
        interface_conductance_W_m2_K: float,
        dirichlet_boundaries: Mapping[BoundaryKey, float] | None = None,
        periodic_axes: tuple[int, ...] = (),
    ) -> None:
        phases = np.asarray(phase_id)
        if phases.ndim != 3 or any(n <= 0 for n in phases.shape):
            raise ValueError("phase_id must be a non-empty three-dimensional array")
        if not np.issubdtype(phases.dtype, np.integer) and phases.dtype != np.bool_:
            raise ValueError("phase_id must contain integer or boolean phase labels")
        spacing = np.asarray(spacing_m, dtype=float)
        if spacing.shape != (3,) or np.any(~np.isfinite(spacing)) or np.any(spacing <= 0.0):
            raise ValueError("spacing_m must contain three finite positive pitches")
        if np.isnan(interface_conductance_W_m2_K) or interface_conductance_W_m2_K < 0.0:
            raise ValueError("interface conductance must be nonnegative or positive infinity")

        shape = phases.shape
        conductivity = self._broadcast_field(conductivity_W_m_K, shape, "conductivity")
        capacity = self._broadcast_field(volumetric_capacity_J_m3_K, shape, "capacity")
        if np.any(~np.isfinite(conductivity)) or np.any(conductivity <= 0.0):
            raise ValueError("all phase conductivities must be finite and positive")
        if np.any(~np.isfinite(capacity)) or np.any(capacity <= 0.0):
            raise ValueError("all phase volumetric capacities must be finite and positive")

        periodic = tuple(sorted(set(int(axis) for axis in periodic_axes)))
        if any(axis not in (0, 1, 2) for axis in periodic):
            raise ValueError("periodic_axes entries must be 0, 1, or 2")
        dirichlet = dict(dirichlet_boundaries or {})
        for key, temperature in dirichlet.items():
            if (
                not isinstance(key, tuple)
                or len(key) != 2
                or key[0] not in (0, 1, 2)
                or key[1] not in ("low", "high")
            ):
                raise ValueError("boundary keys must be (axis, 'low'|'high')")
            if not np.isfinite(temperature) or temperature <= 0.0:
                raise ValueError("Dirichlet temperatures must be finite and positive")
        if any(axis in periodic for axis, _ in dirichlet):
            raise ValueError("an axis cannot be both periodic and Dirichlet")

        self.phase_id = phases.copy()
        self.shape = shape
        self.spacing_m = spacing
        self.conductivity_W_m_K = conductivity
        self.volumetric_capacity_J_m3_K = capacity
        self.interface_conductance_W_m2_K = float(interface_conductance_W_m2_K)
        self.periodic_axes = periodic
        self.dirichlet_boundaries = dirichlet
        self.cell_volume_m3 = float(np.prod(spacing))
        self.n_cells = int(np.prod(shape))

        edge_i: list[np.ndarray] = []
        edge_j: list[np.ndarray] = []
        edge_g: list[np.ndarray] = []
        for axis in range(3):
            self._append_interior_faces(axis, edge_i, edge_j, edge_g)
            if axis in periodic and shape[axis] > 1:
                self._append_periodic_faces(axis, edge_i, edge_j, edge_g)
        self._edge_i = np.concatenate(edge_i) if edge_i else np.empty(0, dtype=np.int64)
        self._edge_j = np.concatenate(edge_j) if edge_j else np.empty(0, dtype=np.int64)
        self._edge_g_W_K = np.concatenate(edge_g) if edge_g else np.empty(0, dtype=float)
        self._interface_edge = (
            self.phase_id.ravel()[self._edge_i] != self.phase_id.ravel()[self._edge_j]
        )

        boundary_cells: list[np.ndarray] = []
        boundary_g: list[np.ndarray] = []
        boundary_t: list[np.ndarray] = []
        for (axis, side), temperature in sorted(dirichlet.items()):
            cells, conductance = self._boundary_face_conductance(axis, side)
            boundary_cells.append(cells)
            boundary_g.append(conductance)
            boundary_t.append(np.full(cells.size, temperature, dtype=float))
        self._boundary_cells = np.concatenate(boundary_cells) if boundary_cells else np.empty(0, dtype=np.int64)
        self._boundary_g_W_K = np.concatenate(boundary_g) if boundary_g else np.empty(0, dtype=float)
        self._boundary_t_K = np.concatenate(boundary_t) if boundary_t else np.empty(0, dtype=float)

        diagonal = np.zeros(self.n_cells, dtype=float)
        np.add.at(diagonal, self._edge_i, self._edge_g_W_K)
        np.add.at(diagonal, self._edge_j, self._edge_g_W_K)
        np.add.at(diagonal, self._boundary_cells, self._boundary_g_W_K)
        self._conductance_diagonal_W_K = diagonal

    @staticmethod
    def _broadcast_field(values: np.ndarray | float, shape: tuple[int, ...], name: str) -> np.ndarray:
        try:
            return np.broadcast_to(np.asarray(values, dtype=float), shape).astype(float, copy=True)
        except ValueError as exc:
            raise ValueError(f"{name} must be scalar or broadcastable to phase_id.shape") from exc

    def _face_conductance(self, axis: int, i: np.ndarray, j: np.ndarray) -> np.ndarray:
        flat_k = self.conductivity_W_m_K.ravel()
        flat_phase = self.phase_id.ravel()
        dx = float(self.spacing_m[axis])
        area = self.cell_volume_m3 / dx
        resistance_area = dx / (2.0 * flat_k[i]) + dx / (2.0 * flat_k[j])
        unlike = flat_phase[i] != flat_phase[j]
        if self.interface_conductance_W_m2_K == 0.0:
            result = area / resistance_area
            result[unlike] = 0.0
            return result
        if np.isfinite(self.interface_conductance_W_m2_K):
            resistance_area[unlike] += 1.0 / self.interface_conductance_W_m2_K
        return area / resistance_area

    def _append_interior_faces(
        self,
        axis: int,
        edge_i: list[np.ndarray],
        edge_j: list[np.ndarray],
        edge_g: list[np.ndarray],
    ) -> None:
        if self.shape[axis] < 2:
            return
        ids = np.arange(self.n_cells, dtype=np.int64).reshape(self.shape)
        left = [slice(None)] * 3
        right = [slice(None)] * 3
        left[axis] = slice(None, -1)
        right[axis] = slice(1, None)
        i = ids[tuple(left)].ravel()
        j = ids[tuple(right)].ravel()
        edge_i.append(i)
        edge_j.append(j)
        edge_g.append(self._face_conductance(axis, i, j))

    def _append_periodic_faces(
        self,
        axis: int,
        edge_i: list[np.ndarray],
        edge_j: list[np.ndarray],
        edge_g: list[np.ndarray],
    ) -> None:
        ids = np.arange(self.n_cells, dtype=np.int64).reshape(self.shape)
        low = [slice(None)] * 3
        high = [slice(None)] * 3
        low[axis] = 0
        high[axis] = -1
        i = ids[tuple(high)].ravel()
        j = ids[tuple(low)].ravel()
        edge_i.append(i)
        edge_j.append(j)
        edge_g.append(self._face_conductance(axis, i, j))

    def _boundary_face_conductance(self, axis: int, side: str) -> tuple[np.ndarray, np.ndarray]:
        ids = np.arange(self.n_cells, dtype=np.int64).reshape(self.shape)
        slab = [slice(None)] * 3
        slab[axis] = 0 if side == "low" else -1
        cells = ids[tuple(slab)].ravel()
        area = self.cell_volume_m3 / float(self.spacing_m[axis])
        k = self.conductivity_W_m_K.ravel()[cells]
        conductance = area / (float(self.spacing_m[axis]) / (2.0 * k))
        return cells, conductance

    @property
    def interface_conductances_W_K(self) -> np.ndarray:
        """Conductances of all unlike-phase faces, in W/K."""
        return self._edge_g_W_K[self._interface_edge].copy()

    def _validated_temperature(self, temperature_K: np.ndarray) -> np.ndarray:
        temperature = np.asarray(temperature_K, dtype=float)
        if temperature.shape != self.shape:
            raise ValueError("temperature field shape does not match phase_id")
        if np.any(~np.isfinite(temperature)) or np.any(temperature <= 0.0):
            raise ValueError("temperatures must be finite and positive Kelvin values")
        return temperature

    def _validated_source(self, source_W_m3: np.ndarray | float | None) -> np.ndarray:
        if source_W_m3 is None:
            return np.zeros(self.shape, dtype=float)
        source = self._broadcast_field(source_W_m3, self.shape, "source")
        if np.any(~np.isfinite(source)):
            raise ValueError("source must be finite")
        return source

    def boundary_power_W(self, temperature_K: np.ndarray) -> float:
        """Net heat entering through all fixed-temperature outer faces."""
        temperature = self._validated_temperature(temperature_K).ravel()
        return float(np.sum(self._boundary_g_W_K * (self._boundary_t_K - temperature[self._boundary_cells])))

    def net_power_W(self, temperature_K: np.ndarray, source_W_m3: np.ndarray | float | None = None) -> np.ndarray:
        """Return net cell power, including sources, internal links, and boundaries."""
        temperature = self._validated_temperature(temperature_K).ravel()
        source = self._validated_source(source_W_m3).ravel()
        net = source * self.cell_volume_m3
        link_power = self._edge_g_W_K * (temperature[self._edge_j] - temperature[self._edge_i])
        np.add.at(net, self._edge_i, link_power)
        np.add.at(net, self._edge_j, -link_power)
        boundary_power = self._boundary_g_W_K * (
            self._boundary_t_K - temperature[self._boundary_cells]
        )
        np.add.at(net, self._boundary_cells, boundary_power)
        return net.reshape(self.shape)

    def energy_J(self, temperature_K: np.ndarray) -> float:
        """Sensible energy for the fixed-capacity prototype (zero at 0 K)."""
        temperature = self._validated_temperature(temperature_K)
        return float(np.sum(self.volumetric_capacity_J_m3_K * temperature) * self.cell_volume_m3)

    def stable_timestep_s(self, safety: float = 0.45) -> float:
        """Sufficient forward-Euler bound for the current fixed conductances."""
        if not np.isfinite(safety) or not (0.0 < safety <= 1.0):
            raise ValueError("safety must lie in (0, 1]")
        active = self._conductance_diagonal_W_K > 0.0
        if not np.any(active):
            return float("inf")
        capacity_J_K = self.volumetric_capacity_J_m3_K.ravel() * self.cell_volume_m3
        return float(safety * np.min(capacity_J_K[active] / self._conductance_diagonal_W_K[active]))

    def step(
        self,
        temperature_K: np.ndarray,
        dt_s: float,
        source_W_m3: np.ndarray | float | None = None,
        stability_safety: float = 0.45,
    ) -> tuple[np.ndarray, StepAccounting]:
        """Advance one explicit step and return a closed energy ledger."""
        if not np.isfinite(dt_s) or dt_s <= 0.0:
            raise ValueError("dt_s must be finite and positive")
        temperature = self._validated_temperature(temperature_K)
        source = self._validated_source(source_W_m3)
        stable_dt = self.stable_timestep_s(stability_safety)
        if dt_s > stable_dt * (1.0 + 1.0e-12):
            raise ValueError(f"dt_s exceeds the explicit stability bound {stable_dt:.9g} s")

        before = self.energy_J(temperature)
        source_power_W = float(np.sum(source) * self.cell_volume_m3)
        boundary_power_W = self.boundary_power_W(temperature)
        net_power = self.net_power_W(temperature, source)
        capacity_J_K = self.volumetric_capacity_J_m3_K * self.cell_volume_m3
        advanced = temperature + dt_s * net_power / capacity_J_K
        if np.any(~np.isfinite(advanced)) or np.any(advanced <= 0.0):
            raise FloatingPointError("finite-volume step produced invalid temperature")
        after = self.energy_J(advanced)
        source_energy = source_power_W * dt_s
        boundary_energy = boundary_power_W * dt_s
        residual = after - before - source_energy - boundary_energy
        accounting = StepAccounting(
            energy_before_J=before,
            energy_after_J=after,
            source_energy_J=source_energy,
            boundary_energy_J=boundary_energy,
            residual_J=residual,
            max_temperature_K=float(advanced.max()),
        )
        return advanced, accounting

    def steady_state(self, source_W_m3: np.ndarray | float | None = None) -> np.ndarray:
        """Solve the linear steady conduction system for fixed properties.

        At least one Dirichlet face is required to remove the constant null
        mode. A nonzero net source in a fully insulated system has no steady
        state and is rejected.
        """
        if self._boundary_cells.size == 0:
            raise ValueError("steady_state requires at least one Dirichlet boundary")
        source = self._validated_source(source_W_m3).ravel()
        boundary_g = self._boundary_g_W_K
        edge_g = self._edge_g_W_K
        edge_i, edge_j = self._edge_i, self._edge_j
        rows = np.concatenate((edge_i, edge_j, edge_i, edge_j, self._boundary_cells))
        cols = np.concatenate((edge_i, edge_j, edge_j, edge_i, self._boundary_cells))
        data = np.concatenate((edge_g, edge_g, -edge_g, -edge_g, boundary_g))
        matrix = coo_matrix((data, (rows, cols)), shape=(self.n_cells, self.n_cells)).tocsr()
        rhs = source * self.cell_volume_m3
        np.add.at(rhs, self._boundary_cells, boundary_g * self._boundary_t_K)
        solution = np.asarray(spsolve(matrix, rhs), dtype=float)
        if solution.shape != (self.n_cells,) or np.any(~np.isfinite(solution)):
            raise RuntimeError("steady conduction solve is singular or non-finite")
        return solution.reshape(self.shape)
