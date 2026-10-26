"""Connectivity-preserving Ag-only continuum TTM backend for native QSGS voxels."""

from __future__ import annotations

from pathlib import Path
import time

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import LinearOperator, cg

from ag_ttm_reference import AgParameters
from qsgs_component_coarsegrain import ComponentGraph, build_component_graph
from qsgs_enthalpy_solver import EV_PER_PS_ANGSTROM2_TO_W_M2
from qsgs_source_integration import integrate_qsgs_source_power_W
from thermophysical_model import AgEnthalpy


EV_J = 1.602176634e-19
ANGSTROM3_TO_M3 = 1.0e-30


class ReusableGraphOperator:
    """Fixed CSR sparsity with per-step edge-weight updates."""

    def __init__(self, graph: ComponentGraph) -> None:
        self.node_count = int(graph.voxel_count.size)
        self.edge_i = np.asarray(graph.edge_i, dtype=np.int32)
        self.edge_j = np.asarray(graph.edge_j, dtype=np.int32)
        self.edge_count = int(self.edge_i.size)
        n = self.node_count
        e = self.edge_count
        rows = np.concatenate((np.arange(n, dtype=np.int32), self.edge_i, self.edge_j))
        cols = np.concatenate((np.arange(n, dtype=np.int32), self.edge_j, self.edge_i))
        # CSR indices need sorted columns. Keep record IDs only long enough to
        # build direct positions for each diagonal and directed edge.
        order = np.lexsort((cols, rows))
        row_counts = np.bincount(rows, minlength=n)
        indptr = np.empty(n + 1, dtype=np.int64)
        indptr[0] = 0
        np.cumsum(row_counts, out=indptr[1:])
        if indptr[-1] > np.iinfo(np.int32).max:
            raise MemoryError("Graph CSR exceeds the SciPy int32 index range")
        indices = cols[order]
        record_to_position = np.empty(n + 2 * e, dtype=np.int64)
        record_to_position[order] = np.arange(n + 2 * e, dtype=np.int64)
        self.diagonal_positions = record_to_position[:n].astype(np.int32)
        self.edge_ij_positions = record_to_position[n : n + e].astype(np.int32)
        self.edge_ji_positions = record_to_position[n + e :].astype(np.int32)
        data = np.zeros(n + 2 * e, dtype=float)
        self.matrix = csr_matrix(
            (data, indices.astype(np.int32, copy=False), indptr.astype(np.int32)),
            shape=(n, n),
            copy=False,
        )
        del rows, cols, order, row_counts, indptr, indices, record_to_position
        self.conductance_sum = np.zeros(n, dtype=float)

    def update_system(self, edge_conductance: np.ndarray, mass_diagonal: np.ndarray, dt_ps: float) -> None:
        conductance = np.asarray(edge_conductance, dtype=float)
        if conductance.shape != (self.edge_count,) or np.any(conductance < 0.0):
            raise ValueError("Edge conductance has invalid shape or values")
        self.conductance_sum.fill(0.0)
        self.conductance_sum += np.bincount(
            self.edge_i, weights=conductance, minlength=self.node_count
        )
        self.conductance_sum += np.bincount(
            self.edge_j, weights=conductance, minlength=self.node_count
        )
        data = self.matrix.data
        data[self.diagonal_positions] = mass_diagonal + dt_ps * self.conductance_sum
        data[self.edge_ij_positions] = -dt_ps * conductance
        data[self.edge_ji_positions] = -dt_ps * conductance

    def inverse_diagonal(self) -> np.ndarray:
        return 1.0 / np.maximum(self.matrix.diagonal(), 1.0e-30)


class QSGSComponentTTMBackend:
    """Dual-temperature Ag system with component nodes and exact voxel source.

    This first production backend supports boundary case A (the scaffold is an
    adiabatic geometric obstacle and does not store heat), x/y periodicity,
    and insulated z faces. It keeps electron/lattice states on component nodes
    and expands fields to the native QSGS voxel grid only when snapshots are
    requested.
    """

    def __init__(
        self,
        source_geometry_npz: str | Path,
        electron_table: str | Path,
        boundary_case: str = "A",
        peak_flux_eV_ps_A2: float = 5.0,
        gaussian_radius_A: float = 30.0,
        absorption_depth_A: float = 20.0,
        tile_xy: int = 1,
        coarse_voxels: tuple[int, int, int] = (2, 1, 1),
        linear_preconditioner: str = "jacobi",
        initial_electron_temperature_K: float = 300.0,
        initial_lattice_temperature_K: float = 300.0,
    ) -> None:
        if boundary_case != "A":
            raise ValueError("The component-graph backend currently supports boundary case A only")
        if not isinstance(tile_xy, (int, np.integer)) or tile_xy <= 0:
            raise ValueError("tile_xy must be a positive integer")
        self.boundary_case = boundary_case
        self.source_geometry_npz = Path(source_geometry_npz)
        self.electron_table = Path(electron_table)
        self.tile_xy = int(tile_xy)
        self.coarse_voxels = tuple(int(v) for v in coarse_voxels)
        if linear_preconditioner not in {"jacobi", "amg"}:
            raise ValueError("linear_preconditioner must be 'jacobi' or 'amg'")
        self.linear_preconditioner = linear_preconditioner
        self.peak_flux_eV_ps_A2 = float(peak_flux_eV_ps_A2)
        self.peak_flux_W_m2 = self.peak_flux_eV_ps_A2 * EV_PER_PS_ANGSTROM2_TO_W_M2
        self.gaussian_radius_A = float(gaussian_radius_A)
        self.absorption_depth_A = float(absorption_depth_A)
        with np.load(self.source_geometry_npz, allow_pickle=False) as data:
            base_ag = ~np.asarray(data["skeleton"], dtype=bool)
            base_size_A = np.asarray(data["size_A"], dtype=float)
        self.base_shape_xyz = tuple(int(v) for v in base_ag.shape)
        self.source_ag = np.tile(base_ag, (self.tile_xy, self.tile_xy, 1))
        self.size_A = base_size_A * np.asarray([self.tile_xy, self.tile_xy, 1.0])
        self.spacing_A = self.size_A / np.asarray(self.source_ag.shape, dtype=float)
        self.voxel_volume_A3 = float(np.prod(self.spacing_A))
        self.graph_build_started_s = time.perf_counter()
        self.graph = build_component_graph(
            self.source_ag,
            self.size_A,
            self.coarse_voxels,
            intrinsic_k=1.0,
            periodic_axes=(0, 1),
        )
        self.graph_build_seconds = time.perf_counter() - self.graph_build_started_s
        self.node_volume_A3 = self.graph.voxel_count.astype(float) * self.voxel_volume_A3
        self.parameters = AgParameters.from_csv(
            self.electron_table, allow_low_temperature_extension=False
        )
        self.lattice = AgEnthalpy()
        if not (
            self.parameters.Tmin <= initial_electron_temperature_K <= self.parameters.Tmax
            and 300.0 <= initial_lattice_temperature_K <= self.lattice.melting_temperature_K
        ):
            raise ValueError("Initial temperature is outside the supported Ag property range")
        source_voxel_W = integrate_qsgs_source_power_W(
            self.source_ag,
            self.source_ag.shape,
            self.size_A,
            self.peak_flux_W_m2,
            self.gaussian_radius_A,
            self.absorption_depth_A,
        )
        self.source_voxel_power_W = source_voxel_W
        ids = self.graph.component_id[self.source_ag]
        self.source_node_W = np.bincount(
            ids,
            weights=source_voxel_W[self.source_ag],
            minlength=self.graph.voxel_count.size,
        )
        self.source_power_W = float(source_voxel_W.sum())
        self.electron_temperature_K = np.full(
            self.graph.voxel_count.size, float(initial_electron_temperature_K)
        )
        self.lattice_temperature_K = np.full(
            self.graph.voxel_count.size, float(initial_lattice_temperature_K)
        )
        self.lattice_liquid_fraction = np.zeros(self.graph.voxel_count.size, dtype=float)
        self.cumulative_coupling_eV = np.zeros(self.graph.voxel_count.size, dtype=float)
        self.operator = ReusableGraphOperator(self.graph)
        self._amg_solver = None
        self._amg_preconditioner = None
        self.amg_setup_seconds = 0.0
        self.time_ps = 0.0
        self.injected_energy_J = 0.0
        self.coupling_energy_J = 0.0
        self.electron_diffusion_energy_J = 0.0
        self.lattice_diffusion_energy_J = 0.0
        self.initial_total_energy_J = self.total_energy_J
        self.field_snapshots: dict[float, dict[str, np.ndarray]] = {}
        self.last_diffusion_solve_seconds = 0.0
        self.last_coupling_solve_seconds = 0.0
        self.last_diffusion_cg_iterations = 0
        self.total_diffusion_solve_seconds = 0.0
        self.total_coupling_solve_seconds = 0.0
        self.total_diffusion_cg_iterations = 0

    @property
    def electron_energy_eV(self) -> np.ndarray:
        return self.node_volume_A3 * self.parameters.ue(self.electron_temperature_K)

    @property
    def lattice_enthalpy_eV(self) -> np.ndarray:
        h_mass = self.lattice.specific_enthalpy_J_kg(
            self.lattice_temperature_K, self.lattice_liquid_fraction
        )
        h_density = (
            self.lattice.reference_density_kg_m3
            * np.asarray(h_mass)
            * ANGSTROM3_TO_M3
            / EV_J
        )
        return self.node_volume_A3 * h_density

    @property
    def electron_energy_J(self) -> float:
        return float(self.electron_energy_eV.sum() * EV_J)

    @property
    def lattice_energy_J(self) -> float:
        return float(self.lattice_enthalpy_eV.sum() * EV_J)

    @property
    def total_energy_J(self) -> float:
        return self.electron_energy_J + self.lattice_energy_J

    @property
    def ag_volume_fraction(self) -> float:
        return float(self.graph.voxel_count.sum() / self.source_ag.size)

    def _h_vol_eV_A3(self, temperature_K: np.ndarray, liquid_fraction: np.ndarray) -> np.ndarray:
        h_mass = self.lattice.specific_enthalpy_J_kg(temperature_K, liquid_fraction)
        return (
            self.lattice.reference_density_kg_m3
            * np.asarray(h_mass)
            * ANGSTROM3_TO_M3
            / EV_J
        )

    def _state_from_h_vol(self, h_vol_eV_A3: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        h_mass = (
            np.asarray(h_vol_eV_A3)
            * EV_J
            / (self.lattice.reference_density_kg_m3 * ANGSTROM3_TO_M3)
        )
        state = self.lattice.state_from_specific_enthalpy(h_mass)
        return np.asarray(state.temperature_K), np.asarray(state.liquid_fraction)

    def _source_half_step(self, dt_ps: float) -> float:
        energy_eV = self.source_node_W * (0.5 * dt_ps * 1.0e-12 / EV_J)
        current = self.electron_energy_eV + energy_eV
        self.electron_temperature_K = self.parameters.temperature_from_ue(
            current / self.node_volume_A3
        )
        return float(energy_eV.sum() * EV_J)

    def _electron_diffusion(self, dt_ps: float) -> float:
        if dt_ps <= 0.0 or self.graph.edge_i.size == 0:
            return 0.0
        started = time.perf_counter()
        t0 = self.electron_temperature_K.copy()
        u0 = self.parameters.ue(t0)
        energy0 = self.node_volume_A3 * u0
        edge_i = self.operator.edge_i
        edge_j = self.operator.edge_j
        ki = self.parameters.ke(t0[edge_i])
        kj = self.parameters.ke(t0[edge_j])
        base = self.graph.edge_conductance
        edge_k = 2.0 * ki * kj / (ki + kj)
        physical_g = base * edge_k
        temperature = t0.copy()
        cg_iterations = 0
        for _ in range(16):
            residual = self.node_volume_A3 * (self.parameters.ue(temperature) - u0)
            residual += dt_ps * (self._apply_laplacian(physical_g, temperature))
            residual_norm = float(np.max(np.abs(residual)))
            if residual_norm <= max(2.0e-13, float(np.max(np.abs(energy0))) * 2.0e-10):
                break
            slopes = self._ue_slope(temperature, self.parameters)
            mass = self.node_volume_A3 * slopes
            self.operator.update_system(physical_g, mass, dt_ps)
            if self.linear_preconditioner == "amg":
                if self._amg_preconditioner is None:
                    try:
                        import pyamg
                    except ImportError as exc:
                        raise RuntimeError(
                            "AMG preconditioning requires the optional pyamg package"
                        ) from exc
                    setup_started = time.perf_counter()
                    precondition_matrix = self.operator.matrix.copy()
                    precondition_matrix.sum_duplicates()
                    precondition_matrix.sort_indices()
                    self._amg_solver = pyamg.smoothed_aggregation_solver(
                        precondition_matrix, symmetry="symmetric"
                    )
                    self._amg_preconditioner = self._amg_solver.aspreconditioner(cycle="V")
                    self.amg_setup_seconds += time.perf_counter() - setup_started
                preconditioner = self._amg_preconditioner
            else:
                inverse_diagonal = self.operator.inverse_diagonal()
                preconditioner = LinearOperator(
                    self.operator.matrix.shape, matvec=lambda vector, inv=inverse_diagonal: inv * vector
                )
            def count_iteration(_: np.ndarray) -> None:
                nonlocal cg_iterations
                cg_iterations += 1

            delta, info = cg(
                self.operator.matrix,
                -residual,
                M=preconditioner,
                rtol=1.0e-8,
                atol=max(1.0e-14, float(np.max(np.abs(energy0))) * 1.0e-12),
                maxiter=min(3000, max(300, 4 * len(temperature))),
                callback=count_iteration,
            )
            if info != 0:
                raise RuntimeError(f"Graph electron diffusion CG failed (info={info})")
            old_norm = float(np.linalg.norm(residual))
            alpha = 1.0
            accepted = False
            while alpha >= 1.0 / 128.0:
                trial = temperature + alpha * delta
                # Round-off in the sparse solve can put unchanged 300 K nodes a
                # few ulps below the tabulated endpoint. Accept and clamp only
                # that negligible undershoot; real range excursions still
                # trigger line-search backtracking.
                bound_tolerance_K = 1.0e-9
                if np.any(trial < self.parameters.Tmin - bound_tolerance_K) or np.any(
                    trial > self.parameters.Tmax + bound_tolerance_K
                ):
                    alpha *= 0.5
                    continue
                trial = np.clip(trial, self.parameters.Tmin, self.parameters.Tmax)
                trial_residual = self.node_volume_A3 * (self.parameters.ue(trial) - u0)
                trial_residual += dt_ps * self._apply_laplacian(physical_g, trial)
                if float(np.linalg.norm(trial_residual)) < old_norm:
                    temperature = trial
                    accepted = True
                    break
                alpha *= 0.5
            if not accepted:
                raise RuntimeError("Graph electron diffusion Newton line search failed")
        else:
            raise RuntimeError("Graph electron diffusion Newton solve did not converge")
        self.electron_temperature_K = temperature
        self.last_diffusion_solve_seconds = time.perf_counter() - started
        self.last_diffusion_cg_iterations = cg_iterations
        self.total_diffusion_solve_seconds += self.last_diffusion_solve_seconds
        self.total_diffusion_cg_iterations += cg_iterations
        return float((self.electron_energy_eV.sum() - energy0.sum()) * EV_J)

    @staticmethod
    def _ue_slope(temperature_K: np.ndarray, parameters: AgParameters) -> np.ndarray:
        t = parameters.temperature_K
        slopes = np.diff(parameters.ue_eV_A3) / np.diff(t)
        index = np.clip(np.searchsorted(t, temperature_K, side="right") - 1, 0, len(slopes) - 1)
        return slopes[index]

    def _apply_laplacian(self, edge_conductance: np.ndarray, temperature_K: np.ndarray) -> np.ndarray:
        i, j = self.operator.edge_i, self.operator.edge_j
        flux = np.asarray(edge_conductance) * (temperature_K[j] - temperature_K[i])
        result = np.bincount(i, weights=-flux, minlength=temperature_K.size)
        result += np.bincount(j, weights=flux, minlength=temperature_K.size)
        return result

    def _implicit_local_coupling(self, dt_ps: float) -> float:
        started = time.perf_counter()
        te0 = self.electron_temperature_K
        tl0 = self.lattice_temperature_K
        fl0 = self.lattice_liquid_fraction
        h0 = self._h_vol_eV_A3(tl0, fl0)
        total_density = self.parameters.ue(te0) + h0
        low = np.minimum(te0, tl0).copy()
        high = np.maximum(te0, tl0).copy()

        def residual_at(te: np.ndarray):
            h_trial = total_density - self.parameters.ue(te)
            tl, fl = self._state_from_h_vol(h_trial)
            residual = self.parameters.ue(te) - self.parameters.ue(te0)
            residual += dt_ps * self.parameters.g(te) * (te - tl)
            return residual, tl, fl

        f_low = residual_at(low)[0]
        f_high = residual_at(high)[0]
        if np.any(f_low > 2.0e-12) or np.any(f_high < -2.0e-12):
            raise RuntimeError("Graph electron-lattice coupling root is not bracketed")
        for _ in range(40):
            middle = 0.5 * (low + high)
            f_middle = residual_at(middle)[0]
            positive = f_middle > 0.0
            high = np.where(positive, middle, high)
            low = np.where(positive, low, middle)
        te1 = 0.5 * (low + high)
        h1 = total_density - self.parameters.ue(te1)
        tl1, fl1 = self._state_from_h_vol(h1)
        transfer_eV = self.node_volume_A3 * (self.parameters.ue(te0) - self.parameters.ue(te1))
        self.electron_temperature_K = te1
        self.lattice_temperature_K = tl1
        self.lattice_liquid_fraction = fl1
        self.cumulative_coupling_eV += transfer_eV
        self.last_coupling_solve_seconds = time.perf_counter() - started
        self.total_coupling_solve_seconds += self.last_coupling_solve_seconds
        return float(transfer_eV.sum() * EV_J)

    def _expand_node_field(self, values: np.ndarray, fill_value: float = 0.0) -> np.ndarray:
        field = np.full(self.source_ag.shape, fill_value, dtype=float)
        active = self.source_ag
        field[active] = values[self.graph.component_id[active]]
        return field

    def capture_snapshot(self, time_ps: float) -> None:
        self.field_snapshots[float(time_ps)] = {
            "electron_temperature_K": self._expand_node_field(self.electron_temperature_K, 300.0),
            "lattice_temperature_K": self._expand_node_field(self.lattice_temperature_K, 300.0),
            "ag_liquid_fraction": self._expand_node_field(self.lattice_liquid_fraction, 0.0),
            "electron_energy_density_eV_A3": self._expand_node_field(
                self.parameters.ue(self.electron_temperature_K), 0.0
            ),
            "lattice_enthalpy_density_J_m3": self._expand_node_field(
                self._h_vol_eV_A3(self.lattice_temperature_K, self.lattice_liquid_fraction)
                * EV_J
                / ANGSTROM3_TO_M3,
                0.0,
            ),
            "electron_to_lattice_cumulative_energy_density_J_m3": self._expand_node_field(
                self.cumulative_coupling_eV / self.node_volume_A3
                * EV_J
                / ANGSTROM3_TO_M3,
                0.0,
            ),
        }

    def advance(
        self,
        dt_ps: float,
        heat_on: bool = True,
        coupling_on: bool = True,
        electron_diffusion_on: bool = True,
    ) -> dict[str, float]:
        if not np.isfinite(dt_ps) or dt_ps <= 0.0:
            raise ValueError("dt_ps must be finite and positive")
        initial = self.total_energy_J
        source_J = 0.0
        diffusion_J = 0.0
        coupling_J = 0.0
        if heat_on:
            source_J += self._source_half_step(dt_ps)
        if electron_diffusion_on:
            diffusion_J += self._electron_diffusion(0.5 * dt_ps)
        if coupling_on:
            coupling_J = self._implicit_local_coupling(dt_ps)
        if electron_diffusion_on:
            diffusion_J += self._electron_diffusion(0.5 * dt_ps)
        if heat_on:
            source_J += self._source_half_step(dt_ps)
        self.time_ps += dt_ps
        self.injected_energy_J += source_J
        self.coupling_energy_J += coupling_J
        self.electron_diffusion_energy_J += diffusion_J
        final = self.total_energy_J
        residual = final - initial - source_J - diffusion_J
        return {
            "time_ps": self.time_ps,
            "dt_ps": dt_ps,
            "electron_substeps": 1.0,
            "max_electron_temperature_K": float(self.electron_temperature_K.max()),
            "mean_electron_temperature_K": float(
                np.average(self.electron_temperature_K, weights=self.node_volume_A3)
            ),
            "max_lattice_temperature_K": float(self.lattice_temperature_K.max()),
            "mean_ag_lattice_temperature_K": float(
                np.average(self.lattice_temperature_K, weights=self.node_volume_A3)
            ),
            "max_ag_liquid_fraction": float(self.lattice_liquid_fraction.max()),
            "ag_liquid_volume_fraction": float(
                np.dot(self.node_volume_A3, self.lattice_liquid_fraction)
                / self.node_volume_A3.sum()
            ),
            "electron_energy_J": self.electron_energy_J,
            "lattice_enthalpy_J": self.lattice_energy_J,
            "total_energy_J": final,
            "source_energy_J": source_J,
            "electron_to_lattice_energy_J": coupling_J,
            "electron_diffusion_net_J": diffusion_J,
            "lattice_diffusion_net_J": 0.0,
            "energy_balance_residual_J": residual,
        }

    def run(
        self,
        duration_ps: float,
        dt_ps: float,
        heat_on: bool,
        sample_interval_ps: float | None = None,
        snapshot_times_ps: tuple[float, ...] | list[float] = (),
        coupling_on: bool = True,
        electron_diffusion_on: bool = True,
    ) -> list[dict[str, float]]:
        if duration_ps <= 0.0 or dt_ps <= 0.0:
            raise ValueError("duration_ps and dt_ps must be positive")
        sample_interval_ps = dt_ps if sample_interval_ps is None else float(sample_interval_ps)
        if sample_interval_ps <= 0.0:
            raise ValueError("sample_interval_ps must be positive")
        end_ps = self.time_ps + duration_ps
        tolerance = max(1.0e-12, dt_ps * 1.0e-9)
        snapshots = sorted(float(value) for value in snapshot_times_ps)
        if any(value < self.time_ps - tolerance or value > end_ps + tolerance for value in snapshots):
            raise ValueError("Snapshot times must lie within this run interval")
        while snapshots and abs(snapshots[0] - self.time_ps) <= tolerance:
            self.capture_snapshot(snapshots.pop(0))
        next_sample = self.time_ps + sample_interval_ps
        sample_start = self.time_ps
        accumulated = {
            key: 0.0
            for key in (
                "electron_substeps",
                "source_energy_J",
                "electron_to_lattice_energy_J",
                "electron_diffusion_net_J",
                "lattice_diffusion_net_J",
                "energy_balance_residual_J",
            )
        }
        rows: list[dict[str, float]] = []
        while self.time_ps < end_ps - max(1.0e-15, duration_ps * 1.0e-13):
            limits = [dt_ps, end_ps - self.time_ps, next_sample - self.time_ps]
            if snapshots:
                limits.append(snapshots[0] - self.time_ps)
            step = min(limits)
            if step <= 0.0:
                raise RuntimeError("Sampling schedule produced a non-positive step")
            row = self.advance(step, heat_on, coupling_on, electron_diffusion_on)
            while snapshots and abs(self.time_ps - snapshots[0]) <= tolerance:
                self.capture_snapshot(snapshots.pop(0))
            for name in accumulated:
                accumulated[name] += row[name]
            if self.time_ps >= next_sample - tolerance or self.time_ps >= end_ps - tolerance:
                sampled = dict(row)
                sampled.update(accumulated)
                sampled["dt_ps"] = self.time_ps - sample_start
                rows.append(sampled)
                accumulated = {key: 0.0 for key in accumulated}
                sample_start = self.time_ps
                while next_sample <= self.time_ps + tolerance:
                    next_sample += sample_interval_ps
        return rows
