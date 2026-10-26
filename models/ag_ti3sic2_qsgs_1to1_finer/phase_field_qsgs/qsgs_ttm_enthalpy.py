"""Coupled Ag electron/lattice enthalpy solver on the projected QSGS grid.

This is a continuum two-temperature reference solver, not TTM-MD. Electron
parameters and their interpolation are reused from ``ag_ttm_reference.py``;
the Ag lattice is represented by the QSGS finite-volume enthalpy model.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import coo_matrix, diags, eye
from scipy.sparse.linalg import cg, gmres

MODEL_DIR = Path(__file__).resolve().parents[1]
if str(MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(MODEL_DIR))

from ag_ttm_reference import AgParameters, AgTTMGrid
from qsgs_enthalpy_solver import DATA_DIR, PHASE_FRACTION_ROUNDOFF_TOL, QSGSEnthalpyGrid


EV_J = 1.602176634e-19
DEFAULT_ELECTRON_TABLE = DATA_DIR.parent / "ttm_data" / "Ag_TTM_hybrid_figure_fit_300_20000K.csv"


class QSGSDoubleTemperatureGrid:
    """Ag-only electron TTM coupled conservatively to a mixture lattice H field.

    Electron source, Ce, ke, and G are all weighted by the Ag cell fraction in
    ``AgTTMGrid``. Lattice enthalpy includes only Ag in boundary A and Ag plus
    the heat-storing/conducting Ti3SiC2 phase in boundary B.
    """

    def __init__(
        self,
        projected_npz: Path,
        boundary_case: str,
        electron_table: Path = DEFAULT_ELECTRON_TABLE,
        peak_flux_eV_ps_A2: float = 5.0,
        gaussian_radius_A: float = 30.0,
        absorption_depth_A: float = 20.0,
        initial_lattice_temperature_K: float = 300.0,
        initial_electron_temperature_K: float | np.ndarray = 300.0,
        ag_lattice_conductivity_W_m_K: float = 0.0,
        bottom_sink_temperature_K: float | None = None,
        face_contacts_npz: Path | None = None,
        tile_xy: int = 1,
        source_geometry_npz: Path | None = None,
    ) -> None:
        self.lattice = QSGSEnthalpyGrid(
            projected_npz,
            boundary_case,
            peak_flux_eV_ps_A2,
            gaussian_radius_A,
            absorption_depth_A,
            tile_xy=tile_xy,
            source_geometry_npz=source_geometry_npz,
        )
        if not np.isfinite(ag_lattice_conductivity_W_m_K) or ag_lattice_conductivity_W_m_K < 0.0:
            raise ValueError("Ag lattice conductivity must be finite and non-negative")
        # The single-temperature enthalpy model retains total Ag conductivity.
        # In this two-temperature model, total conductivity must not be reused
        # in the lattice equation: Ag electron conduction is already explicit.
        self.ag_lattice_conductivity_W_m_K = float(ag_lattice_conductivity_W_m_K)
        if bottom_sink_temperature_K is not None and (
            not np.isfinite(bottom_sink_temperature_K) or bottom_sink_temperature_K <= 0.0
        ):
            raise ValueError("Bottom sink temperature must be finite and positive")
        self.bottom_sink_temperature_K = (
            None if bottom_sink_temperature_K is None else float(bottom_sink_temperature_K)
        )
        self.face_contacts_npz = None if face_contacts_npz is None else Path(face_contacts_npz)
        self.electron_face_contact_fractions: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None
        if self.face_contacts_npz is not None:
            if self.bottom_sink_temperature_K is not None:
                raise ValueError(
                    "Exact face contacts currently support insulated z boundaries only; "
                    "a bottom-sink surface contact fraction is not provided"
                )
            with np.load(projected_npz, allow_pickle=False) as projected, np.load(
                self.face_contacts_npz, allow_pickle=False
            ) as contacts:
                base_ag = np.asarray(projected["ag_fraction"], dtype=float)
                contact_ag = np.asarray(contacts["ag_fraction"], dtype=float)
                contact_spacing = np.asarray(contacts["spacing_A"], dtype=float)
                contact_size = np.asarray(contacts["size_A"], dtype=float)
                if contact_ag.shape != base_ag.shape or not np.allclose(
                    contact_ag, base_ag, rtol=0.0, atol=2.0e-12
                ):
                    raise ValueError("Face-contact file does not match the projected Ag geometry")
                if not np.allclose(contact_spacing, self.lattice.spacing_A, rtol=0.0, atol=1.0e-10):
                    raise ValueError("Face-contact grid spacing does not match the projected geometry")
                if not np.allclose(
                    contact_size,
                    np.asarray(projected["size_A"], dtype=float),
                    rtol=0.0,
                    atol=1.0e-8,
                ):
                    raise ValueError("Face-contact domain size does not match the projected geometry")
                base_contacts = tuple(
                    np.asarray(contacts[name], dtype=float)
                    for name in ("contact_x", "contact_y", "contact_z")
                )
            if any(values.shape != base_ag.shape for values in base_contacts):
                raise ValueError("Each face-contact field must match the base projected grid")
            if any(
                np.any(~np.isfinite(values))
                or np.any((values < -1.0e-12) | (values > 1.0 + 1.0e-12))
                for values in base_contacts
            ):
                raise ValueError("Face-contact fractions must be finite and lie in [0,1]")
            base_contacts = tuple(np.clip(values, 0.0, 1.0) for values in base_contacts)
            self.electron_face_contact_fractions = tuple(
                np.tile(values, (self.lattice.tile_xy, self.lattice.tile_xy, 1))
                for values in base_contacts
            )
            if self.ag_lattice_conductivity_W_m_K > 0.0:
                raise ValueError(
                    "Geometry-corrected contacts currently apply to electron conduction only; "
                    "set Ag lattice conductivity to zero for this comparison"
                )
        if not (300.0 <= initial_lattice_temperature_K <= self.lattice.melting_temperature_K):
            raise ValueError("Initial lattice temperature must be within [300 K, Ag Tm]")
        initial_lattice = np.full(
            self.lattice.temperature_shape, initial_lattice_temperature_K, dtype=float
        )
        self.lattice.enthalpy_density_J_m3 = self.lattice.enthalpy_density_from_state(
            initial_lattice, np.zeros_like(initial_lattice)
        )
        self.lattice.state = self.lattice.state_from_enthalpy_density(
            self.lattice.enthalpy_density_J_m3
        )
        self.parameters = AgParameters.from_csv(electron_table, allow_low_temperature_extension=False)
        self.electron = AgTTMGrid(
            self.parameters,
            self.lattice.ag_fraction,
            tuple(float(v) for v in self.lattice.spacing_A),
            source_power_eV_ps=self._source_eV_ps(),
            initial_temperature_K=self.parameters.Tmin,
        )
        initial_te = np.broadcast_to(
            np.asarray(initial_electron_temperature_K, dtype=float), self.lattice.temperature_shape
        ).copy()
        self.parameters._checked_temperature(initial_te[self.lattice.ag_fraction > 0.0])
        self.electron.temperature_K[self.lattice.ag_fraction > 0.0] = initial_te[
            self.lattice.ag_fraction > 0.0
        ]
        self.time_ps = 0.0
        self.injected_energy_J = 0.0
        self.coupling_energy_J = 0.0
        # Signed local cumulative transfer: positive means electron -> lattice.
        self.coupling_energy_density_J_m3 = np.zeros_like(
            self.lattice.enthalpy_density_J_m3
        )
        self.electron_diffusion_energy_J = 0.0
        self.lattice_diffusion_energy_J = 0.0
        self.initial_total_energy_J = self.total_energy_J
        self.field_snapshots: dict[float, dict[str, np.ndarray]] = {}
        # Diagnostic bound for the actual two-temperature lattice conductivity
        # split. The implicit update does not subcycle to satisfy this limit.
        self.lattice_diffusion_timestep_ps = self._explicit_lattice_stability_limit_ps()

    @staticmethod
    def _finite_volume_laplacian(
        conductivity: np.ndarray,
        spacing: np.ndarray,
        cell_volume: float,
        bottom_sink_temperature_K: float | None = None,
        face_contact_fractions: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None,
    ):
        """Return a symmetric conservative graph Laplacian for x/y-periodic, z-insulated cells.

        The internal-face operator is positive semidefinite with zero row sum.
        An optional z=0 Dirichlet sink adds boundary conductance to the diagonal
        and its fixed-temperature contribution to the returned power vector.
        Thus integrated conductive power is ``-K @ T + boundary_power``.
        """
        k = np.asarray(conductivity, dtype=float)
        shape = k.shape
        if face_contact_fractions is not None:
            if bottom_sink_temperature_K is not None:
                raise ValueError("Face-contact fractions are not defined for a z-boundary sink")
            if len(face_contact_fractions) != 3 or any(
                np.asarray(contact).shape != shape for contact in face_contact_fractions
            ):
                raise ValueError("Each face-contact fraction field must match conductivity shape")
        ncell = int(k.size)
        rows: list[np.ndarray] = []
        cols: list[np.ndarray] = []
        values: list[np.ndarray] = []
        indices = np.arange(ncell, dtype=np.int64).reshape(shape)
        for axis in range(3):
            if shape[axis] <= 1:
                continue
            if axis < 2:
                left = indices
                right = np.roll(indices, -1, axis=axis)
                kr = np.roll(k, -1, axis=axis)
                denom = k + kr
                face_k = np.divide(2.0 * k * kr, denom, out=np.zeros_like(k), where=denom > 0.0)
                conductance = face_k * cell_volume / spacing[axis] ** 2
                if face_contact_fractions is not None:
                    conductance *= face_contact_fractions[axis]
                i = left.ravel()
                j = right.ravel()
                g = conductance.ravel()
            else:
                low = [slice(None)] * 3
                high = [slice(None)] * 3
                low[axis], high[axis] = slice(0, -1), slice(1, None)
                low_t, high_t = tuple(low), tuple(high)
                kl, kh = k[low_t], k[high_t]
                denom = kl + kh
                face_k = np.divide(2.0 * kl * kh, denom, out=np.zeros_like(kl), where=denom > 0.0)
                conductance = face_k * cell_volume / spacing[axis] ** 2
                if face_contact_fractions is not None:
                    conductance *= face_contact_fractions[axis][low_t]
                i = indices[low_t].ravel()
                j = indices[high_t].ravel()
                g = conductance.ravel()
            valid = (g > 0.0) & (i != j)
            i, j, g = i[valid], j[valid], g[valid]
            rows.extend((i, i, j, j))
            cols.extend((i, j, i, j))
            values.extend((g, -g, -g, g))
        if not rows:
            from scipy.sparse import csr_matrix
            operator = csr_matrix((ncell, ncell), dtype=float)
            boundary_power = np.zeros(ncell, dtype=float)
            if bottom_sink_temperature_K is not None and shape[2] > 0:
                bottom_ids = indices[:, :, 0].ravel()
                bottom_k = k[:, :, 0].ravel()
                conductance = 2.0 * bottom_k * cell_volume / spacing[2] ** 2
                valid = conductance > 0.0
                if np.any(valid):
                    ids = bottom_ids[valid]
                    g = conductance[valid]
                    operator = coo_matrix(
                        (g, (ids, ids)), shape=(ncell, ncell)
                    ).tocsr()
                    boundary_power[ids] = g * float(bottom_sink_temperature_K)
            return operator, boundary_power
        operator = coo_matrix(
            (np.concatenate(values), (np.concatenate(rows), np.concatenate(cols))),
            shape=(ncell, ncell),
        ).tocsr()
        boundary_power = np.zeros(ncell, dtype=float)
        if bottom_sink_temperature_K is not None and shape[2] > 0:
            # Cell centres are half a cell from the z=0 isothermal boundary.
            bottom_ids = indices[:, :, 0].ravel()
            bottom_k = k[:, :, 0].ravel()
            conductance = 2.0 * bottom_k * cell_volume / spacing[2] ** 2
            valid = conductance > 0.0
            if np.any(valid):
                ids = bottom_ids[valid]
                g = conductance[valid]
                operator = operator + coo_matrix(
                    (g, (ids, ids)), shape=(ncell, ncell)
                ).tocsr()
                boundary_power[ids] = g * float(bottom_sink_temperature_K)
        return operator, boundary_power

    def _implicit_electron_diffusion(self, dt_ps: float) -> float:
        """Backward-Euler electron diffusion in tabulated internal-energy variables.

        A conservative finite-volume matrix and Newton/CG solve are used. The
        returned value is the net electron-energy change (zero for insulated
        boundaries, up to the nonlinear-solver tolerance).
        """
        active = self.lattice.ag_fraction.ravel() > 0.0
        if not np.any(active):
            return 0.0
        te0 = self.electron.temperature_K.ravel().copy()
        e0 = self.electron.electron_energy_eV.ravel()[active]
        fvol = (
            self.lattice.ag_fraction.ravel()[active]
            * self.electron.cell_volume_A3
        )
        if self.electron_face_contact_fractions is None:
            full_k = self.lattice.ag_fraction * self.parameters.ke(self.electron.temperature_K)
        else:
            # Face areas already encode the exact QSGS Ag contact fraction.
            # Conductivity is therefore intrinsic ke on Ag-bearing cells; using
            # phi*ke here as well would count the volume fraction twice.
            full_k = np.where(
                self.lattice.ag_fraction > 0.0,
                self.parameters.ke(self.electron.temperature_K),
                0.0,
            )
        operator = self._finite_volume_laplacian(
            full_k,
            self.electron.spacing_A,
            self.electron.cell_volume_A3,
            self.bottom_sink_temperature_K,
            self.electron_face_contact_fractions,
        )
        operator_full, boundary_power_full = operator
        operator = operator_full[active][:, active].tocsr()
        boundary_power = boundary_power_full[active]
        if operator.nnz == 0:
            return 0.0
        params = self.parameters
        temperatures = params.temperature_K
        ue_values = params.ue_eV_A3

        def energy_and_capacity(t: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
            u = params.ue(t)
            ix = np.clip(np.searchsorted(temperatures, t, side="right") - 1, 0, len(temperatures) - 2)
            slope = (ue_values[ix + 1] - ue_values[ix]) / (temperatures[ix + 1] - temperatures[ix])
            return fvol * u, fvol * slope

        t = te0[active].copy()
        scale = max(1.0e-12, float(np.max(np.abs(e0))))
        for _ in range(20):
            energy, capacity = energy_and_capacity(t)
            residual = energy - e0 + dt_ps * (operator @ t - boundary_power)
            if float(np.max(np.abs(residual))) <= max(2.0e-13, scale * 2.0e-10):
                break
            jacobian = diags(capacity) + dt_ps * operator
            diagonal = np.asarray(jacobian.diagonal())
            preconditioner = diags(1.0 / np.maximum(diagonal, 1.0e-30))
            delta, info = cg(
                jacobian,
                -residual,
                M=preconditioner,
                rtol=1.0e-10,
                atol=max(1.0e-14, scale * 1.0e-12),
                maxiter=min(3000, max(200, 4 * len(t))),
            )
            if info != 0:
                raise RuntimeError(f"Implicit electron diffusion CG did not converge (info={info})")
            old_norm = float(np.linalg.norm(residual))
            accepted = False
            alpha = 1.0
            while alpha >= 1.0 / 128.0:
                trial = t + alpha * delta
                # Round-off in the sparse solve can put an unchanged 300 K
                # cell a few ulps below the tabulated lower bound. Accept and
                # clamp only that negligible undershoot; larger excursions
                # still trigger a shorter Newton line-search step.
                bound_tolerance_K = 1.0e-9
                if np.any(trial < params.Tmin - bound_tolerance_K) or np.any(
                    trial > params.Tmax + bound_tolerance_K
                ):
                    alpha *= 0.5
                    continue
                trial = np.clip(trial, params.Tmin, params.Tmax)
                trial_energy, _ = energy_and_capacity(trial)
                trial_residual = trial_energy - e0 + dt_ps * (operator @ trial - boundary_power)
                if float(np.linalg.norm(trial_residual)) < old_norm:
                    t = trial
                    accepted = True
                    break
                alpha *= 0.5
            if not accepted:
                raise RuntimeError("Implicit electron diffusion Newton line search failed")
        else:
            raise RuntimeError("Implicit electron diffusion Newton solve did not converge")
        self.electron.temperature_K.ravel()[active] = t
        e1 = self.electron.electron_energy_eV.ravel()[active]
        return float((e1 - e0).sum() * EV_J)

    def _set_electron_energy(self, energy_eV: np.ndarray) -> None:
        active = self.electron.f_ag > 0.0
        fvol = self.electron.f_ag[active] * self.electron.cell_volume_A3
        density = np.asarray(energy_eV, dtype=float)[active] / fvol
        self.electron.temperature_K[active] = self.parameters.temperature_from_ue(density)

    def _active_lattice_temperature_from_enthalpy(self, h_density: np.ndarray, active: np.ndarray) -> np.ndarray:
        """Invert the cell enthalpy law for Ag-bearing cells without full-grid allocation."""
        h = np.asarray(h_density, dtype=float)
        af = self.lattice.ag_fraction.ravel()[active]
        sf = self.lattice.sc_fraction.ravel()[active]
        solid_tm = self.lattice._solid_at_tm_J_m3.ravel()[active]
        latent = self.lattice._latent_volumetric_J_m3.ravel()[active]
        c0 = self.lattice._solid_capacity_at_300_J_m3_K.ravel()[active]
        slope = self.lattice._solid_capacity_slope_J_m3_K2.ravel()[active]
        temp = np.full_like(h, 300.0)
        below = h < 0.0
        temp[below] = 300.0 + h[below] / c0[below]
        solid = (h >= 0.0) & (h < solid_tm) & (c0 > 0.0)
        if np.any(solid):
            target = h[solid]
            c = c0[solid]
            sl = slope[solid]
            delta = target / c
            curved = np.abs(sl) > 1.0e-18
            delta[curved] = 2.0 * target[curved] / (
                c[curved] + np.sqrt(c[curved] ** 2 + 2.0 * sl[curved] * target[curved])
            )
            temp[solid] = 300.0 + delta
        liquidus = solid_tm + latent
        plateau = (h >= solid_tm) & (h <= liquidus) & (latent > 0.0)
        temp[plateau] = self.lattice.melting_temperature_K
        above = (h > liquidus) & (c0 > 0.0)
        if np.any(above):
            d_h = h[above] - liquidus[above]
            c_tm = af[above] * self.lattice.ag.reference_density_kg_m3 * float(
                self.lattice.ag.liquid_cp_J_kg_K[0]
            )
            liquid_slope = np.zeros_like(c_tm)
            if self.lattice.boundary_case == "B":
                sc_cp_tm = float(self.lattice.sc.cp_J_kg_K(self.lattice.melting_temperature_K))
                sc_slope = (
                    self.lattice.sc.cp_mass_J_kg_K[1] - self.lattice.sc.cp_mass_J_kg_K[0]
                ) / (self.lattice.sc.temperature_K[1] - self.lattice.sc.temperature_K[0])
                c_tm += sf[above] * self.lattice.sc.density_kg_m3 * sc_cp_tm
                liquid_slope = sf[above] * self.lattice.sc.density_kg_m3 * sc_slope
            delta = d_h / c_tm
            curved = np.abs(liquid_slope) > 1.0e-18
            delta[curved] = 2.0 * d_h[curved] / (
                c_tm[curved] + np.sqrt(c_tm[curved] ** 2 + 2.0 * liquid_slope[curved] * d_h[curved])
            )
            temp[above] = self.lattice.melting_temperature_K + delta
        if self.lattice.boundary_case == "B" and np.any(
            temp > self.lattice.sc.temperature_K[-1] + 1.0e-9
        ):
            raise ValueError("Ti3SiC2 candidate table exceeds 1473.15 K")
        return temp

    def _implicit_local_coupling(self, dt_ps: float) -> float:
        """Solve each Ag-cell electron/lattice exchange implicitly and conservatively."""
        if dt_ps <= 0.0:
            return 0.0
        active_mask = self.lattice.ag_fraction > 0.0
        if not np.any(active_mask):
            return 0.0
        active = active_mask.ravel()
        fvol = self.lattice.ag_fraction.ravel()[active] * self.electron.cell_volume_A3
        te0 = self.electron.temperature_K.ravel()[active].copy()
        ee0 = self.electron.electron_energy_eV.ravel()[active]
        h0 = self.lattice.enthalpy_density_J_m3.ravel()[active] * self.lattice.cell_volume_m3
        total_J = ee0 * EV_J + h0
        tl0 = self.lattice.state.temperature_K.ravel()[active]
        low = np.minimum(te0, tl0)
        high = np.maximum(te0, tl0)
        active_ids = np.flatnonzero(active)

        def residual_at(te: np.ndarray) -> np.ndarray:
            ee = fvol * self.parameters.ue(te)
            h_trial = (total_J - ee * EV_J) / self.lattice.cell_volume_m3
            tl = self._active_lattice_temperature_from_enthalpy(h_trial, active)
            return ee - ee0 + dt_ps * fvol * self.parameters.g(te) * (te - tl)

        f_low = residual_at(low)
        f_high = residual_at(high)
        if np.any(f_low > 1.0e-11) or np.any(f_high < -1.0e-11):
            raise RuntimeError("Implicit electron-lattice coupling root is not bracketed")
        for _ in range(40):
            mid = 0.5 * (low + high)
            f_mid = residual_at(mid)
            positive = f_mid > 0.0
            high = np.where(positive, mid, high)
            low = np.where(positive, low, mid)
            if float(np.max(high - low)) <= 1.0e-5:
                break
        te1 = 0.5 * (low + high)
        ee1 = fvol * self.parameters.ue(te1)
        h1 = total_J - ee1 * EV_J
        self.electron.temperature_K.ravel()[active] = te1
        h_full = self.lattice.enthalpy_density_J_m3.ravel().copy()
        h_full[active_ids] = h1 / self.lattice.cell_volume_m3
        self.lattice.enthalpy_density_J_m3 = h_full.reshape(self.lattice.temperature_shape)
        local_transfer_J = (ee0 - ee1) * EV_J
        self.coupling_energy_density_J_m3.ravel()[active_ids] += (
            local_transfer_J / self.lattice.cell_volume_m3
        )
        self.lattice.state = self.lattice.state_from_enthalpy_density(
            self.lattice.enthalpy_density_J_m3
        )
        return float(local_transfer_J.sum())

    def _lattice_conductivity(self) -> np.ndarray:
        state = self.lattice.state
        k = self.lattice.ag_fraction * self.ag_lattice_conductivity_W_m_K
        if self.lattice.boundary_case == "B":
            k = k + self.lattice.sc_fraction * np.asarray(
                self.lattice.sc.conductivity_W_m_K(state.temperature_K)
            )
        return np.asarray(k, dtype=float)

    def _explicit_lattice_stability_limit_ps(self) -> float:
        k = self._lattice_conductivity()
        capacity = self.lattice.sensible_capacity_J_m3_K()
        inverse_spacing2 = 1.0 / (self.lattice.spacing_A * 1.0e-10) ** 2
        face_sum = np.zeros(self.lattice.temperature_shape, dtype=float)
        for axis in range(3):
            if axis < 2:
                right_k = np.roll(k, -1, axis=axis)
                left_k = np.roll(k, 1, axis=axis)
            else:
                right_k = np.concatenate((k[:, :, 1:], np.zeros_like(k[:, :, :1])), axis=2)
                left_k = np.concatenate((np.zeros_like(k[:, :, :1]), k[:, :, :-1]), axis=2)
            right = np.divide(2.0 * k * right_k, k + right_k, out=np.zeros_like(k), where=(k + right_k) > 0.0)
            left = np.divide(2.0 * k * left_k, k + left_k, out=np.zeros_like(k), where=(k + left_k) > 0.0)
            face_sum += (right + left) * inverse_spacing2[axis]
        if self.bottom_sink_temperature_K is not None:
            face_sum[:, :, 0] += 2.0 * k[:, :, 0] * inverse_spacing2[2]
        rate = np.divide(face_sum, capacity, out=np.zeros_like(face_sum), where=capacity > 0.0)
        max_rate = float(rate.max())
        return np.inf if max_rate <= 0.0 else 0.45 / max_rate * 1.0e12

    def _implicit_lattice_diffusion(self, dt_s: float) -> float:
        """Backward-Euler enthalpy diffusion with a conservative face operator."""
        if dt_s <= 0.0:
            return 0.0
        volume = self.lattice.cell_volume_m3
        h0 = self.lattice.enthalpy_density_J_m3.ravel().copy() * volume
        k = self._lattice_conductivity()
        operator = self._finite_volume_laplacian(
            k,
            self.lattice.spacing_A * 1.0e-10,
            volume,
            self.bottom_sink_temperature_K,
        )
        operator, boundary_power = operator
        if operator.nnz == 0:
            return 0.0
        h = h0.copy()
        n = h.size
        identity = eye(n, format="csr")
        scale = max(1.0e-24, float(np.max(np.abs(h0))))

        def evaluate(h_cell: np.ndarray):
            h_density = (h_cell / volume).reshape(self.lattice.temperature_shape)
            state = self.lattice.state_from_enthalpy_density(h_density)
            temp = state.temperature_K.ravel()
            residual = h_cell - h0 + dt_s * (operator @ temp - boundary_power)
            capacity = self.lattice.sensible_capacity_J_m3_K(state).ravel() * volume
            solid_tm = self.lattice._solid_at_tm_J_m3.ravel()
            latent = self.lattice._latent_volumetric_J_m3.ravel()
            h_density_flat = h_cell / volume
            plateau = (latent > 0.0) & (h_density_flat > solid_tm) & (h_density_flat < solid_tm + latent)
            dtdh = np.divide(1.0, capacity, out=np.zeros_like(capacity), where=capacity > 0.0)
            dtdh[plateau] = 0.0
            return residual, temp, dtdh

        for _ in range(20):
            residual, temp, dtdh = evaluate(h)
            if float(np.max(np.abs(residual))) <= max(1.0e-30, scale * 1.0e-10):
                break
            jacobian = identity + dt_s * (operator @ diags(dtdh))
            diagonal = np.maximum(np.abs(jacobian.diagonal()), 1.0e-12)
            preconditioner = diags(1.0 / diagonal)
            delta, info = gmres(
                jacobian,
                -residual,
                M=preconditioner,
                rtol=1.0e-9,
                atol=max(1.0e-32, scale * 1.0e-12),
                restart=50,
                maxiter=300,
            )
            if info != 0:
                raise RuntimeError(f"Implicit lattice diffusion GMRES did not converge (info={info})")
            old_norm = float(np.linalg.norm(residual))
            accepted = False
            alpha = 1.0
            while alpha >= 1.0 / 128.0:
                trial = h + alpha * delta
                try:
                    trial_residual, _, _ = evaluate(trial)
                except ValueError:
                    alpha *= 0.5
                    continue
                if float(np.linalg.norm(trial_residual)) < old_norm:
                    h = trial
                    accepted = True
                    break
                alpha *= 0.5
            if not accepted:
                raise RuntimeError("Implicit lattice diffusion Newton line search failed")
        else:
            raise RuntimeError("Implicit lattice diffusion Newton solve did not converge")
        self.lattice.enthalpy_density_J_m3 = (h / volume).reshape(self.lattice.temperature_shape)
        self.lattice.state = self.lattice.state_from_enthalpy_density(
            self.lattice.enthalpy_density_J_m3
        )
        return float(h.sum() - h0.sum())

    def _source_eV_ps(self) -> np.ndarray:
        return self.lattice.source_W_m3 * self.lattice.cell_volume_m3 * 1.0e-12 / EV_J

    @property
    def actual_source_power_W(self) -> float:
        return float(self.electron.source_power.sum() * EV_J / 1.0e-12)

    @property
    def electron_energy_J(self) -> float:
        return float(self.electron.electron_energy_eV.sum() * EV_J)

    @property
    def lattice_energy_J(self) -> float:
        return float(self.lattice.enthalpy_density_J_m3.sum() * self.lattice.cell_volume_m3)

    @property
    def total_energy_J(self) -> float:
        return self.electron_energy_J + self.lattice_energy_J

    def advance(
        self,
        dt_ps: float,
        heat_on: bool = True,
        coupling_on: bool = True,
    ) -> dict[str, float]:
        if not np.isfinite(dt_ps) or dt_ps <= 0.0:
            raise ValueError("dt_ps must be finite and positive")
        initial_total = self.total_energy_J
        source_rate = self.electron.source_power if heat_on else np.zeros_like(self.electron.source_power)
        source_half_eV = source_rate * (0.5 * dt_ps)
        source_increment_J = float(source_rate.sum() * dt_ps * EV_J)
        diffusion_increment_J = 0.0
        lattice_diffusion_increment_J = 0.0
        coupling_increment_J = 0.0

        # Symmetric operator splitting: source/diffusion half steps surround
        # a fully implicit, cell-conservative electron-lattice exchange step.
        self._set_electron_energy(self.electron.electron_energy_eV + source_half_eV)
        diffusion_increment_J += self._implicit_electron_diffusion(0.5 * dt_ps)
        lattice_diffusion_increment_J += self._implicit_lattice_diffusion(0.5 * dt_ps * 1.0e-12)
        if coupling_on:
            coupling_increment_J = self._implicit_local_coupling(dt_ps)
        lattice_diffusion_increment_J += self._implicit_lattice_diffusion(0.5 * dt_ps * 1.0e-12)
        diffusion_increment_J += self._implicit_electron_diffusion(0.5 * dt_ps)
        self._set_electron_energy(self.electron.electron_energy_eV + source_half_eV)
        self.time_ps += dt_ps
        self.electron.time_ps += dt_ps
        substeps = 1

        self.injected_energy_J += source_increment_J
        self.coupling_energy_J += coupling_increment_J
        self.electron_diffusion_energy_J += diffusion_increment_J
        self.lattice_diffusion_energy_J += lattice_diffusion_increment_J
        final_total = self.total_energy_J
        energy_residual = (
            final_total
            - initial_total
            - source_increment_J
            - diffusion_increment_J
            - lattice_diffusion_increment_J
        )
        active = self.lattice.ag_fraction > 0.0
        ag_weight_sum = float(self.lattice.ag_fraction.sum())
        max_electron_temperature = (
            float(self.electron.temperature_K[active].max()) if np.any(active) else self.parameters.Tmin
        )
        mean_electron_temperature = (
            float(np.sum(self.lattice.ag_fraction * self.electron.temperature_K) / ag_weight_sum)
            if ag_weight_sum > 0.0 else self.parameters.Tmin
        )
        return {
            "time_ps": self.time_ps,
            "dt_ps": dt_ps,
            "electron_substeps": float(substeps),
            "max_electron_temperature_K": max_electron_temperature,
            "mean_electron_temperature_K": mean_electron_temperature,
            "max_lattice_temperature_K": float(self.lattice.state.temperature_K.max()),
            "mean_ag_lattice_temperature_K": self.lattice._ag_weighted_mean_temperature(),
            "max_ag_liquid_fraction": float(self.lattice.state.ag_liquid_fraction.max()),
            "ag_liquid_volume_fraction": self.lattice._liquid_ag_volume_fraction(),
            "electron_energy_J": self.electron_energy_J,
            "lattice_enthalpy_J": self.lattice_energy_J,
            "total_energy_J": final_total,
            "source_energy_J": source_increment_J,
            "electron_to_lattice_energy_J": coupling_increment_J,
            "electron_diffusion_net_J": diffusion_increment_J,
            "lattice_diffusion_net_J": lattice_diffusion_increment_J,
            "energy_balance_residual_J": energy_residual,
        }

    def run(
        self,
        duration_ps: float,
        dt_ps: float,
        heat_on: bool,
        sample_interval_ps: float | None = None,
        snapshot_times_ps: tuple[float, ...] | list[float] = (),
    ) -> list[dict[str, float]]:
        if duration_ps <= 0.0 or dt_ps <= 0.0:
            raise ValueError("duration_ps and dt_ps must be positive")
        sample_interval_ps = dt_ps if sample_interval_ps is None else sample_interval_ps
        if sample_interval_ps <= 0.0:
            raise ValueError("sample_interval_ps must be positive")
        end_ps = self.time_ps + duration_ps
        snapshots = sorted(float(value) for value in snapshot_times_ps)
        tolerance = max(1.0e-12, dt_ps * 1.0e-9)
        if any(not np.isfinite(value) for value in snapshots):
            raise ValueError("snapshot_times_ps must contain only finite values")
        if any(value < self.time_ps - tolerance or value > end_ps + tolerance for value in snapshots):
            raise ValueError("Snapshot times must lie within this run interval")

        def capture_snapshot(target_ps: float) -> None:
            self.field_snapshots[float(target_ps)] = {
                "electron_temperature_K": self.electron.temperature_K.copy(),
                "lattice_temperature_K": self.lattice.state.temperature_K.copy(),
                "ag_liquid_fraction": self.lattice.state.ag_liquid_fraction.copy(),
                "electron_energy_density_eV_A3": (
                    self.electron.electron_energy_eV / self.electron.cell_volume_A3
                ),
                "lattice_enthalpy_density_J_m3": self.lattice.enthalpy_density_J_m3.copy(),
                "electron_to_lattice_cumulative_energy_density_J_m3": (
                    self.coupling_energy_density_J_m3.copy()
                ),
            }

        while snapshots and abs(snapshots[0] - self.time_ps) <= tolerance:
            capture_snapshot(snapshots.pop(0))
        next_sample = self.time_ps + sample_interval_ps
        sample_start_ps = self.time_ps
        accumulated: dict[str, float] = {
            "electron_substeps": 0.0,
            "source_energy_J": 0.0,
            "electron_to_lattice_energy_J": 0.0,
            "electron_diffusion_net_J": 0.0,
            "lattice_diffusion_net_J": 0.0,
            "energy_balance_residual_J": 0.0,
        }
        rows: list[dict[str, float]] = []
        while self.time_ps < end_ps - max(1.0e-15, duration_ps * 1.0e-13):
            step_limits = [dt_ps, end_ps - self.time_ps, next_sample - self.time_ps]
            if snapshots:
                step_limits.append(snapshots[0] - self.time_ps)
            step = min(step_limits)
            if step <= tolerance * 1.0e-6:
                raise RuntimeError("Snapshot schedule produced a non-positive integration step")
            row = self.advance(step, heat_on=heat_on)
            while snapshots and abs(self.time_ps - snapshots[0]) <= tolerance:
                capture_snapshot(snapshots.pop(0))
            for name in accumulated:
                accumulated[name] += row[name]
            if self.time_ps >= next_sample - max(1.0e-15, duration_ps * 1.0e-13) or self.time_ps >= end_ps - 1.0e-15:
                sampled = dict(row)
                sampled.update(accumulated)
                sampled["dt_ps"] = self.time_ps - sample_start_ps
                rows.append(sampled)
                accumulated = {name: 0.0 for name in accumulated}
                sample_start_ps = self.time_ps
                while next_sample <= self.time_ps + max(1.0e-15, duration_ps * 1.0e-13):
                    next_sample += sample_interval_ps
        return rows


def write_csv(path: Path, rows: list[dict[str, float]]) -> None:
    if not rows:
        raise ValueError("Cannot write an empty time history")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, quantile: float) -> float:
    valid = (weights > 0.0) & np.isfinite(values)
    v = np.asarray(values[valid], dtype=float)
    w = np.asarray(weights[valid], dtype=float)
    if v.size == 0:
        return 300.0
    order = np.argsort(v)
    cumulative = np.cumsum(w[order])
    index = min(int(np.searchsorted(cumulative, quantile * cumulative[-1])), v.size - 1)
    return float(v[order[index]])


def field_snapshot_metrics(
    time_ps: float,
    fields: dict[str, np.ndarray],
    ag_fraction: np.ndarray,
    spacing_A: np.ndarray,
    size_A: np.ndarray,
    roi_xyz: tuple[slice, slice, slice] | None = None,
    gaussian_radius_A: float = 30.0,
    absorption_depth_A: float = 20.0,
) -> dict[str, float]:
    """Summarize Ag-weighted temperatures and half-rise hotspot geometry."""
    ag = np.asarray(ag_fraction, dtype=float)
    active = ag >= PHASE_FRACTION_ROUNDOFF_TOL
    total_ag_weight = float(ag[active].sum())
    x = (np.arange(ag.shape[0]) + 0.5) * spacing_A[0]
    y = (np.arange(ag.shape[1]) + 0.5) * spacing_A[1]
    z = (np.arange(ag.shape[2]) + 0.5) * spacing_A[2]
    cx, cy = size_A[0] / 2.0, size_A[1] / 2.0
    xrel = (x - cx + size_A[0] / 2.0) % size_A[0] - size_A[0] / 2.0
    yrel = (y - cy + size_A[1] / 2.0) % size_A[1] - size_A[1] / 2.0
    radius2 = xrel[:, None, None] ** 2 + yrel[None, :, None] ** 2
    depth = size_A[2] - z[None, None, :]
    source_zone = (radius2 <= gaussian_radius_A**2) & (depth <= absorption_depth_A)
    summary: dict[str, float] = {"time_ps": float(time_ps)}
    for label, key in (("electron", "electron_temperature_K"), ("lattice", "lattice_temperature_K")):
        temperature = np.asarray(fields[key], dtype=float)
        peak = float(temperature[active].max()) if np.any(active) else 300.0
        mean = (
            float(np.sum(ag[active] * temperature[active]) / total_ag_weight)
            if total_ag_weight > 0.0
            else 300.0
        )
        threshold = 300.0 + 0.5 * max(peak - 300.0, 0.0)
        hot = active & (temperature >= threshold) if peak > 300.0 else np.zeros_like(active)
        hot_weight = float(ag[hot].sum())
        if hot_weight > 0.0:
            radius_rms = float(np.sqrt(np.sum(ag * radius2 * hot) / hot_weight))
            hot_depth = np.broadcast_to(depth, ag.shape)[hot]
            depth_min = float(hot_depth.min())
            depth_max = float(hot_depth.max())
            depth_extent = depth_max - depth_min
        else:
            radius_rms = depth_min = depth_max = depth_extent = 0.0
        summary[f"{label}_peak_K"] = peak
        summary[f"{label}_ag_weighted_mean_K"] = mean
        for quantile in (0.5, 0.9, 0.99, 0.999):
            summary[f"{label}_ag_weighted_q{quantile:g}_K"] = _weighted_quantile(
                temperature, ag, quantile
            )
        for delta_K in (1.0, 5.0, 10.0, 50.0, 100.0):
            raised = active & (temperature >= 300.0 + delta_K)
            summary[f"{label}_ag_volume_fraction_above_300_plus_{delta_K:g}K"] = (
                float(ag[raised].sum() / total_ag_weight) if total_ag_weight > 0.0 else 0.0
            )
        summary[f"{label}_half_rise_threshold_K"] = threshold
        summary[f"{label}_half_rise_hotspot_diameter_2rms_A"] = 2.0 * radius_rms
        summary[f"{label}_half_rise_hotspot_depth_extent_A"] = depth_extent
        summary[f"{label}_half_rise_hotspot_min_depth_A"] = depth_min
        summary[f"{label}_half_rise_hotspot_max_depth_A"] = depth_max
        summary[f"{label}_half_rise_hotspot_ag_volume_fraction"] = (
            hot_weight / total_ag_weight if total_ag_weight > 0.0 else 0.0
        )
        outside = active & ~source_zone
        outside_weight = float(ag[outside].sum())
        summary[f"{label}_outside_source_ag_weighted_mean_K"] = (
            float(np.sum(ag[outside] * temperature[outside]) / outside_weight)
            if outside_weight > 0.0
            else 300.0
        )
        summary[f"{label}_outside_source_ag_volume_fraction_above_5K_rise"] = (
            float(ag[outside & (temperature >= 305.0)].sum() / outside_weight)
            if outside_weight > 0.0
            else 0.0
        )
        if roi_xyz is not None:
            roi_ag = ag[roi_xyz]
            roi_temp = temperature[roi_xyz]
            roi_total = float(roi_ag.sum())
            roi_active = roi_ag >= PHASE_FRACTION_ROUNDOFF_TOL
            roi_hot = roi_active & (roi_temp >= threshold) if peak > 300.0 else np.zeros_like(roi_active)
            summary[f"{label}_central_tile_half_rise_hotspot_ag_volume_fraction"] = (
                float(roi_ag[roi_hot].sum() / roi_total) if roi_total > 0.0 else 0.0
            )
    liquid = np.asarray(fields["ag_liquid_fraction"], dtype=float)
    summary["max_ag_liquid_fraction"] = float(liquid[active].max()) if np.any(active) else 0.0
    summary["ag_liquid_volume_fraction"] = (
        float(np.sum(ag[active] * liquid[active]) / total_ag_weight)
        if total_ag_weight > 0.0
        else 0.0
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spacing", choices=("5A", "10A"), default="10A")
    parser.add_argument("--boundary", choices=("A", "B"), required=True)
    parser.add_argument("--heat-ps", type=float, default=0.01)
    parser.add_argument("--cool-ps", type=float, default=0.01)
    parser.add_argument("--dt-fs", type=float, default=0.5,
                        help="outer operator-splitting step; diffusion and local coupling are implicit")
    parser.add_argument("--sample-ps", type=float, default=0.005)
    parser.add_argument("--bottom-sink", action="store_true",
                        help="apply a 300 K ideal isothermal sink at the bottom z face")
    parser.add_argument("--face-contacts", type=Path, default=None,
                        help="use exact QSGS Ag contact fractions for electron conduction links")
    parser.add_argument("--source-geometry", type=Path,
                        default=DATA_DIR.parent / "qsgs_morphology.npz",
                        help="binary QSGS geometry used to integrate heat input over actual Ag voxels")
    parser.add_argument("--tile-xy", type=int, default=1,
                        help="repeat the projected QSGS geometry this many times in periodic x/y (odd values center the source on a tile)")
    parser.add_argument("--snapshot-ps", default="",
                        help="comma-separated heating times at which to save 3D Te, Tl, and Ag liquid fields")
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    if args.heat_ps <= 0.0 or args.cool_ps < 0.0:
        parser.error("--heat-ps must be positive and --cool-ps must be non-negative")
    try:
        snapshot_times_ps = tuple(
            float(value.strip()) for value in args.snapshot_ps.split(",") if value.strip()
        )
    except ValueError:
        parser.error("--snapshot-ps must be a comma-separated list of numbers")
    if any(value < 0.0 or value > args.heat_ps for value in snapshot_times_ps):
        parser.error("--snapshot-ps values must lie within the heating interval")
    projected = DATA_DIR / f"qsgs_phase_fractions_{args.spacing}.npz"
    model = QSGSDoubleTemperatureGrid(
        projected,
        args.boundary,
        bottom_sink_temperature_K=300.0 if args.bottom_sink else None,
        face_contacts_npz=args.face_contacts,
        tile_xy=args.tile_xy,
        source_geometry_npz=args.source_geometry,
    )
    dt_ps = args.dt_fs * 1.0e-3
    rows = model.run(
        args.heat_ps,
        dt_ps,
        heat_on=True,
        sample_interval_ps=args.sample_ps,
        snapshot_times_ps=snapshot_times_ps,
    )
    if args.cool_ps > 0.0:
        rows += model.run(args.cool_ps, dt_ps, heat_on=False, sample_interval_ps=args.sample_ps)
    output = args.output or DATA_DIR / f"qsgs_ttm_enthalpy_{args.spacing}_BC{args.boundary}.csv"
    write_csv(output, rows)
    field_output = None
    snapshot_metrics = []
    if model.field_snapshots:
        field_output = output.with_name(f"{output.stem}_fields.npz")
        arrays: dict[str, np.ndarray] = {
            "ag_fraction": model.lattice.ag_fraction,
            "skeleton_fraction": model.lattice.sc_fraction,
            "electron_source_power_W_per_cell": model.lattice.source_cell_power_W,
            "electron_source_density_W_m3": model.lattice.source_W_m3,
            "spacing_A": model.lattice.spacing_A,
            "size_A": model.lattice.size_A,
            "x_centers_A": (np.arange(model.lattice.temperature_shape[0]) + 0.5) * model.lattice.spacing_A[0],
            "y_centers_A": (np.arange(model.lattice.temperature_shape[1]) + 0.5) * model.lattice.spacing_A[1],
            "z_centers_A": (np.arange(model.lattice.temperature_shape[2]) + 0.5) * model.lattice.spacing_A[2],
        }
        for snapshot_time, fields in sorted(model.field_snapshots.items()):
            label = f"t_{snapshot_time:.9g}ps"
            for name, field in fields.items():
                arrays[f"{label}_{name}"] = field
            roi = None
            if model.lattice.tile_xy > 1:
                nx, ny = np.load(projected, allow_pickle=False)["ag_fraction"].shape[:2]
                center = model.lattice.tile_xy // 2
                roi = (slice(center * nx, (center + 1) * nx), slice(center * ny, (center + 1) * ny), slice(None))
            snapshot_metrics.append(
                field_snapshot_metrics(
                    snapshot_time,
                    fields,
                    model.lattice.ag_fraction,
                    model.lattice.spacing_A,
                    model.lattice.size_A,
                    roi_xyz=roi,
                    gaussian_radius_A=model.lattice.gaussian_radius_A,
                    absorption_depth_A=model.lattice.absorption_depth_A,
                )
            )
        np.savez_compressed(field_output, **arrays)
    summary = {
        "status": "CONTINUUM_TWO_TEMPERATURE_ENTHALPY_SHORT_PILOT",
        "integration": "conservative semi-implicit finite-volume diffusion with implicit local electron-lattice exchange",
        "ag_lattice_conductivity_W_m_K": model.ag_lattice_conductivity_W_m_K,
        "boundary_case": args.boundary,
        "bottom_sink_temperature_K": model.bottom_sink_temperature_K,
        "grid_spacing_label": args.spacing,
        "dt_fs": args.dt_fs,
        "grid_shape_xyz": list(model.lattice.temperature_shape),
        "domain_size_A": model.lattice.size_A.tolist(),
        "xy_tile_count": model.lattice.tile_xy,
        "source_peak_flux_W_m2": model.lattice.peak_flux_W_m2,
        "source_gaussian_radius_A": model.lattice.gaussian_radius_A,
        "source_absorption_depth_A": model.lattice.absorption_depth_A,
        "electron_conduction_geometry_correction": (
            "exact source-QSGS Ag face-contact area fractions" if model.electron_face_contact_fractions is not None else "legacy volume-fraction harmonic links"
        ),
        "face_contacts_npz": model.face_contacts_npz.name if model.face_contacts_npz is not None else None,
        "source_geometry_npz": model.lattice.source_geometry_npz.name if model.lattice.source_geometry_npz is not None else None,
        "source_integration_method": model.lattice.source_integration_method,
        "source_power_in_roundoff_only_cells_W": model.lattice.source_power_in_roundoff_only_cells_W,
        "source_power_grid_integral_W": model.actual_source_power_W,
        "heat_duration_ps": args.heat_ps,
        "cool_duration_ps": args.cool_ps,
        "total_time_ps": model.time_ps,
        "electron_table_range_K": [float(model.parameters.temperature_K[0]), float(model.parameters.temperature_K[-1])],
        "electron_energy_final_J": model.electron_energy_J,
        "lattice_enthalpy_final_J": model.lattice_energy_J,
        "total_energy_final_J": model.total_energy_J,
        "source_energy_total_J": model.injected_energy_J,
        "electron_to_lattice_cumulative_J": model.coupling_energy_J,
        "electron_diffusion_net_cumulative_J": model.electron_diffusion_energy_J,
        "lattice_diffusion_net_cumulative_J": model.lattice_diffusion_energy_J,
        "energy_balance_residual_J": (
            model.total_energy_J
            - model.initial_total_energy_J
            - model.injected_energy_J
            - model.electron_diffusion_energy_J
            - model.lattice_diffusion_energy_J
        ),
        "electron_temperature_max_K": float(model.electron.temperature_K[model.lattice.ag_fraction > 0.0].max()),
        "lattice_temperature_max_K": float(model.lattice.state.temperature_K.max()),
        "ag_liquid_volume_fraction": model.lattice._liquid_ag_volume_fraction(),
        "electron_temperature_peak_over_run_K": max(
            row["max_electron_temperature_K"] for row in rows
        ),
        "lattice_temperature_peak_over_run_K": max(
            row["max_lattice_temperature_K"] for row in rows
        ),
        "max_ag_liquid_fraction_over_run": max(
            row["max_ag_liquid_fraction"] for row in rows
        ),
        "limitations": [
            "This is a continuum two-temperature enthalpy reference solver, not atomistic TTM-MD.",
            "Ag electron properties are cross-paper candidate fits over 300–20000 K; low-temperature extension is disabled.",
            "Ti3SiC2 boundary cases A/B are limiting assumptions, not measured interface heat transfer.",
            "x/y are periodic; z is insulated unless --bottom-sink requests an ideal 300 K bottom bath.",
            "No stable-pool controller is enabled in this short numerical pilot.",
        ],
        "time_history_csv": output.name,
        "field_snapshots_npz": field_output.name if field_output is not None else None,
        "field_snapshot_times_ps": sorted(model.field_snapshots),
        "field_snapshot_spatial_metrics": snapshot_metrics,
    }
    output.with_suffix(".json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
