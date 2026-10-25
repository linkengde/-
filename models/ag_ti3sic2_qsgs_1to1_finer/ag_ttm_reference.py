"""Energy-based reference solver for an Ag-only TTM grid.

This is a numerical reference/prototype, not a LAMMPS fix. It uses metal-unit
electron properties, x/y periodic finite-volume conduction, zero-flux z faces,
and a deterministic mean-field electron-lattice exchange. The caller applies
the returned lattice energy to Newton-layer Ag velocities with
``rescale_lattice_kinetic``. This exact energy operator omits Langevin noise.
"""

from __future__ import annotations

import csv
import itertools
from dataclasses import dataclass
from pathlib import Path

import numpy as np


KB_EV_K = 8.617333262145e-5
M_VV2E_METAL = 1.0364269e-4  # eV / (g mol^-1 * (Angstrom/ps)^2)


class TTMRangeError(ValueError):
    """The electronic temperature or energy left the tabulated range."""


@dataclass(frozen=True)
class AgParameters:
    temperature_K: np.ndarray
    ce_eV_A3_K: np.ndarray
    g_eV_ps_A3_K: np.ndarray
    ke_eV_ps_A_K: np.ndarray
    ue_eV_A3: np.ndarray
    allow_low_temperature_extension: bool = False

    @classmethod
    def from_csv(
        cls,
        path: str | Path,
        allow_low_temperature_extension: bool = False,
    ) -> "AgParameters":
        columns: dict[str, list[float]] = {
            "Te_K": [],
            "Ce_eV_A-3_K-1": [],
            "G_eV_ps-1_A-3_K-1": [],
            "ke_eV_ps-1_A-1_K-1": [],
            "Ue_delta_from_300K_eV_A-3": [],
        }
        with Path(path).open("r", newline="", encoding="utf-8-sig") as stream:
            reader = csv.DictReader(stream)
            missing = set(columns) - set(reader.fieldnames or ())
            if missing:
                raise ValueError(f"Missing Ag parameter columns: {sorted(missing)}")
            for row in reader:
                for name in columns:
                    columns[name].append(float(row[name]))

        temp, ce, g, ke, ue = (np.asarray(columns[k], dtype=float) for k in columns)
        if temp.size < 2 or not np.all(np.diff(temp) > 0.0):
            raise ValueError("Ag parameter temperatures must strictly increase")
        if not (np.all(np.isfinite(ce)) and np.all(np.isfinite(g)) and np.all(np.isfinite(ke))):
            raise ValueError("Ag Ce, G, and ke must be finite")
        if np.any(ce <= 0.0) or np.any(g <= 0.0) or np.any(ke <= 0.0):
            raise ValueError("Ag Ce, G, and ke must be positive")
        if not np.all(np.isfinite(ue)) or not np.all(np.diff(ue) > 0.0):
            raise ValueError("Ag electron internal energy must be finite and strictly increasing")
        return cls(temp, ce, g, ke, ue, allow_low_temperature_extension)

    @property
    def Tmin(self) -> float:
        return 1.0 if self.allow_low_temperature_extension else float(self.temperature_K[0])

    @property
    def Tmax(self) -> float:
        return float(self.temperature_K[-1])

    def _checked_temperature(self, temperature_K: np.ndarray | float) -> np.ndarray:
        t = np.asarray(temperature_K, dtype=float)
        if not np.all(np.isfinite(t)):
            raise TTMRangeError("Non-finite electronic temperature")
        if np.any(t < self.Tmin) or np.any(t > self.Tmax):
            lo = float(np.min(t))
            hi = float(np.max(t))
            raise TTMRangeError(
                f"Electronic temperature [{lo:g}, {hi:g}] K outside table "
                f"[{self.Tmin:g}, {self.Tmax:g}] K; extrapolation is disabled"
            )
        return t

    def _interpolate(self, values: np.ndarray, temperature_K: np.ndarray | float) -> np.ndarray:
        t = self._checked_temperature(temperature_K)
        return np.interp(np.maximum(t, self.temperature_K[0]), self.temperature_K, values)

    def ce(self, temperature_K: np.ndarray | float) -> np.ndarray:
        t = self._checked_temperature(temperature_K)
        value = self._interpolate(self.ce_eV_A3_K, t)
        if self.allow_low_temperature_extension:
            value = np.where(t < self.temperature_K[0], self.ce_eV_A3_K[0] * t / self.temperature_K[0], value)
        return value

    def g(self, temperature_K: np.ndarray | float) -> np.ndarray:
        t = self._checked_temperature(temperature_K)
        value = self._interpolate(self.g_eV_ps_A3_K, t)
        if self.allow_low_temperature_extension:
            value = np.where(t < self.temperature_K[0], self.g_eV_ps_A3_K[0], value)
        return value

    def ke(self, temperature_K: np.ndarray | float) -> np.ndarray:
        t = self._checked_temperature(temperature_K)
        value = self._interpolate(self.ke_eV_ps_A_K, t)
        if self.allow_low_temperature_extension:
            value = np.where(t < self.temperature_K[0], self.ke_eV_ps_A_K[0], value)
        return value

    def ue(self, temperature_K: np.ndarray | float) -> np.ndarray:
        """Tabulated electron energy density relative to 300 K."""
        t = self._checked_temperature(temperature_K)
        value = self._interpolate(self.ue_eV_A3, t)
        if self.allow_low_temperature_extension:
            t0 = self.temperature_K[0]
            gamma = self.ce_eV_A3_K[0] / t0
            low_value = self.ue_eV_A3[0] + 0.5 * gamma * (t * t - t0 * t0)
            value = np.where(t < t0, low_value, value)
        return value

    def temperature_from_ue(self, energy_density_eV_A3: np.ndarray | float) -> np.ndarray:
        energy = np.asarray(energy_density_eV_A3, dtype=float)
        if not np.all(np.isfinite(energy)):
            raise TTMRangeError("Non-finite electron energy density")
        high = self.ue_eV_A3[-1]
        if self.allow_low_temperature_extension:
            t0 = self.temperature_K[0]
            gamma = self.ce_eV_A3_K[0] / t0
            low = self.ue_eV_A3[0] + 0.5 * gamma * (self.Tmin**2 - t0**2)
        else:
            low = self.ue_eV_A3[0]
        if np.any(energy < low - 1.0e-14) or np.any(energy > high + 1.0e-14):
            raise TTMRangeError(
                f"Electron energy density outside [{low:g}, {high:g}] eV/A^3"
            )
        clipped = np.clip(energy, low, high)
        if self.allow_low_temperature_extension:
            t0 = self.temperature_K[0]
            gamma = self.ce_eV_A3_K[0] / t0
            low_region = clipped < self.ue_eV_A3[0]
            temperature = np.interp(
                np.maximum(clipped, self.ue_eV_A3[0]), self.ue_eV_A3, self.temperature_K
            )
            low_temperature = np.sqrt(
                np.maximum(0.0, t0 * t0 + 2.0 * (clipped - self.ue_eV_A3[0]) / gamma)
            )
            return np.where(low_region, low_temperature, temperature)
        return np.interp(clipped, self.ue_eV_A3, self.temperature_K)


@dataclass(frozen=True)
class StepAccounting:
    elapsed_ps: float
    source_to_electron_eV: float
    electron_to_lattice_eV: float
    diffusion_net_eV: float
    conservation_residual_eV: float
    substeps: int
    lattice_energy_per_cell_eV: np.ndarray


@dataclass(frozen=True)
class AtomGridMap:
    """Conservative atom-to-cell interpolation weights, flattened cell IDs."""

    cell_ids: np.ndarray  # shape (N_atoms, N_weights)
    weights: np.ndarray  # same shape; each atom's weights sum to one
    grid_shape: tuple[int, int, int]

    def __post_init__(self) -> None:
        ids = np.asarray(self.cell_ids, dtype=int)
        weights = np.asarray(self.weights, dtype=float)
        ncell = int(np.prod(self.grid_shape))
        if ids.ndim != 2 or weights.shape != ids.shape or ids.shape[1] == 0:
            raise ValueError("AtomGridMap arrays must be matching non-empty 2D arrays")
        if any(n <= 0 for n in self.grid_shape) or np.any(ids < 0) or np.any(ids >= ncell):
            raise ValueError("AtomGridMap contains invalid cell IDs or grid dimensions")
        if not np.all(np.isfinite(weights)) or np.any(weights < 0.0):
            raise ValueError("AtomGridMap weights must be finite and non-negative")
        if not np.allclose(weights.sum(axis=1), 1.0, rtol=0.0, atol=1.0e-12):
            raise ValueError("Each atom's grid interpolation weights must sum to one")

    @classmethod
    def nearest_cell(cls, cell_ids: np.ndarray, grid_shape: tuple[int, int, int]) -> "AtomGridMap":
        ids = np.asarray(cell_ids, dtype=int).reshape(-1, 1)
        return cls(ids, np.ones_like(ids, dtype=float), grid_shape)


def cloud_in_cell_atom_map(
    xyz_A: np.ndarray,
    grid_shape: tuple[int, int, int],
    spacing_A: tuple[float, float, float],
) -> AtomGridMap:
    """Build trilinear atom-to-cell weights with periodic x/y and open z."""
    xyz = np.asarray(xyz_A, dtype=float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or not np.all(np.isfinite(xyz)):
        raise ValueError("xyz_A must be a finite array with shape (N,3)")
    n = np.asarray(grid_shape, dtype=int)
    spacing = np.asarray(spacing_A, dtype=float)
    if n.shape != (3,) or spacing.shape != (3,) or np.any(n <= 0) or np.any(spacing <= 0.0):
        raise ValueError("grid_shape and spacing_A must each contain three positive values")
    u = xyz / spacing - 0.5
    lower = np.floor(u).astype(int)
    frac = u - lower
    ids: list[np.ndarray] = []
    weights: list[np.ndarray] = []
    for bx, by, bz in ((x, y, z) for x in (0, 1) for y in (0, 1) for z in (0, 1)):
        ix = (lower[:, 0] + bx) % n[0]
        iy = (lower[:, 1] + by) % n[1]
        iz = np.clip(lower[:, 2] + bz, 0, n[2] - 1)
        wx = frac[:, 0] if bx else 1.0 - frac[:, 0]
        wy = frac[:, 1] if by else 1.0 - frac[:, 1]
        wz = frac[:, 2] if bz else 1.0 - frac[:, 2]
        ids.append(np.ravel_multi_index((ix, iy, iz), tuple(n)))
        weights.append(wx * wy * wz)
    return AtomGridMap(np.column_stack(ids), np.column_stack(weights), tuple(int(v) for v in n))


@dataclass(frozen=True)
class AgGridFields:
    fraction_ag: np.ndarray
    source_power_eV_ps: np.ndarray
    atom_map: AtomGridMap
    raw_overfilled_cells: int
    raw_max_fill_fraction: float
    redistributed_ag_volume_A3: float


def build_ag_grid_fields(
    xyz_A: np.ndarray,
    q_eV_ps_A3: np.ndarray,
    grid_shape: tuple[int, int, int],
    spacing_A: tuple[float, float, float],
    ag_atomic_volume_A3: float,
) -> AgGridFields:
    """Build bounded Ag fill and source grids using a conservative CIC map."""
    xyz = np.asarray(xyz_A, dtype=float)
    q = np.asarray(q_eV_ps_A3, dtype=float)
    if xyz.shape != (len(q), 3) or not np.all(np.isfinite(q)) or np.any(q < 0.0):
        raise ValueError("q_eV_ps_A3 must be finite, non-negative, and match xyz_A")
    if not np.isfinite(ag_atomic_volume_A3) or ag_atomic_volume_A3 <= 0.0:
        raise ValueError("Ag atomic volume must be finite and positive")
    mapping = cloud_in_cell_atom_map(xyz, grid_shape, spacing_A)
    ncell = int(np.prod(grid_shape))
    flat_ids = mapping.cell_ids.ravel()
    flat_weights = mapping.weights.ravel()
    atom_volumes = np.full(len(xyz), ag_atomic_volume_A3)
    raw_volume = np.bincount(
        flat_ids,
        weights=np.repeat(atom_volumes, mapping.cell_ids.shape[1]) * flat_weights,
        minlength=ncell,
    ).reshape(grid_shape)
    raw_power = np.bincount(
        flat_ids,
        weights=np.repeat(atom_volumes * q, mapping.cell_ids.shape[1]) * flat_weights,
        minlength=ncell,
    ).reshape(grid_shape)
    cell_volume = float(np.prod(spacing_A))
    raw_fill = raw_volume / cell_volume
    raw_overfilled = int(np.count_nonzero(raw_fill > 1.0))
    raw_max_fill = float(raw_fill.max())

    volume = raw_volume.ravel().copy()
    source = raw_power.ravel().copy()
    shape = tuple(int(v) for v in grid_shape)
    nx, ny, nz = shape
    spacing = np.asarray(spacing_A, dtype=float)
    donors: list[tuple[int, float, float]] = []
    for key in np.flatnonzero(volume > cell_volume):
        excess = volume[key] - cell_volume
        power_per_volume = source[key] / volume[key]
        donors.append((int(key), float(excess), float(power_per_volume)))
        volume[key] = cell_volume
        source[key] = cell_volume * power_per_volume

    tolerance = max(cell_volume * 1.0e-12, 1.0e-12)
    moved = 0.0
    offsets_by_radius: dict[int, list[tuple[int, int, int]]] = {}
    for donor, excess, power_per_volume in sorted(donors, key=lambda item: -item[1]):
        ix = donor // (ny * nz)
        iy = (donor // nz) % ny
        iz = donor % nz
        visited: set[int] = set()
        for radius in range(1, max(shape) // 2 + 1):
            if radius not in offsets_by_radius:
                offsets = list(itertools.product(range(-radius, radius + 1), repeat=3))
                offsets.sort(
                    key=lambda d: (d[0] * spacing[0]) ** 2
                    + (d[1] * spacing[1]) ** 2
                    + (d[2] * spacing[2]) ** 2
                )
                offsets_by_radius[radius] = offsets
            for di, dj, dk in offsets_by_radius[radius]:
                if max(abs(di), abs(dj), abs(dk)) != radius:
                    continue
                z = iz + dk
                if z < 0 or z >= nz:
                    continue
                key = int(np.ravel_multi_index(((ix + di) % nx, (iy + dj) % ny, z), shape))
                if key in visited or key == donor:
                    continue
                visited.add(key)
                if volume[key] <= 0.0 or volume[key] >= cell_volume - tolerance:
                    continue
                amount = min(excess, cell_volume - volume[key])
                volume[key] += amount
                source[key] += amount * power_per_volume
                excess -= amount
                moved += amount
                if excess <= tolerance:
                    break
            if excess <= tolerance:
                break
        if excess > tolerance:
            raise RuntimeError(f"Could not locally redistribute {excess:.6g} A^3 from cell {donor}")

    volume = volume.reshape(shape)
    source = source.reshape(shape)
    fraction = volume / cell_volume
    expected_volume = len(xyz) * ag_atomic_volume_A3
    if abs(float(volume.sum()) - expected_volume) > max(expected_volume * 1.0e-12, 1.0e-8):
        raise RuntimeError("CIC Ag volume mapping is not conservative")
    if float(fraction.max()) > 1.0 + 1.0e-12 or np.any(fraction < 0.0):
        raise RuntimeError("CIC Ag volume mapping produced an invalid fill fraction")
    expected_power = float(np.dot(atom_volumes, q))
    if abs(float(source.sum()) - expected_power) > max(abs(expected_power) * 1.0e-12, 1.0e-10):
        raise RuntimeError("CIC Ag source mapping is not conservative")
    return AgGridFields(
        fraction_ag=fraction,
        source_power_eV_ps=source,
        atom_map=mapping,
        raw_overfilled_cells=raw_overfilled,
        raw_max_fill_fraction=raw_max_fill,
        redistributed_ag_volume_A3=moved,
    )


class AgTTMGrid:
    """Finite-volume electron grid containing Ag fraction only.

    ``fraction_ag`` is a bounded cell occupancy. ``source_power_eV_ps`` is
    already integrated over each cell, not a volumetric density. The x and y
    axes are periodic; both z faces have zero normal electron heat flux. Cells
    with zero Ag occupancy have no electron energy, source, or conduction.
    """

    def __init__(
        self,
        parameters: AgParameters,
        fraction_ag: np.ndarray,
        spacing_A: tuple[float, float, float],
        source_power_eV_ps: np.ndarray | None = None,
        initial_temperature_K: float = 300.0,
        max_substep_deltaT_K: float = 50.0,
        stability_safety: float = 0.4,
    ) -> None:
        self.parameters = parameters
        self.f_ag = np.asarray(fraction_ag, dtype=float).copy()
        if self.f_ag.ndim != 3 or any(n <= 0 for n in self.f_ag.shape):
            raise ValueError("fraction_ag must be a non-empty 3D array")
        if not np.all(np.isfinite(self.f_ag)) or np.any(self.f_ag < 0.0) or np.any(self.f_ag > 1.0 + 1e-12):
            raise ValueError("Ag cell fractions must be finite and in [0, 1]")
        if len(spacing_A) != 3 or any((not np.isfinite(d) or d <= 0.0) for d in spacing_A):
            raise ValueError("All three grid spacings must be finite and positive")
        self.spacing_A = np.asarray(spacing_A, dtype=float)
        self.cell_volume_A3 = float(np.prod(self.spacing_A))
        if source_power_eV_ps is None:
            self.source_power = np.zeros_like(self.f_ag)
        else:
            self.source_power = np.asarray(source_power_eV_ps, dtype=float).copy()
            if self.source_power.shape != self.f_ag.shape:
                raise ValueError("source_power_eV_ps shape must match fraction_ag")
            if not np.all(np.isfinite(self.source_power)) or np.any(self.source_power < 0.0):
                raise ValueError("Source powers must be finite and non-negative")
        if np.any(self.source_power[self.f_ag == 0.0] != 0.0):
            raise ValueError("A cell without Ag cannot receive the Ag-selective source")
        if not (0.0 < stability_safety < 0.5):
            raise ValueError("stability_safety must be between zero and 0.5")
        if max_substep_deltaT_K <= 0.0:
            raise ValueError("max_substep_deltaT_K must be positive")
        self.max_substep_deltaT_K = float(max_substep_deltaT_K)
        self.stability_safety = float(stability_safety)
        self.temperature_K = np.full(self.f_ag.shape, parameters.Tmin, dtype=float)
        if np.any((self.f_ag > 0.0)):
            if initial_temperature_K < parameters.Tmin or initial_temperature_K > parameters.Tmax:
                raise TTMRangeError("Initial electronic temperature is outside the parameter table")
            self.temperature_K[self.f_ag > 0.0] = initial_temperature_K
        self.time_ps = 0.0

    @property
    def electron_energy_eV(self) -> np.ndarray:
        energy = self.f_ag * self.cell_volume_A3 * self.parameters.ue(self.temperature_K)
        energy[self.f_ag == 0.0] = 0.0
        return energy

    def _harmonic_face_conductance(self, te: np.ndarray, axis: int) -> tuple[np.ndarray, bool]:
        """Return face conductance (eV/(ps K)) and whether axis is periodic."""
        k_eff = self.f_ag * self.parameters.ke(te)
        if axis < 2:
            other = np.roll(k_eff, -1, axis=axis)
            denom = k_eff + other
            k_face = np.divide(2.0 * k_eff * other, denom, out=np.zeros_like(k_eff), where=denom > 0.0)
            face_area = self.cell_volume_A3 / self.spacing_A[axis]
            return k_face * face_area / self.spacing_A[axis], True

        low = [slice(None)] * 3
        high = [slice(None)] * 3
        low[axis] = slice(0, -1)
        high[axis] = slice(1, None)
        kl, kh = k_eff[tuple(low)], k_eff[tuple(high)]
        denom = kl + kh
        k_face = np.divide(2.0 * kl * kh, denom, out=np.zeros_like(kl), where=denom > 0.0)
        face_area = self.cell_volume_A3 / self.spacing_A[axis]
        return k_face * face_area / self.spacing_A[axis], False

    def _diffusion_rate(self, te: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        rate = np.zeros_like(te)
        conductance_sum = np.zeros_like(te)
        for axis in range(3):
            if te.shape[axis] == 1:
                continue
            conductance, periodic = self._harmonic_face_conductance(te, axis)
            if periodic:
                other_t = np.roll(te, -1, axis=axis)
                flow = conductance * (other_t - te)
                rate += flow - np.roll(flow, 1, axis=axis)
                conductance_sum += conductance + np.roll(conductance, 1, axis=axis)
            else:
                low = [slice(None)] * 3
                high = [slice(None)] * 3
                low[axis] = slice(0, -1)
                high[axis] = slice(1, None)
                low_t, high_t = tuple(low), tuple(high)
                flow = conductance * (te[high_t] - te[low_t])
                rate[low_t] += flow
                rate[high_t] -= flow
                conductance_sum[low_t] += conductance
                conductance_sum[high_t] += conductance
        return rate, conductance_sum

    def advance(
        self,
        dt_ps: float,
        lattice_temperature_K: np.ndarray | float = 300.0,
        heat_on: bool = True,
        coupling_on: bool = True,
    ) -> StepAccounting:
        """Advance TTM fields, adapting explicit substeps and returning energy.

        Positive ``electron_to_lattice_eV`` is the energy the caller must add
        to Newton-layer Ag kinetic energy. Negative values remove the same
        amount from Ag. Heat-source power is zero when ``heat_on`` is false.
        """
        if not np.isfinite(dt_ps) or dt_ps <= 0.0:
            raise ValueError("dt_ps must be finite and positive")
        tl = np.broadcast_to(np.asarray(lattice_temperature_K, dtype=float), self.f_ag.shape)
        active_cells = self.f_ag > 0.0
        if not np.all(np.isfinite(tl[active_cells])) or np.any(tl[active_cells] <= 0.0):
            raise ValueError("Lattice temperatures must be finite and positive")
        tl = np.where(active_cells, tl, self.parameters.Tmin)

        initial_total = float(self.electron_energy_eV.sum())
        source_total = 0.0
        lattice_total = 0.0
        diffusion_total = 0.0
        lattice_energy_grid = np.zeros_like(self.f_ag)
        elapsed = 0.0
        substeps = 0
        active = self.f_ag > 0.0
        max_ue = self.parameters.ue(self.parameters.Tmax)
        min_ue = self.parameters.ue(self.parameters.Tmin)

        while elapsed < dt_ps - max(1.0e-15, dt_ps * 1.0e-13):
            te = self.temperature_K
            ce = self.parameters.ce(te)
            g = self.parameters.g(te)
            current_energy = self.electron_energy_eV
            diffusion_rate, gsum = self._diffusion_rate(te)
            heat_capacity_cell = self.f_ag * self.cell_volume_A3 * ce
            limits: list[float] = [dt_ps - elapsed]
            diff_active = active & (gsum > 0.0)
            if np.any(diff_active):
                limits.append(
                    self.stability_safety
                    * float(np.min(heat_capacity_cell[diff_active] / gsum[diff_active]))
                )
            coupling_active = active & (g > 0.0)
            if np.any(coupling_active):
                limits.append(
                    self.stability_safety
                    * float(np.min(ce[coupling_active] / g[coupling_active]))
                )
            if heat_on:
                source_rate = self.source_power
                source_active = active & (source_rate > 0.0)
                if np.any(source_active):
                    heating_rate_K_ps = source_rate[source_active] / heat_capacity_cell[source_active]
                    limits.append(self.max_substep_deltaT_K / float(np.max(heating_rate_K_ps)))
            else:
                source_rate = np.zeros_like(self.source_power)
            subdt = min(limits)
            if not np.isfinite(subdt) or subdt <= 0.0:
                raise RuntimeError("TTM substep selection produced a non-positive step")
            if coupling_on:
                lattice_rate = self.f_ag * self.cell_volume_A3 * g * (te - tl)
            else:
                lattice_rate = np.zeros_like(te)
            emax = self.f_ag * self.cell_volume_A3 * max_ue
            emin = self.f_ag * self.cell_volume_A3 * min_ue
            while True:
                lattice_energy = lattice_rate * subdt
                new_energy = current_energy + subdt * (diffusion_rate + source_rate) - lattice_energy
                new_energy[~active] = 0.0
                bad_low = active & (new_energy < emin - np.maximum(1.0e-13, np.abs(emin) * 1.0e-12))
                bad_high = active & (new_energy > emax + np.maximum(1.0e-13, emax * 1.0e-12))
                if not np.any(bad_low) and not np.any(bad_high):
                    break
                if subdt <= max(1.0e-14, dt_ps * 1.0e-12):
                    reason = "below 300 K table floor" if np.any(bad_low) else "above 20000 K table ceiling"
                    raise TTMRangeError(f"TTM energy update leaves the tabulated range ({reason})")
                subdt *= 0.5

            low_roundoff = active & (new_energy < emin)
            new_energy[low_roundoff] = emin[low_roundoff]
            energy_density = np.zeros_like(new_energy)
            energy_density[active] = new_energy[active] / (
                self.f_ag[active] * self.cell_volume_A3
            )
            new_te = self.temperature_K.copy()
            new_te[active] = self.parameters.temperature_from_ue(energy_density[active])

            # Recompute actual accepted increments after an optional range retry.
            accepted_dt = subdt
            source_energy = float((source_rate * accepted_dt).sum())
            coupling_energy = float(lattice_energy.sum())
            diffusion_energy = float((diffusion_rate * accepted_dt).sum())
            self.temperature_K = new_te
            self.time_ps += accepted_dt
            elapsed += accepted_dt
            source_total += source_energy
            lattice_total += coupling_energy
            diffusion_total += diffusion_energy
            lattice_energy_grid += lattice_energy
            substeps += 1

        final_total = float(self.electron_energy_eV.sum())
        residual = (final_total - initial_total) - source_total + lattice_total - diffusion_total
        return StepAccounting(
            elapsed_ps=elapsed,
            source_to_electron_eV=source_total,
            electron_to_lattice_eV=lattice_total,
            diffusion_net_eV=diffusion_total,
            conservation_residual_eV=residual,
            substeps=substeps,
            lattice_energy_per_cell_eV=lattice_energy_grid,
        )


def rescale_lattice_kinetic(
    velocities_A_ps: np.ndarray,
    masses_g_mol: np.ndarray,
    atom_grid_map: AtomGridMap,
    energy_to_lattice_eV: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply equal/opposite cell energy by deterministic Ag velocity scaling.

    This preserves exact cell kinetic-energy changes but not local momentum and
    does not add Langevin fluctuations. It is a mean-field transfer operator.
    Callers must pass only the Newton-layer Ag atoms represented in the TTM grid.
    """
    velocity = np.asarray(velocities_A_ps, dtype=float).copy()
    mass = np.asarray(masses_g_mol, dtype=float)
    mapping = atom_grid_map
    target = np.asarray(energy_to_lattice_eV, dtype=float)
    if velocity.ndim != 2 or velocity.shape[1] != 3:
        raise ValueError("velocities must have shape (N,3)")
    if mass.shape != (len(velocity),) or mapping.cell_ids.shape[0] != len(velocity):
        raise ValueError("mass and atom-grid map must match atom count")
    if not np.all(np.isfinite(velocity)) or np.any(mass <= 0.0):
        raise ValueError("Velocities must be finite and masses positive")
    if tuple(target.shape) != mapping.grid_shape:
        raise ValueError("Energy-transfer grid shape must match AtomGridMap")
    if not np.all(np.isfinite(target)):
        raise ValueError("Cell energy transfers must be finite")

    flat_cells = mapping.cell_ids.ravel()
    flat_weights = mapping.weights.ravel()
    cell_weight_sum = np.bincount(
        flat_cells, weights=flat_weights, minlength=target.size
    )
    normalized_reverse = np.divide(
        mapping.weights,
        cell_weight_sum[mapping.cell_ids],
        out=np.zeros_like(mapping.weights),
        where=cell_weight_sum[mapping.cell_ids] > 0.0,
    )
    atom_transfer = np.sum(
        normalized_reverse * target.ravel()[mapping.cell_ids], axis=1
    )
    kinetic = 0.5 * M_VV2E_METAL * mass * np.einsum("ij,ij->i", velocity, velocity)
    new_kinetic = kinetic + atom_transfer
    active_transfer = atom_transfer != 0.0
    negative = np.flatnonzero(active_transfer & (new_kinetic < -1.0e-12))
    if negative.size:
        raise ValueError(
            f"Atom {int(negative[0])} would have negative lattice kinetic energy"
        )
    zero_kinetic = np.flatnonzero(active_transfer & (kinetic <= 0.0))
    if zero_kinetic.size:
        raise ValueError(
            f"Atom {int(zero_kinetic[0])} has zero kinetic energy; cannot rescale"
        )
    scale_squared = np.ones_like(kinetic)
    np.divide(
        np.maximum(new_kinetic, 0.0),
        kinetic,
        out=scale_squared,
        where=active_transfer,
    )
    velocity *= np.sqrt(scale_squared)[:, None]

    before_atom = 0.5 * M_VV2E_METAL * mass * np.sum(
        np.asarray(velocities_A_ps, dtype=float) ** 2, axis=1
    )
    after_atom = 0.5 * M_VV2E_METAL * mass * np.sum(velocity**2, axis=1)
    actual_total = float(np.sum(after_atom - before_atom))
    target_total = float(np.sum(target))
    if not np.isclose(actual_total, target_total, rtol=1.0e-10, atol=1.0e-12):
        raise RuntimeError("Applied lattice kinetic-energy transfer does not match target")
    # Per-cell reports use the same conservative adjoint weights as the map.
    actual_flat = np.bincount(
        mapping.cell_ids.ravel(),
        weights=np.repeat(after_atom - before_atom, mapping.cell_ids.shape[1]) * mapping.weights.ravel(),
        minlength=target.size,
    )
    return velocity, actual_flat.reshape(target.shape)


def lattice_temperature_by_cell(
    velocities_A_ps: np.ndarray,
    masses_g_mol: np.ndarray,
    atom_grid_map: AtomGridMap,
    smoothing_radius_cells: int = 1,
) -> np.ndarray:
    """Coarse-grained kinetic temperature of Newton-layer Ag atoms per cell.

    A small box filter is applied to both kinetic energy and represented atom
    count. This reduces noisy temperatures where a thin Ag filament occupies
    only a small part of a TTM cell. x/y use periodic smoothing; z uses an
    edge-clamped filter, consistent with the zero-flux free-surface treatment.
    """
    velocity = np.asarray(velocities_A_ps, dtype=float)
    mass = np.asarray(masses_g_mol, dtype=float)
    mapping = atom_grid_map
    cell_count = int(np.prod(mapping.grid_shape))
    if velocity.shape != (len(mass), 3) or mapping.cell_ids.shape[0] != len(mass):
        raise ValueError("Velocity, mass, and atom-grid map have inconsistent shapes")
    if not np.all(np.isfinite(velocity)) or np.any(mass <= 0.0):
        raise ValueError("Velocities must be finite and masses positive")
    kinetic = 0.5 * M_VV2E_METAL * mass * np.sum(velocity**2, axis=1)
    flat_ids = mapping.cell_ids.ravel()
    flat_weights = mapping.weights.ravel()
    kinetic_sum = np.bincount(
        flat_ids,
        weights=np.repeat(kinetic, mapping.cell_ids.shape[1]) * flat_weights,
        minlength=cell_count,
    )
    atom_count = np.bincount(flat_ids, weights=flat_weights, minlength=cell_count)
    kinetic_grid = kinetic_sum.reshape(mapping.grid_shape)
    count_grid = atom_count.reshape(mapping.grid_shape)

    if smoothing_radius_cells < 0:
        raise ValueError("smoothing_radius_cells cannot be negative")
    for _ in range(smoothing_radius_cells):
        for axis in (0, 1):
            kinetic_grid = (
                np.roll(kinetic_grid, 1, axis=axis)
                + kinetic_grid
                + np.roll(kinetic_grid, -1, axis=axis)
            ) / 3.0
            count_grid = (
                np.roll(count_grid, 1, axis=axis)
                + count_grid
                + np.roll(count_grid, -1, axis=axis)
            ) / 3.0
        for array_name in ("kinetic", "count"):
            array = kinetic_grid if array_name == "kinetic" else count_grid
            lower_indices = np.maximum(np.arange(array.shape[2]) - 1, 0)
            upper_indices = np.minimum(np.arange(array.shape[2]) + 1, array.shape[2] - 1)
            smoothed = (
                np.take(array, lower_indices, axis=2)
                + array
                + np.take(array, upper_indices, axis=2)
            ) / 3.0
            if array_name == "kinetic":
                kinetic_grid = smoothed
            else:
                count_grid = smoothed

    temperature = np.full(cell_count, np.nan, dtype=float)
    active = count_grid.ravel() > 0.0
    temperature[active] = 2.0 * kinetic_grid.ravel()[active] / (
        3.0 * count_grid.ravel()[active] * KB_EV_K
    )
    return temperature.reshape(mapping.grid_shape)
