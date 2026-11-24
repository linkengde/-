"""Ag–W extension of the repository's Python finite-volume TTM reference core.

This module reuses the atom-grid mapping, electron finite-volume operator,
energy table handling contract, and lattice kinetic-energy handoff from
``models/ag_ti3sic2_qsgs_1to1_finer/ag_ttm_reference.py``.  The legacy Ag–Ti3SiC2
files remain untouched.  Ag and W have separate electron grids; direct laser
source is applied to Ag only, while a provisional face-based electronic
interface conductance exchanges energy between the two grids.

This remains an exploratory Python reference solver, not a compiled LAMMPS
TTM fix.  Interface exchange is resolved only across unlike-phase cell faces;
sub-cell interface area is not inferred.
"""
from __future__ import annotations

import csv
import importlib.util
from dataclasses import dataclass
from pathlib import Path
import sys

import numpy as np


def _load_legacy_core():
    for parent in Path(__file__).resolve().parents:
        source = parent / "models/ag_ti3sic2_qsgs_1to1_finer/ag_ttm_reference.py"
        if source.is_file():
            spec = importlib.util.spec_from_file_location("_ag_ttm_reference_core", source)
            if spec is None or spec.loader is None:
                raise ImportError(f"Cannot load existing TTM core at {source}")
            module = importlib.util.module_from_spec(spec)
            sys.modules[spec.name] = module
            spec.loader.exec_module(module)
            return module
    raise ImportError("The repository's existing Ag TTM reference core was not found")


_CORE = _load_legacy_core()
AgTTMGrid = _CORE.AgTTMGrid
AtomGridMap = _CORE.AtomGridMap
M_VV2E_METAL = _CORE.M_VV2E_METAL
build_ag_grid_fields = _CORE.build_ag_grid_fields
cloud_in_cell_atom_map = _CORE.cloud_in_cell_atom_map
lattice_temperature_by_cell = _CORE.lattice_temperature_by_cell
rescale_lattice_kinetic = _CORE.rescale_lattice_kinetic

PACK = next(
    parent for parent in Path(__file__).resolve().parents
    if (parent / "ttm_data/AgW_TTM_baseline_v1").is_dir()
)
DATA = PACK / "ttm_data"
BASELINE = DATA / "AgW_TTM_baseline_v1"

J_M3_TO_EV_A3 = 6.241509074e-12
W_M3_TO_EV_PS_A3 = 6.241509074e-24
W_M_TO_EV_PS_A = 6.241509074e-4
W_M2_TO_EV_PS_A2 = 6.241509074e-14
KB_EV_K = 8.617333262145e-5
AG_ATOMIC_VOLUME_A3 = 4.09**3 / 4.0
W_ATOMIC_VOLUME_A3 = 3.16**3 / 2.0


class AgWParameterError(ValueError):
    """Input property or temperature is outside the selected baseline."""


class BaselineMetalParameters:
    """Duck-typed property interface consumed by the existing AgTTMGrid core."""

    def __init__(self, element: str):
        if element not in {"Ag", "W"}:
            raise ValueError("element must be Ag or W")
        self.element = element
        if element == "Ag":
            path = BASELINE / "Ag_TTM_Lin2008_digitized_300_20000K.csv"
            ce_col = "Ce_Lin2008_digitized_J_m-3_K-1"
            g_col = "G_Lin2008_digitized_W_m-3_K-1"
        else:
            path = DATA / "W_TTM_Lin2008_digitized_300_20000K.csv"
            ce_col = "Ce_Lin2008_J_m-3_K-1"
            g_col = "G_Lin2008_digitized_W_m-3_K-1"
        temperatures, ce_si, g_si = [], [], []
        with path.open("r", newline="", encoding="utf-8") as stream:
            for row in csv.DictReader(stream):
                temperatures.append(float(row["Te_K"]))
                ce_si.append(float(row[ce_col]))
                g_si.append(float(row[g_col]))
        self.temperature_K = np.asarray(temperatures)
        if (len(temperatures) < 2 or np.any(np.diff(self.temperature_K) <= 0.0)
                or not np.all(np.isfinite(ce_si)) or not np.all(np.isfinite(g_si))
                or np.any(np.asarray(ce_si) <= 0.0) or np.any(np.asarray(g_si) <= 0.0)):
            raise ValueError(f"Invalid Lin 2008 {element} Ce/G table")
        if self.temperature_K[0] != 300.0:
            raise ValueError("Baseline parameter tables must start at 300 K")
        self.ce_eV_A3_K = np.asarray(ce_si) * J_M3_TO_EV_A3
        self.g_eV_ps_A3_K = np.asarray(g_si) * W_M3_TO_EV_PS_A3
        self.Tmin = 1.0  # explicit, documented low-Te extension below 300 K
        self.Tmax = float(self.temperature_K[-1])
        self._ce_gamma = float(self.ce_eV_A3_K[0] / 300.0)
        self.ue_eV_A3 = np.zeros_like(self.ce_eV_A3_K)
        self.ue_eV_A3[1:] = np.cumsum(
            0.5 * (self.ce_eV_A3_K[:-1] + self.ce_eV_A3_K[1:])
            * np.diff(self.temperature_K)
        )
        self.lattice_temperature_K: np.ndarray | float = 300.0
        if element == "Ag":
            import json
            model = json.loads((BASELINE / "baseline_parameters.json").read_text(encoding="utf-8"))
            self.a_over_b_K_inv = float(model["material_parameters"]["Ag"]["ke_Te_Tl"]["A_over_B_K_minus1"])
            self.ke0_W_m_K = float(model["material_parameters"]["Ag"]["ke_Te_Tl"]["ke0_W_m_K"])
        else:
            import json
            model = json.loads((BASELINE / "baseline_parameters.json").read_text(encoding="utf-8"))
            self.ke0_W_m_K = float(model["material_parameters"]["W"]["ke_Te_Tl"]["ke0_W_m_K"])
            self.ke_exponent = float(model["material_parameters"]["W"]["ke_Te_Tl"]["n"])

    def set_lattice_temperature(self, temperature_K: np.ndarray | float) -> None:
        value = np.asarray(temperature_K, dtype=float)
        if not np.all(np.isfinite(value)) or np.any(value <= 0.0):
            raise AgWParameterError("Lattice temperature must be finite and positive")
        self.lattice_temperature_K = value

    def _temperature(self, temperature_K: np.ndarray | float) -> np.ndarray:
        te = np.asarray(temperature_K, dtype=float)
        if not np.all(np.isfinite(te)) or np.any(te < self.Tmin) or np.any(te > self.Tmax):
            lo = float(np.min(te)) if te.size else float(te)
            hi = float(np.max(te)) if te.size else float(te)
            raise AgWParameterError(
                f"{self.element} Te range [{lo:g},{hi:g}] K outside [1,{self.Tmax:g}] K"
            )
        return te

    def _table(self, values: np.ndarray, te: np.ndarray) -> np.ndarray:
        return np.interp(np.maximum(te, 300.0), self.temperature_K, values)

    def ce(self, temperature_K: np.ndarray | float) -> np.ndarray:
        te = self._temperature(temperature_K)
        value = self._table(self.ce_eV_A3_K, te)
        return np.where(te < 300.0, self._ce_gamma * te, value)

    def g(self, temperature_K: np.ndarray | float) -> np.ndarray:
        te = self._temperature(temperature_K)
        value = self._table(self.g_eV_ps_A3_K, te)
        # Explicit convenience continuation below the table, inherited from
        # the prior reference solver: constant G below 300 K.
        return np.where(te < 300.0, self.g_eV_ps_A3_K[0], value)

    def ke(self, temperature_K: np.ndarray | float) -> np.ndarray:
        te = self._temperature(temperature_K)
        tl = np.broadcast_to(np.asarray(self.lattice_temperature_K, dtype=float), te.shape)
        if not np.all(np.isfinite(tl)) or np.any(tl <= 0.0):
            raise AgWParameterError("Invalid lattice temperature in conductivity closure")
        if self.element == "Ag":
            a = self.a_over_b_K_inv
            ke_si = self.ke0_W_m_K * (self.ce(te) / self.ce(300.0)) * (
                (300.0 + a * 300.0**2) / (tl + a * te**2)
            )
        else:
            ke_si = self.ke0_W_m_K * (self.ce(te) / self.ce(300.0)) * (
                300.0 / tl
            ) ** self.ke_exponent
        return ke_si * W_M_TO_EV_PS_A

    def ue(self, temperature_K: np.ndarray | float) -> np.ndarray:
        te = self._temperature(temperature_K)
        value = self._table(self.ue_eV_A3, te)
        low = self.ue_eV_A3[0] + 0.5 * self._ce_gamma * (te**2 - 300.0**2)
        return np.where(te < 300.0, low, value)

    def temperature_from_ue(self, energy_density_eV_A3: np.ndarray | float) -> np.ndarray:
        energy = np.asarray(energy_density_eV_A3, dtype=float)
        if not np.all(np.isfinite(energy)):
            raise AgWParameterError("Non-finite electron energy density")
        low_energy = self.ue(self.Tmin)
        if np.any(energy < low_energy - 1.0e-13) or np.any(energy > self.ue_eV_A3[-1] + 1.0e-13):
            raise AgWParameterError(f"{self.element} electron energy outside the baseline table")
        clipped = np.clip(energy, low_energy, self.ue_eV_A3[-1])
        low_region = clipped < self.ue_eV_A3[0]
        ordinary = np.interp(clipped, self.ue_eV_A3, self.temperature_K)
        low_te = np.sqrt(np.maximum(0.0, 300.0**2 + 2.0 * clipped / self._ce_gamma))
        return np.where(low_region, low_te, ordinary)


@dataclass(frozen=True)
class AgWStepAccounting:
    elapsed_ps: float
    source_to_electron_eV: float
    electron_to_lattice_eV: float
    diffusion_net_eV: float
    interface_net_eV: float
    conservation_residual_eV: float
    substeps: int
    lattice_energy_ag_eV: np.ndarray
    lattice_energy_w_eV: np.ndarray


class AgWTTMGrid:
    """Two phase electron grids using the existing reference TTM operator.

    The phase grids share one finite-volume mesh but store separate `Te` and
    electron energy fields. Ag receives the spatial source; W receives none.
    Electron diffusion is within each phase. Ag/W electron exchange is an
    explicitly parameterized face conductance, equal and opposite in energy.
    """

    def __init__(
        self,
        fraction_ag: np.ndarray,
        fraction_w: np.ndarray,
        spacing_A: tuple[float, float, float],
        source_power_ag_eV_ps: np.ndarray | None = None,
        interface_conductance_W_m2_K: float = 1.0e9,
        initial_temperature_K: float = 300.0,
        interface_safety: float = 0.35,
    ) -> None:
        f_ag = np.asarray(fraction_ag, dtype=float)
        f_w = np.asarray(fraction_w, dtype=float)
        if f_ag.shape != f_w.shape or f_ag.ndim != 3 or not f_ag.size:
            raise ValueError("Ag and W fractions must be matching non-empty 3D arrays")
        if np.any(f_ag < 0.0) or np.any(f_w < 0.0) or np.any(f_ag > 1.0) or np.any(f_w > 1.0):
            raise ValueError("Phase fractions must lie in [0,1]")
        if np.any(f_ag + f_w > 1.0 + 1e-8):
            raise ValueError("Ag+W mapped volume exceeds cell capacity; resolve mapping overfill explicitly")
        if not np.isfinite(interface_conductance_W_m2_K) or interface_conductance_W_m2_K < 0.0:
            raise ValueError("Interface conductance must be finite and non-negative")
        if not (0.0 < interface_safety < 0.5):
            raise ValueError("interface_safety must be between zero and 0.5")
        self.fraction_ag = f_ag.copy()
        self.fraction_w = f_w.copy()
        self.spacing_A = tuple(float(v) for v in spacing_A)
        self.cell_volume_A3 = float(np.prod(self.spacing_A))
        self.ag_parameters = BaselineMetalParameters("Ag")
        self.w_parameters = BaselineMetalParameters("W")
        source = np.zeros_like(f_ag) if source_power_ag_eV_ps is None else np.asarray(source_power_ag_eV_ps, dtype=float)
        if source.shape != f_ag.shape or np.any(~np.isfinite(source)) or np.any(source < 0.0):
            raise ValueError("Ag source power must be a finite non-negative array matching the mesh")
        if np.any(source[f_ag == 0.0] != 0.0):
            raise ValueError("Cells without Ag cannot receive the source")
        self.source_power_ag_base = source.copy()
        self.source_envelope = 1.0
        self.interface_h_metal = interface_conductance_W_m2_K * W_M2_TO_EV_PS_A2
        self.interface_safety = interface_safety
        self.ag = AgTTMGrid(
            self.ag_parameters, f_ag, self.spacing_A, source,
            initial_temperature_K=initial_temperature_K,
        )
        self.w = AgTTMGrid(
            self.w_parameters, f_w, self.spacing_A, np.zeros_like(f_w),
            initial_temperature_K=initial_temperature_K,
        )

    @property
    def electron_energy_ag_eV(self) -> np.ndarray:
        return self.ag.electron_energy_eV

    @property
    def electron_energy_w_eV(self) -> np.ndarray:
        return self.w.electron_energy_eV

    @property
    def temperature_ag_K(self) -> np.ndarray:
        return self.ag.temperature_K

    @property
    def temperature_w_K(self) -> np.ndarray:
        return self.w.temperature_K

    def set_source_envelope(self, value: float) -> None:
        if not np.isfinite(value) or value < 0.0:
            raise ValueError("Source envelope S(t) must be finite and non-negative")
        self.source_envelope = float(value)
        self.ag.source_power = self.source_power_ag_base * self.source_envelope

    def _face_specs(self):
        """Yield each +face once with its two possible unlike-phase contacts."""
        nx, ny, nz = self.fraction_ag.shape
        for axis in range(3):
            if self.fraction_ag.shape[axis] <= 1:
                continue
            area = self.cell_volume_A3 / self.spacing_A[axis]
            if axis < 2:
                low = [slice(None)] * 3
                high = [slice(None)] * 3
                f_ag_i = self.fraction_ag
                f_w_i = self.fraction_w
                f_ag_j = np.roll(self.fraction_ag, -1, axis=axis)
                f_w_j = np.roll(self.fraction_w, -1, axis=axis)
                yield axis, None, None, area, f_ag_i, f_w_j, f_ag_j, f_w_i
            else:
                low = [slice(None)] * 3
                high = [slice(None)] * 3
                low[axis], high[axis] = slice(0, -1), slice(1, None)
                sl, sh = tuple(low), tuple(high)
                yield axis, sl, sh, area, self.fraction_ag[sl], self.fraction_w[sh], self.fraction_ag[sh], self.fraction_w[sl]

    def _interface_conductance_sums(self) -> tuple[np.ndarray, np.ndarray]:
        ag_sum = np.zeros_like(self.fraction_ag)
        w_sum = np.zeros_like(self.fraction_w)
        for axis, sl, sh, area, ag_i, w_j, ag_j, w_i in self._face_specs():
            g1 = self.interface_h_metal * area * ag_i * w_j
            g2 = self.interface_h_metal * area * ag_j * w_i
            if axis < 2:
                ag_sum += g1 + np.roll(g2, 1, axis=axis)
                w_sum += np.roll(g1, 1, axis=axis) + g2
            else:
                ag_sum[sl] += g1
                w_sum[sh] += g1
                ag_sum[sh] += g2
                w_sum[sl] += g2
        return ag_sum, w_sum

    def _apply_interface_exchange(self, dt_ps: float) -> float:
        if self.interface_h_metal == 0.0:
            return 0.0
        d_ag = np.zeros_like(self.fraction_ag)
        d_w = np.zeros_like(self.fraction_w)
        te_ag, te_w = self.ag.temperature_K, self.w.temperature_K
        for axis, sl, sh, area, ag_i, w_j, ag_j, w_i in self._face_specs():
            if axis < 2:
                te_ag_i = te_ag
                te_w_i = te_w
                te_ag_j = np.roll(te_ag, -1, axis=axis)
                te_w_j = np.roll(te_w, -1, axis=axis)
                e1 = self.interface_h_metal * area * ag_i * w_j * (te_ag_i - te_w_j) * dt_ps
                e2 = self.interface_h_metal * area * ag_j * w_i * (te_ag_j - te_w_i) * dt_ps
                d_ag -= e1
                d_w += np.roll(e1, 1, axis=axis)
                d_ag -= np.roll(e2, 1, axis=axis)
                d_w += e2
            else:
                e1 = self.interface_h_metal * area * ag_i * w_j * (te_ag[sl] - te_w[sh]) * dt_ps
                e2 = self.interface_h_metal * area * ag_j * w_i * (te_ag[sh] - te_w[sl]) * dt_ps
                d_ag[sl] -= e1
                d_w[sh] += e1
                d_ag[sh] -= e2
                d_w[sl] += e2

        for phase_grid, fraction, delta, params in (
            (self.ag, self.fraction_ag, d_ag, self.ag_parameters),
            (self.w, self.fraction_w, d_w, self.w_parameters),
        ):
            active = fraction > 0.0
            if not np.any(active):
                continue
            energy = phase_grid.electron_energy_eV + delta
            density = np.zeros_like(energy)
            density[active] = energy[active] / (fraction[active] * self.cell_volume_A3)
            new_t = phase_grid.temperature_K.copy()
            new_t[active] = params.temperature_from_ue(density[active])
            phase_grid.temperature_K = new_t
        return float(d_ag.sum() + d_w.sum())

    def advance(
        self,
        dt_ps: float,
        lattice_temperature_ag_K: np.ndarray,
        lattice_temperature_w_K: np.ndarray,
        heat_on: bool = True,
    ) -> AgWStepAccounting:
        if not np.isfinite(dt_ps) or dt_ps <= 0.0:
            raise ValueError("dt_ps must be finite and positive")
        tl_ag = np.broadcast_to(np.asarray(lattice_temperature_ag_K, dtype=float), self.fraction_ag.shape)
        tl_w = np.broadcast_to(np.asarray(lattice_temperature_w_K, dtype=float), self.fraction_w.shape)
        active_ag, active_w = self.fraction_ag > 0.0, self.fraction_w > 0.0
        if (not np.all(np.isfinite(tl_ag[active_ag])) or np.any(tl_ag[active_ag] <= 0.0)
                or not np.all(np.isfinite(tl_w[active_w])) or np.any(tl_w[active_w] <= 0.0)):
            raise AgWParameterError("Active-cell lattice temperatures must be finite and positive")
        # The reused single-phase operator requires finite full-sized arrays
        # although it ignores lattice temperatures in empty phase cells.
        tl_ag = np.where(active_ag, tl_ag, 300.0)
        tl_w = np.where(active_w, tl_w, 300.0)
        self.ag_parameters.set_lattice_temperature(tl_ag)
        self.w_parameters.set_lattice_temperature(tl_w)
        e0 = float(self.ag.electron_energy_eV.sum() + self.w.electron_energy_eV.sum())
        source_total = lattice_total = diffusion_total = interface_total = 0.0
        lattice_ag = np.zeros_like(self.fraction_ag)
        lattice_w = np.zeros_like(self.fraction_w)
        g_ag, g_w = self._interface_conductance_sums()
        max_steps = 100000
        elapsed = 0.0
        substeps = 0
        while elapsed < dt_ps - max(1.0e-15, dt_ps * 1.0e-13):
            limits = [dt_ps - elapsed]
            for grid, fraction, params, gsum in (
                (self.ag, self.fraction_ag, self.ag_parameters, g_ag),
                (self.w, self.fraction_w, self.w_parameters, g_w),
            ):
                active = (fraction > 0.0) & (gsum > 0.0)
                if np.any(active):
                    cap = fraction * self.cell_volume_A3 * params.ce(grid.temperature_K)
                    limits.append(self.interface_safety * float(np.min(cap[active] / gsum[active])))
            subdt = min(limits)
            if not np.isfinite(subdt) or subdt <= 0.0:
                raise RuntimeError("Ag–W interface substep selection failed")
            acc_ag = self.ag.advance(subdt, tl_ag, heat_on=heat_on)
            acc_w = self.w.advance(subdt, tl_w, heat_on=False)
            iface = self._apply_interface_exchange(subdt)
            source_total += acc_ag.source_to_electron_eV + acc_w.source_to_electron_eV
            lattice_total += acc_ag.electron_to_lattice_eV + acc_w.electron_to_lattice_eV
            diffusion_total += acc_ag.diffusion_net_eV + acc_w.diffusion_net_eV
            interface_total += iface
            lattice_ag += acc_ag.lattice_energy_per_cell_eV
            lattice_w += acc_w.lattice_energy_per_cell_eV
            elapsed += subdt
            substeps += 1
            if substeps > max_steps:
                raise RuntimeError("Ag–W interface exchange exceeded the substep safety limit")
        e1 = float(self.ag.electron_energy_eV.sum() + self.w.electron_energy_eV.sum())
        residual = (e1 - e0) - source_total + lattice_total - diffusion_total
        return AgWStepAccounting(
            elapsed_ps=elapsed,
            source_to_electron_eV=source_total,
            electron_to_lattice_eV=lattice_total,
            diffusion_net_eV=diffusion_total,
            interface_net_eV=interface_total,
            conservation_residual_eV=residual,
            substeps=substeps,
            lattice_energy_ag_eV=lattice_ag,
            lattice_energy_w_eV=lattice_w,
        )


def gaussian_ag_source_power(
    xyz_absolute_A: np.ndarray,
    ag_atomic_volume_A3: float = AG_ATOMIC_VOLUME_A3,
    i0_eV_ps_A2: float = 5.0,
    radius_A: float = 30.0,
    skin_depth_A: float = 20.0,
    center_xy_A: tuple[float, float] = (100.0, 100.0),
    surface_z_A: float = 200.0,
    allow_above_surface: bool = False,
) -> np.ndarray:
    """Per-atom Ag source quadrature at S(t)=1, in eV/ps/A^3.

    Above-surface atoms are rejected by default. The gated short pilot can opt
    into a fixed nominal z=surface plane, in which case their source depth is
    clamped to zero; this does not model a moving free surface.
    """
    xyz = np.asarray(xyz_absolute_A, dtype=float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or not np.all(np.isfinite(xyz)):
        raise ValueError("xyz_absolute_A must be a finite (N,3) array")
    if min(i0_eV_ps_A2, radius_A, skin_depth_A, ag_atomic_volume_A3) <= 0.0:
        raise ValueError("Source and atomic-volume parameters must be positive")
    if not allow_above_surface and np.any(xyz[:, 2] > surface_z_A + 1e-9):
        raise ValueError("Atom is above the source surface; moving-surface handling is undefined")
    dxy = xyz[:, :2] - np.asarray(center_xy_A, dtype=float)
    dxy = (dxy + 100.0) % 200.0 - 100.0
    depth = np.maximum(0.0, surface_z_A - xyz[:, 2])
    return (i0_eV_ps_A2 / skin_depth_A) * np.exp(
        -np.sum(dxy**2, axis=1) / radius_A**2 - depth / skin_depth_A
    )


def map_agw_newton_layer(
    xyz_absolute_A: np.ndarray,
    types: np.ndarray,
    grid_shape: tuple[int, int, int],
    spacing_A: tuple[float, float, float],
    source_kwargs: dict | None = None,
    newton_mask: np.ndarray | None = None,
) -> dict:
    """Map Newton-layer atom volumes/source to separate Ag/W phase grids.

    Coordinates are assumed in the 200×200×300 Å geometry with Newton atoms
    restricted to 40≤z≤200 Å. An explicit mask can preserve the original
    LAMMPS Newton-group membership through a short equilibration. The two phase
    fills are checked together; volume overfill is never silently clipped.
    """
    xyz = np.asarray(xyz_absolute_A, dtype=float)
    typ = np.asarray(types, dtype=int)
    if xyz.shape != (len(typ), 3):
        raise ValueError("types and coordinates have inconsistent shapes")
    if newton_mask is None:
        newton = (xyz[:, 2] >= 40.0 - 1e-9) & (xyz[:, 2] <= 200.0 + 1e-9)
    else:
        newton = np.asarray(newton_mask, dtype=bool)
        if newton.shape != (len(xyz),):
            raise ValueError("newton_mask must match the atom count")
    idx_ag = np.flatnonzero(newton & (typ == 2))
    idx_w = np.flatnonzero(newton & (typ == 1))
    if not len(idx_ag) or not len(idx_w):
        raise ValueError("The mapped Newton layer must contain both Ag and W atoms")
    origin = np.array([0.0, 0.0, 40.0])
    xyz_ag = xyz[idx_ag] - origin
    xyz_w = xyz[idx_w] - origin
    kwargs = dict(source_kwargs or {})
    q_ag = gaussian_ag_source_power(xyz[idx_ag], **kwargs)
    ag_fields = build_ag_grid_fields(xyz_ag, q_ag, grid_shape, spacing_A, AG_ATOMIC_VOLUME_A3)
    w_fields = build_ag_grid_fields(xyz_w, np.zeros(len(idx_w)), grid_shape, spacing_A, W_ATOMIC_VOLUME_A3)
    total_fill = ag_fields.fraction_ag + w_fields.fraction_ag
    if float(total_fill.max()) > 1.0 + 1e-8:
        raise ValueError(
            f"Mapped Ag+W fill peaks at {float(total_fill.max()):.6g}; refine the grid or repair geometry before TTM"
        )
    return {
        "fraction_ag": ag_fields.fraction_ag,
        "fraction_w": w_fields.fraction_ag,
        "source_power_ag_eV_ps": ag_fields.source_power_eV_ps,
        "atom_map_ag": ag_fields.atom_map,
        "atom_map_w": w_fields.atom_map,
        "atom_indices_ag": idx_ag,
        "atom_indices_w": idx_w,
        "mapped_ag_volume_A3": float(ag_fields.fraction_ag.sum() * np.prod(spacing_A)),
        "mapped_w_volume_A3": float(w_fields.fraction_ag.sum() * np.prod(spacing_A)),
        "max_total_fill_fraction": float(total_fill.max()),
        "source_power_on_state_eV_ps": float(ag_fields.source_power_eV_ps.sum()),
    }
