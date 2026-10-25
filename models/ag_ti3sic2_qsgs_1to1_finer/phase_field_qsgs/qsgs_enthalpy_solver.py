"""First-pass finite-volume enthalpy solver on the projected QSGS geometry.

This is a single-temperature continuum baseline. It is intentionally separate
from the repository's TTM-MD work: heat is deposited directly into the Ag
lattice enthalpy, so this module does not model electron/lattice nonequilibrium.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from thermophysical_model import AgEnthalpy, Ti3SiC2ThermalProperties


DATA_DIR = Path(__file__).resolve().parent
ANGSTROM3_TO_M3 = 1.0e-30
EV_PER_PS_ANGSTROM2_TO_W_M2 = 1.602176634e-19 / 1.0e-12 / 1.0e-20


@dataclass(frozen=True)
class ThermalState:
    temperature_K: np.ndarray
    ag_liquid_fraction: np.ndarray


class QSGSEnthalpyGrid:
    """Conservative cell-centered enthalpy method for the two limiting BCs."""

    def __init__(
        self,
        projected_npz: Path,
        boundary_case: str,
        peak_flux_eV_ps_A2: float = 5.0,
        gaussian_radius_A: float = 30.0,
        absorption_depth_A: float = 20.0,
        source_center_xy_A: tuple[float, float] | None = None,
    ) -> None:
        if boundary_case not in {"A", "B"}:
            raise ValueError("boundary_case must be A (adiabatic scaffold) or B (ideal contact)")
        if peak_flux_eV_ps_A2 < 0 or gaussian_radius_A <= 0 or absorption_depth_A <= 0:
            raise ValueError("Source flux must be nonnegative and lengths must be positive")
        with np.load(projected_npz) as archive:
            self.ag_fraction = np.asarray(archive["ag_fraction"], dtype=float)
            self.sc_fraction = np.asarray(archive["skeleton_fraction"], dtype=float)
            self.size_A = np.asarray(archive["size_A"], dtype=float)
            self.spacing_A = np.asarray(archive["spacing_A"], dtype=float)
        if self.ag_fraction.ndim != 3 or self.sc_fraction.shape != self.ag_fraction.shape:
            raise ValueError("Projected phase fractions must be matching three-dimensional arrays")
        if not np.allclose(self.ag_fraction + self.sc_fraction, 1.0, atol=2.0e-12):
            raise ValueError("Projected Ag and scaffold fractions must sum to one")
        if np.any(self.ag_fraction < 0.0) or np.any(self.sc_fraction < 0.0):
            raise ValueError("Phase fractions cannot be negative")

        self.boundary_case = boundary_case
        self.ag = AgEnthalpy()
        self.sc = Ti3SiC2ThermalProperties()
        self.peak_flux_W_m2 = peak_flux_eV_ps_A2 * EV_PER_PS_ANGSTROM2_TO_W_M2
        self.gaussian_radius_A = float(gaussian_radius_A)
        self.absorption_depth_A = float(absorption_depth_A)
        self.cell_volume_m3 = float(np.prod(self.spacing_A) * ANGSTROM3_TO_M3)
        self.temperature_shape = self.ag_fraction.shape
        x = (np.arange(self.temperature_shape[0]) + 0.5) * self.spacing_A[0]
        y = (np.arange(self.temperature_shape[1]) + 0.5) * self.spacing_A[1]
        z = (np.arange(self.temperature_shape[2]) + 0.5) * self.spacing_A[2]
        cx, cy = source_center_xy_A or (self.size_A[0] / 2.0, self.size_A[1] / 2.0)
        rr2 = (x[:, None, None] - cx) ** 2 + (y[None, :, None] - cy) ** 2
        depth = self.size_A[2] - z[None, None, :]
        spatial = np.exp(-rr2 / self.gaussian_radius_A**2) * np.exp(-depth / self.absorption_depth_A)
        self.source_W_m3 = (
            self.ag_fraction
            * self.peak_flux_W_m2
            / (self.absorption_depth_A * 1.0e-10)
            * spatial
        )
        self.actual_source_power_W = float(self.source_W_m3.sum() * self.cell_volume_m3)
        self._solid_at_tm_J_m3 = self._solid_volumetric_enthalpy(
            np.full(self.temperature_shape, self.ag.melting_temperature_K)
        )
        self._latent_volumetric_J_m3 = (
            self.ag_fraction * self.ag.reference_density_kg_m3 * self.ag.latent_heat_J_kg
        )
        self._solid_capacity_at_300_J_m3_K = (
            self.ag_fraction * self.ag.reference_density_kg_m3 * self.ag.solid_cp_J_kg_K[0]
        )
        self._solid_capacity_slope_J_m3_K2 = self.ag_fraction * self.ag.reference_density_kg_m3 * (
            (self.ag.solid_cp_J_kg_K[1] - self.ag.solid_cp_J_kg_K[0])
            / (self.ag.solid_temperature_K[1] - self.ag.solid_temperature_K[0])
        )
        if self.boundary_case == "B":
            sc_cp0 = self.sc.cp_mass_J_kg_K[0]
            sc_slope = (self.sc.cp_mass_J_kg_K[1] - sc_cp0) / (
                self.sc.temperature_K[1] - self.sc.temperature_K[0]
            )
            self._solid_capacity_at_300_J_m3_K += self.sc_fraction * self.sc.density_kg_m3 * sc_cp0
            self._solid_capacity_slope_J_m3_K2 += self.sc_fraction * self.sc.density_kg_m3 * sc_slope
        self.state = self.state_from_enthalpy_density(np.zeros(self.temperature_shape, dtype=float))
        self.enthalpy_density_J_m3 = np.zeros(self.temperature_shape, dtype=float)
        self.time_s = 0.0
        self.injected_energy_J = 0.0

    @property
    def melting_temperature_K(self) -> float:
        return self.ag.melting_temperature_K

    def _solid_volumetric_enthalpy(self, temperature_K: np.ndarray) -> np.ndarray:
        """Mixture sensible enthalpy relative to 300 K, below Ag melting."""
        temp = np.asarray(temperature_K, dtype=float)
        ag_t = np.minimum(temp, self.ag.melting_temperature_K)
        h_ag = np.asarray(self.ag._integrate_solid_cp(ag_t))
        ag_part = self.ag_fraction * self.ag.reference_density_kg_m3 * h_ag
        if self.boundary_case == "A":
            return ag_part

        cp0 = self.sc.cp_mass_J_kg_K[0]
        slope = (self.sc.cp_mass_J_kg_K[1] - cp0) / (self.sc.temperature_K[1] - self.sc.temperature_K[0])
        delta = temp - self.sc.temperature_K[0]
        h_sc = cp0 * delta + 0.5 * slope * delta**2
        h_sc = np.where(temp < self.sc.temperature_K[0], cp0 * delta, h_sc)
        return ag_part + self.sc_fraction * self.sc.density_kg_m3 * h_sc

    def enthalpy_density_from_state(self, temperature_K: np.ndarray, liquid_fraction: np.ndarray) -> np.ndarray:
        temp, liquid = np.broadcast_arrays(np.asarray(temperature_K, dtype=float), np.asarray(liquid_fraction, dtype=float))
        if np.any((liquid < 0.0) | (liquid > 1.0)):
            raise ValueError("Ag liquid fraction must lie in [0,1]")
        tm = self.melting_temperature_K
        solid_at_tm = self._solid_at_tm_J_m3
        latent_vol = self._latent_volumetric_J_m3
        h_solid = self._solid_volumetric_enthalpy(temp)
        h_plateau = solid_at_tm + latent_vol * liquid
        h_liquid = solid_at_tm + latent_vol + self.ag_fraction * self.ag.reference_density_kg_m3 * float(self.ag.liquid_cp_J_kg_K[0]) * (temp - tm)
        if self.boundary_case == "B":
            dtemp = temp - tm
            slope_sc = (self.sc.cp_mass_J_kg_K[1] - self.sc.cp_mass_J_kg_K[0]) / (self.sc.temperature_K[1] - self.sc.temperature_K[0])
            cp_sc_tm = float(self.sc.cp_J_kg_K(tm))
            h_liquid += self.sc_fraction * self.sc.density_kg_m3 * (cp_sc_tm * dtemp + 0.5 * slope_sc * dtemp**2)
        out = np.where(temp < tm, h_solid, np.where(temp > tm, h_liquid, h_plateau))
        return out

    def state_from_enthalpy_density(self, enthalpy_density_J_m3: np.ndarray) -> ThermalState:
        h = np.asarray(enthalpy_density_J_m3, dtype=float)
        if h.shape != self.temperature_shape:
            h = np.broadcast_to(h, self.temperature_shape)
        temp = np.full(self.temperature_shape, 300.0, dtype=float)
        liquid = np.zeros(self.temperature_shape, dtype=float)
        ag_rho = self.ag.reference_density_kg_m3
        c0 = self._solid_capacity_at_300_J_m3_K
        slope = self._solid_capacity_slope_J_m3_K2

        below_300 = h < 0.0
        has_capacity = c0 > 0.0
        temp[below_300 & has_capacity] = 300.0 + h[below_300 & has_capacity] / c0[below_300 & has_capacity]
        solid_at_tm = self._solid_at_tm_J_m3
        solid = (h >= 0.0) & (h < solid_at_tm) & has_capacity
        target = h[solid]
        c = c0[solid]
        sl = slope[solid]
        if target.size:
            delta = target / c
            curved = np.abs(sl) > 1.0e-18
            delta[curved] = 2.0 * target[curved] / (c[curved] + np.sqrt(c[curved] ** 2 + 2.0 * sl[curved] * target[curved]))
            temp[solid] = 300.0 + delta

        latent_vol = self._latent_volumetric_J_m3
        liquidus_enthalpy = solid_at_tm + latent_vol
        plateau = (h >= solid_at_tm) & (h <= liquidus_enthalpy) & (latent_vol > 0.0)
        temp[plateau] = self.melting_temperature_K
        liquid[plateau] = np.clip((h[plateau] - solid_at_tm[plateau]) / latent_vol[plateau], 0.0, 1.0)
        no_ag_at_tm = (latent_vol <= 0.0) & (solid_at_tm > 0.0) & np.isclose(
            h, solid_at_tm, rtol=2.0e-14, atol=1.0e-8
        )
        temp[no_ag_at_tm] = self.melting_temperature_K

        above = (h > liquidus_enthalpy) & has_capacity
        if np.any(above):
            dH = h[above] - liquidus_enthalpy[above]
            cp_liquid = float(self.ag.liquid_cp_J_kg_K[0])
            c_tm = self.ag_fraction[above] * ag_rho * cp_liquid
            slope_liquid = np.zeros_like(c_tm)
            if self.boundary_case == "B":
                sc_cp_tm = float(self.sc.cp_J_kg_K(self.melting_temperature_K))
                sc_slope = (self.sc.cp_mass_J_kg_K[1] - self.sc.cp_mass_J_kg_K[0]) / (self.sc.temperature_K[1] - self.sc.temperature_K[0])
                c_tm += self.sc_fraction[above] * self.sc.density_kg_m3 * sc_cp_tm
                slope_liquid = self.sc_fraction[above] * self.sc.density_kg_m3 * sc_slope
            delta = dH / c_tm
            curved = np.abs(slope_liquid) > 1.0e-18
            delta[curved] = 2.0 * dH[curved] / (c_tm[curved] + np.sqrt(c_tm[curved] ** 2 + 2.0 * slope_liquid[curved] * dH[curved]))
            temp[above] = self.melting_temperature_K + delta
            liquid[above] = 1.0

        if self.boundary_case == "B" and np.any(temp[self.sc_fraction > 0.0] > self.sc.temperature_K[-1] + 1.0e-9):
            raise ValueError("Ti3SiC2 candidate table exceeds 1473.15 K; add high-temperature data before continuing")
        return ThermalState(temp, liquid)

    def conductivity_W_m_K(self, state: ThermalState | None = None) -> np.ndarray:
        state = self.state if state is None else state
        ag_k = np.asarray(self.ag.conductivity_W_m_K(state.temperature_K, state.ag_liquid_fraction))
        result = self.ag_fraction * ag_k
        if self.boundary_case == "B":
            sc_k = np.asarray(self.sc.conductivity_W_m_K(state.temperature_K))
            result = result + self.sc_fraction * sc_k
        return result

    def sensible_capacity_J_m3_K(self, state: ThermalState | None = None) -> np.ndarray:
        state = self.state if state is None else state
        ag_solid = np.asarray(self.ag.cp_solid_J_kg_K(np.minimum(state.temperature_K, self.melting_temperature_K)))
        ag_liquid = np.asarray(self.ag.cp_liquid_J_kg_K(np.maximum(state.temperature_K, self.melting_temperature_K)))
        ag_cp = (1.0 - state.ag_liquid_fraction) * ag_solid + state.ag_liquid_fraction * ag_liquid
        capacity = self.ag_fraction * self.ag.reference_density_kg_m3 * ag_cp
        if self.boundary_case == "B":
            capacity += self.sc_fraction * self.sc.density_kg_m3 * np.asarray(self.sc.cp_J_kg_K(state.temperature_K))
        return capacity

    def stable_timestep_s(self, safety: float = 0.45) -> float:
        """Conservative explicit diffusion limit using the current state."""
        k = self.conductivity_W_m_K()
        capacity = self.sensible_capacity_J_m3_K()
        inverse_spacing2 = 1.0 / (self.spacing_A * 1.0e-10) ** 2
        face_sum = np.zeros(self.temperature_shape, dtype=float)
        for axis in range(3):
            right_k = np.roll(k, -1, axis=axis) if axis < 2 else np.concatenate((k[:, :, 1:], np.zeros_like(k[:, :, :1])), axis=2)
            left_k = np.roll(k, 1, axis=axis) if axis < 2 else np.concatenate((np.zeros_like(k[:, :, :1]), k[:, :, :-1]), axis=2)
            right = np.divide(2.0 * k * right_k, k + right_k, out=np.zeros_like(k), where=(k + right_k) > 0.0)
            left = np.divide(2.0 * k * left_k, k + left_k, out=np.zeros_like(k), where=(k + left_k) > 0.0)
            face_sum += (right + left) * inverse_spacing2[axis]
        rate = np.divide(face_sum, capacity, out=np.zeros_like(face_sum), where=capacity > 0.0)
        max_rate = float(rate.max())
        if max_rate <= 0.0:
            return np.inf
        return safety / max_rate

    def diffusion_power_density_W_m3(self, temperature_K: np.ndarray, conductivity_W_m_K: np.ndarray) -> np.ndarray:
        """Pairwise-conservative periodic x/y, insulated z diffusion operator."""
        result = np.zeros(self.temperature_shape, dtype=float)
        for axis in (0, 1):
            dx = self.spacing_A[axis] * 1.0e-10
            t_right = np.roll(temperature_K, -1, axis=axis)
            k_right = np.roll(conductivity_W_m_K, -1, axis=axis)
            kface = np.divide(
                2.0 * conductivity_W_m_K * k_right,
                conductivity_W_m_K + k_right,
                out=np.zeros_like(conductivity_W_m_K),
                where=(conductivity_W_m_K + k_right) > 0.0,
            )
            result += kface * (t_right - temperature_K) / dx**2
            result += np.roll(kface, 1, axis=axis) * (np.roll(temperature_K, 1, axis=axis) - temperature_K) / dx**2
        dz = self.spacing_A[2] * 1.0e-10
        kface = np.divide(
            2.0 * conductivity_W_m_K[:, :, :-1] * conductivity_W_m_K[:, :, 1:],
            conductivity_W_m_K[:, :, :-1] + conductivity_W_m_K[:, :, 1:],
            out=np.zeros_like(conductivity_W_m_K[:, :, :-1]),
            where=(conductivity_W_m_K[:, :, :-1] + conductivity_W_m_K[:, :, 1:]) > 0.0,
        )
        flux_term = kface * (temperature_K[:, :, 1:] - temperature_K[:, :, :-1]) / dz**2
        result[:, :, :-1] += flux_term
        result[:, :, 1:] -= flux_term
        return result

    def advance(
        self,
        dt_s: float,
        source_on: bool = True,
        enforce_stability: bool = True,
    ) -> dict[str, float]:
        if dt_s <= 0.0 or not np.isfinite(dt_s):
            raise ValueError("dt_s must be finite and positive")
        if enforce_stability and dt_s > self.stable_timestep_s() * (1.0 + 1.0e-10):
            raise ValueError("Requested timestep exceeds the current explicit diffusion stability limit")
        k = self.conductivity_W_m_K()
        diffusion = self.diffusion_power_density_W_m3(self.state.temperature_K, k)
        source = self.source_W_m3 if source_on else 0.0
        self.enthalpy_density_J_m3 = self.enthalpy_density_J_m3 + dt_s * (diffusion + source)
        self.time_s += dt_s
        if source_on:
            self.injected_energy_J += self.actual_source_power_W * dt_s
        self.state = self.state_from_enthalpy_density(self.enthalpy_density_J_m3)
        total_energy = float(self.enthalpy_density_J_m3.sum() * self.cell_volume_m3)
        error = total_energy - self.injected_energy_J
        return {
            "time_ps": self.time_s * 1.0e12,
            "dt_fs": dt_s * 1.0e15,
            "max_temperature_K": float(self.state.temperature_K.max()),
            "mean_ag_temperature_K": self._ag_weighted_mean_temperature(),
            "max_ag_liquid_fraction": float(self.state.ag_liquid_fraction.max()),
            "ag_liquid_volume_fraction": self._liquid_ag_volume_fraction(),
            "total_energy_J": total_energy,
            "injected_energy_J": self.injected_energy_J,
            "energy_balance_error_J": error,
        }

    def _ag_weighted_mean_temperature(self) -> float:
        denominator = float(self.ag_fraction.sum())
        if denominator <= 0.0:
            return 300.0
        return float(np.sum(self.ag_fraction * self.state.temperature_K) / denominator)

    def _liquid_ag_volume_fraction(self) -> float:
        total = float(self.ag_fraction.sum())
        return float(np.sum(self.ag_fraction * self.state.ag_liquid_fraction) / total) if total else 0.0

    def run(self, duration_ps: float, source_on: bool, sample_interval_ps: float = 0.1) -> list[dict[str, float]]:
        if duration_ps <= 0.0 or sample_interval_ps <= 0.0:
            raise ValueError("Run duration and sample interval must be positive")
        end_s = self.time_s + duration_ps * 1.0e-12
        sample_interval_s = sample_interval_ps * 1.0e-12
        next_sample = self.time_s
        # Over the candidate range k(T)/Cp(T) decreases from its 300 K value
        # in both phases; a fixed initial stable step is therefore conservative
        # and avoids recomputing the full-grid diffusion bound every substep.
        stable_dt = self.stable_timestep_s()
        rows: list[dict[str, float]] = []
        while self.time_s < end_s - 1.0e-24:
            dt = min(stable_dt, end_s - self.time_s, next_sample + sample_interval_s - self.time_s)
            if dt <= 1.0e-25:
                next_sample += sample_interval_s
                continue
            row = self.advance(dt, source_on=source_on, enforce_stability=False)
            if self.time_s >= next_sample + sample_interval_s - 1.0e-24 or self.time_s >= end_s - 1.0e-24:
                rows.append(row)
                next_sample += sample_interval_s
        return rows


def write_csv(path: Path, rows: list[dict[str, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spacing", choices=("5A", "10A"), default="10A")
    parser.add_argument("--boundary", choices=("A", "B"), required=True)
    parser.add_argument("--heat-ps", type=float, default=1.0)
    parser.add_argument("--cool-ps", type=float, default=1.0)
    parser.add_argument("--sample-ps", type=float, default=0.1)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    path = DATA_DIR / f"qsgs_phase_fractions_{args.spacing}.npz"
    model = QSGSEnthalpyGrid(path, args.boundary)
    rows = model.run(args.heat_ps, source_on=True, sample_interval_ps=args.sample_ps)
    rows += model.run(args.cool_ps, source_on=False, sample_interval_ps=args.sample_ps)
    output = args.output or DATA_DIR / f"qsgs_enthalpy_{args.spacing}_BC{args.boundary}.csv"
    write_csv(output, rows)
    summary = {
        "status": "SINGLE_TEMPERATURE_ENTHALPY_EXPLORATORY_RUN",
        "boundary_case": args.boundary,
        "spacing_label": args.spacing,
        "grid_shape_xyz": list(model.temperature_shape),
        "actual_spacing_A": model.spacing_A.tolist(),
        "heat_duration_ps": args.heat_ps,
        "cool_duration_ps": args.cool_ps,
        "total_time_ps": model.time_s * 1.0e12,
        "actual_source_power_W": model.actual_source_power_W,
        "injected_energy_J": model.injected_energy_J,
        "final_total_enthalpy_J": float(model.enthalpy_density_J_m3.sum() * model.cell_volume_m3),
        "energy_balance_error_J": float(model.enthalpy_density_J_m3.sum() * model.cell_volume_m3 - model.injected_energy_J),
        "max_temperature_K": float(model.state.temperature_K.max()),
        "ag_liquid_volume_fraction": model._liquid_ag_volume_fraction(),
        "max_ag_liquid_fraction": float(model.state.ag_liquid_fraction.max()),
        "limitations": [
            "Heat is deposited directly into lattice enthalpy; this is not TTM or TTM-MD.",
            "Candidate properties are first-pass values pending source-table archival and citation review.",
            "Boundary A/B are limiting assumptions, not measured Ag/Ti3SiC2 interface heat transfer.",
            "External z faces are insulated and x/y are periodic; no free-surface radiation or convection is included.",
        ],
        "timeseries_csv": str(output.name),
    }
    summary_path = output.with_suffix(".json")
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
