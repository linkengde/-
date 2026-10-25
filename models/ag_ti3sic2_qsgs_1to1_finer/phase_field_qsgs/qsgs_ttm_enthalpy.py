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

MODEL_DIR = Path(__file__).resolve().parents[1]
if str(MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(MODEL_DIR))

from ag_ttm_reference import AgParameters, AgTTMGrid
from qsgs_enthalpy_solver import DATA_DIR, QSGSEnthalpyGrid


EV_J = 1.602176634e-19
J_M3_TO_EV_A3 = 1.0e-30 / EV_J
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
        coupling_safety: float = 0.1,
    ) -> None:
        if not 0.0 < coupling_safety < 0.5:
            raise ValueError("coupling_safety must lie in (0, 0.5)")
        self.lattice = QSGSEnthalpyGrid(
            projected_npz,
            boundary_case,
            peak_flux_eV_ps_A2,
            gaussian_radius_A,
            absorption_depth_A,
        )
        # This is the cold-state, worst-case explicit lattice conduction bound
        # for the implemented k/Cp curves; it remains conservative as the
        # material heats because conductivity falls while heat capacity rises.
        self.lattice_diffusion_timestep_ps = self.lattice.stable_timestep_s() * 1.0e12
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
        self.max_g_eV_ps_A3_K = float(self.parameters.g(self.parameters.temperature_K).max())
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
        self.coupling_safety = coupling_safety
        self.time_ps = 0.0
        self.injected_energy_J = 0.0
        self.coupling_energy_J = 0.0
        self.electron_diffusion_energy_J = 0.0
        self.lattice_diffusion_energy_J = 0.0
        self.initial_total_energy_J = self.total_energy_J

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

    def _maximum_coupling_step_ps(self, coupling_on: bool) -> float:
        if not coupling_on:
            return np.inf
        active = self.lattice.ag_fraction > 0.0
        g_max = self.max_g_eV_ps_A3_K
        if g_max <= 0.0:
            return np.inf
        capacity_eV_A3_K = self.lattice.sensible_capacity_J_m3_K() * J_M3_TO_EV_A3
        rate = np.divide(
            self.lattice.ag_fraction * g_max,
            capacity_eV_A3_K,
            out=np.zeros_like(capacity_eV_A3_K),
            where=capacity_eV_A3_K > 0.0,
        )
        max_rate = float(rate.max())
        return np.inf if max_rate <= 0.0 else self.coupling_safety / max_rate

    def advance(
        self,
        dt_ps: float,
        heat_on: bool = True,
        coupling_on: bool = True,
    ) -> dict[str, float]:
        if not np.isfinite(dt_ps) or dt_ps <= 0.0:
            raise ValueError("dt_ps must be finite and positive")
        initial_total = self.total_energy_J
        source_increment_J = 0.0
        coupling_increment_J = 0.0
        diffusion_increment_J = 0.0
        lattice_diffusion_increment_J = 0.0
        t_end = self.time_ps + dt_ps
        substeps = 0
        while self.time_ps < t_end - max(1.0e-15, dt_ps * 1.0e-13):
            subdt = min(
                t_end - self.time_ps,
                self._maximum_coupling_step_ps(coupling_on),
                self.lattice_diffusion_timestep_ps,
            )
            if not np.isfinite(subdt) or subdt <= 0.0:
                raise RuntimeError("Coupled TTM step selection produced a non-positive timestep")
            lattice_temperature = self.lattice.state.temperature_K
            record = self.electron.advance(
                subdt,
                lattice_temperature_K=lattice_temperature,
                heat_on=heat_on,
                coupling_on=coupling_on,
            )
            lattice_diffusion = self.lattice.diffusion_power_density_W_m3(
                self.lattice.state.temperature_K,
                self.lattice.conductivity_W_m_K(),
            )
            transfer_J_m3 = (
                record.lattice_energy_per_cell_eV * EV_J / self.lattice.cell_volume_m3
            )
            accepted_dt_s = record.elapsed_ps * 1.0e-12
            self.lattice.enthalpy_density_J_m3 += accepted_dt_s * lattice_diffusion + transfer_J_m3
            self.lattice.state = self.lattice.state_from_enthalpy_density(
                self.lattice.enthalpy_density_J_m3
            )
            self.time_ps += record.elapsed_ps
            source_increment_J += record.source_to_electron_eV * EV_J
            coupling_increment_J += record.electron_to_lattice_eV * EV_J
            diffusion_increment_J += record.diffusion_net_eV * EV_J
            lattice_diffusion_increment_J += float(
                lattice_diffusion.sum() * accepted_dt_s * self.lattice.cell_volume_m3
            )
            substeps += record.substeps

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
        return {
            "time_ps": self.time_ps,
            "dt_ps": dt_ps,
            "electron_substeps": float(substeps),
            "max_electron_temperature_K": float(self.electron.temperature_K[active].max()),
            "mean_electron_temperature_K": float(
                np.sum(self.lattice.ag_fraction * self.electron.temperature_K) / ag_weight_sum
            ),
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
    ) -> list[dict[str, float]]:
        if duration_ps <= 0.0 or dt_ps <= 0.0:
            raise ValueError("duration_ps and dt_ps must be positive")
        sample_interval_ps = dt_ps if sample_interval_ps is None else sample_interval_ps
        if sample_interval_ps <= 0.0:
            raise ValueError("sample_interval_ps must be positive")
        end_ps = self.time_ps + duration_ps
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
            step = min(dt_ps, end_ps - self.time_ps, next_sample - self.time_ps)
            row = self.advance(step, heat_on=heat_on)
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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spacing", choices=("5A", "10A"), default="10A")
    parser.add_argument("--boundary", choices=("A", "B"), required=True)
    parser.add_argument("--heat-ps", type=float, default=0.01)
    parser.add_argument("--cool-ps", type=float, default=0.01)
    parser.add_argument("--dt-fs", type=float, default=5.0,
                        help="outer coupling step; lattice and electron diffusion subcycle internally")
    parser.add_argument("--sample-ps", type=float, default=0.005)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    projected = DATA_DIR / f"qsgs_phase_fractions_{args.spacing}.npz"
    model = QSGSDoubleTemperatureGrid(projected, args.boundary)
    dt_ps = args.dt_fs * 1.0e-3
    rows = model.run(args.heat_ps, dt_ps, heat_on=True, sample_interval_ps=args.sample_ps)
    rows += model.run(args.cool_ps, dt_ps, heat_on=False, sample_interval_ps=args.sample_ps)
    output = args.output or DATA_DIR / f"qsgs_ttm_enthalpy_{args.spacing}_BC{args.boundary}.csv"
    write_csv(output, rows)
    summary = {
        "status": "CONTINUUM_TWO_TEMPERATURE_ENTHALPY_SHORT_PILOT",
        "boundary_case": args.boundary,
        "grid_spacing_label": args.spacing,
        "grid_shape_xyz": list(model.lattice.temperature_shape),
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
        "limitations": [
            "This is a continuum two-temperature enthalpy reference solver, not atomistic TTM-MD.",
            "Ag electron properties are cross-paper candidate fits over 300–20000 K; low-temperature extension is disabled.",
            "Ti3SiC2 boundary cases A/B are limiting assumptions, not measured interface heat transfer.",
            "x/y are periodic and z faces insulated; no external substrate heat sink, convection, or radiation is included.",
            "No stable-pool controller is enabled in this short numerical pilot.",
        ],
        "time_history_csv": output.name,
    }
    output.with_suffix(".json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
