#!/usr/bin/env python3
"""Matched full-domain 20 fs comparison of native and component Ag graphs."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
import time
import resource

import numpy as np

MODEL_DIR = Path(__file__).resolve().parents[1]
PHASE_DIR = Path(__file__).resolve().parent
if str(MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(MODEL_DIR))
if str(PHASE_DIR) not in sys.path:
    sys.path.insert(0, str(PHASE_DIR))

from qsgs_component_ttm_backend import QSGSComponentTTMBackend  # noqa: E402


GEOMETRY = MODEL_DIR / "qsgs_morphology.npz"
ELECTRON_TABLE = MODEL_DIR / "ttm_data" / "Ag_TTM_hybrid_figure_fit_300_20000K.csv"
DEFAULT_OUTPUT = PHASE_DIR / "nested_grid_diagnostics" / "qsgs_component_full_domain_20fs_match.json"
COARSE_GROUP = (2, 1, 1)
NATIVE_GROUP = (1, 1, 1)


def _fields_on_native_voxels(backend: QSGSComponentTTMBackend) -> dict[str, np.ndarray]:
    ag = backend.source_ag
    node_ids = backend.graph.component_id[ag]
    return {
        "electron_temperature_K": backend.electron_temperature_K[node_ids].copy(),
        "lattice_temperature_K": backend.lattice_temperature_K[node_ids].copy(),
        "electron_excess_energy_eV_A3": backend.parameters.ue(
            backend.electron_temperature_K[node_ids]
        ),
        "lattice_enthalpy_eV_A3": backend._h_vol_eV_A3(
            backend.lattice_temperature_K[node_ids],
            backend.lattice_liquid_fraction[node_ids],
        ),
        "electron_to_lattice_energy_eV_A3": (
            backend.cumulative_coupling_eV[node_ids]
            / backend.node_volume_A3[node_ids]
        ),
        "ag_liquid_fraction": backend.lattice_liquid_fraction[node_ids].copy(),
    }


def _relative_l1(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sum(np.abs(a - b)) / max(float(np.sum(np.abs(b))), 1.0e-30))


def _spatial_comparison(coarse: dict[str, np.ndarray], native: dict[str, np.ndarray]) -> dict[str, float]:
    result: dict[str, float] = {}
    for name in (
        "electron_excess_energy_eV_A3",
        "lattice_enthalpy_eV_A3",
        "electron_to_lattice_energy_eV_A3",
    ):
        result[f"{name}_relative_L1"] = _relative_l1(coarse[name], native[name])
    for name in ("electron_temperature_K", "lattice_temperature_K"):
        delta = coarse[name] - native[name]
        result[f"{name}_weighted_rms_K"] = float(np.sqrt(np.mean(delta**2)))
        result[f"{name}_max_abs_K"] = float(np.max(np.abs(delta)))
    result["ag_liquid_fraction_absolute_mean_difference"] = float(
        np.mean(np.abs(coarse["ag_liquid_fraction"] - native["ag_liquid_fraction"]))
    )
    return result


def _temporal_comparison(current: dict[str, np.ndarray], previous: dict[str, np.ndarray]) -> dict[str, float]:
    result: dict[str, float] = {}
    for name in (
        "electron_excess_energy_eV_A3",
        "lattice_enthalpy_eV_A3",
        "electron_to_lattice_energy_eV_A3",
    ):
        result[f"{name}_relative_L1"] = _relative_l1(current[name], previous[name])
    for name in ("electron_temperature_K", "lattice_temperature_K"):
        result[f"{name}_peak_change_K"] = float(
            abs(current[name].max() - previous[name].max())
        )
        result[f"{name}_peak_relative_change"] = result[f"{name}_peak_change_K"] / max(
            abs(float(previous[name].max())), 1.0
        )
    result["ag_liquid_fraction_mean_absolute_change"] = float(
        np.mean(np.abs(current["ag_liquid_fraction"] - previous["ag_liquid_fraction"]))
    )
    return result


def _run_one(
    geometry: Path,
    electron_table: Path,
    tile_xy: int,
    group: tuple[int, int, int],
    preconditioner: str,
    dt_ps: float,
    duration_ps: float,
    label: str,
) -> tuple[dict, dict[str, np.ndarray]]:
    started = time.perf_counter()
    backend = QSGSComponentTTMBackend(
        geometry,
        electron_table,
        boundary_case="A",
        peak_flux_eV_ps_A2=5.0,
        gaussian_radius_A=30.0,
        absorption_depth_A=20.0,
        tile_xy=tile_xy,
        coarse_voxels=group,
        linear_preconditioner=preconditioner,
    )
    graph_init_seconds = time.perf_counter() - started
    graph_info = {
        "node_count": int(backend.graph.voxel_count.size),
        "edge_count": int(backend.graph.edge_i.size),
        "edge_counts_xyz": np.bincount(backend.graph.edge_axis, minlength=3).tolist(),
        "native_ag_voxels": int(backend.graph.voxel_count.sum()),
        "node_to_native_voxel_ratio": float(backend.graph.voxel_count.size / backend.graph.voxel_count.sum()),
        "source_power_W": backend.source_power_W,
        "graph_build_seconds": backend.graph_build_seconds,
        "initialize_seconds": graph_init_seconds,
    }
    print(
        f"[{label}] initialized nodes={graph_info['node_count']:,}, edges={graph_info['edge_count']:,}, "
        f"P={backend.source_power_W:.12g} W; starting {duration_ps:g} ps",
        flush=True,
    )
    steps_float = duration_ps / dt_ps
    steps = int(round(steps_float))
    if not math.isclose(steps_float, steps, rel_tol=0.0, abs_tol=1.0e-10):
        raise ValueError("duration must be an integer multiple of the selected time step")
    peak_te = float(backend.electron_temperature_K.max())
    peak_tl = float(backend.lattice_temperature_K.max())
    row = None
    progress_stride = max(1, steps // 10)
    for step_index in range(steps):
        row = backend.advance(dt_ps, heat_on=True, coupling_on=True, electron_diffusion_on=True)
        peak_te = max(peak_te, float(backend.electron_temperature_K.max()))
        peak_tl = max(peak_tl, float(backend.lattice_temperature_K.max()))
        if (step_index + 1) % progress_stride == 0 or step_index + 1 == steps:
            print(
                f"[{label}] {step_index + 1}/{steps} steps, t={backend.time_ps:.6g} ps, "
                f"Te,max={peak_te:.2f} K, Tl,max={peak_tl:.2f} K, "
                f"elapsed={time.perf_counter() - started:.1f} s",
                flush=True,
            )
    assert row is not None
    fields = _fields_on_native_voxels(backend)
    energy_residual = (
        backend.total_energy_J
        - backend.initial_total_energy_J
        - backend.injected_energy_J
        - backend.electron_diffusion_energy_J
        - backend.lattice_diffusion_energy_J
    )
    graph_info.update(
        {
            "time_ps": backend.time_ps,
            "dt_fs": dt_ps * 1000.0,
            "steps": steps,
            "max_electron_temperature_K": float(backend.electron_temperature_K.max()),
            "peak_electron_temperature_over_run_K": peak_te,
            "mean_electron_temperature_K": float(
                np.average(backend.electron_temperature_K, weights=backend.node_volume_A3)
            ),
            "max_lattice_temperature_K": float(backend.lattice_temperature_K.max()),
            "peak_lattice_temperature_over_run_K": peak_tl,
            "mean_lattice_temperature_K": float(
                np.average(backend.lattice_temperature_K, weights=backend.node_volume_A3)
            ),
            "max_ag_liquid_fraction": float(backend.lattice_liquid_fraction.max()),
            "ag_liquid_volume_fraction": float(
                np.dot(backend.node_volume_A3, backend.lattice_liquid_fraction)
                / backend.node_volume_A3.sum()
            ),
            "electron_energy_J": backend.electron_energy_J,
            "lattice_enthalpy_J": backend.lattice_energy_J,
            "source_energy_J": backend.injected_energy_J,
            "electron_to_lattice_cumulative_J": backend.coupling_energy_J,
            "energy_residual_J": energy_residual,
            "energy_residual_relative_to_source": abs(energy_residual)
            / max(abs(backend.injected_energy_J), 1.0e-30),
            "last_step_source_energy_J": row["source_energy_J"],
            "last_step_electron_to_lattice_energy_J": row["electron_to_lattice_energy_J"],
            "last_diffusion_solve_seconds": backend.last_diffusion_solve_seconds,
            "last_diffusion_cg_iterations": backend.last_diffusion_cg_iterations,
            "total_diffusion_cg_iterations": backend.total_diffusion_cg_iterations,
            "total_diffusion_solve_seconds": backend.total_diffusion_solve_seconds,
            "total_coupling_solve_seconds": backend.total_coupling_solve_seconds,
            "amg_setup_seconds": backend.amg_setup_seconds,
            "linear_preconditioner": backend.linear_preconditioner,
            "run_seconds": time.perf_counter() - started,
            "process_peak_rss_GiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0 / 1024.0,
        }
    )
    del backend
    return graph_info, fields


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry", type=Path, default=GEOMETRY)
    parser.add_argument("--electron-table", type=Path, default=ELECTRON_TABLE)
    parser.add_argument("--tile-xy", type=int, default=3)
    parser.add_argument("--duration-ps", type=float, default=0.02)
    parser.add_argument("--time-steps-fs", default="0.3125,0.15625")
    parser.add_argument(
        "--preconditioner",
        choices=("jacobi", "amg"),
        default="amg",
        help="CG preconditioner; AMG requires the optional pyamg package",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    dt_fs_values = tuple(float(v.strip()) for v in args.time_steps_fs.split(","))
    if len(dt_fs_values) != 2 or any(not np.isfinite(v) or v <= 0.0 for v in dt_fs_values):
        parser.error("--time-steps-fs must specify two positive comma-separated values")

    source_power_by_group: dict[str, float] = {}
    matched_by_dt: dict[str, dict] = {}
    previous_fields: dict[str, dict[str, np.ndarray]] = {}
    previous_dt_label: str | None = None
    for dt_fs in dt_fs_values:
        dt_ps = dt_fs * 1.0e-3
        dt_label = f"{dt_fs:g}fs"
        group_results = {}
        group_fields = {}
        for name, group in (("component_2x1x1", COARSE_GROUP), ("native_voxel", NATIVE_GROUP)):
            metrics, fields = _run_one(
                args.geometry,
                args.electron_table,
                args.tile_xy,
                group,
                args.preconditioner,
                dt_ps,
                args.duration_ps,
                f"{dt_label} {name}",
            )
            source_power_by_group[name] = metrics["source_power_W"]
            group_results[name] = metrics
            group_fields[name] = fields
        component = group_fields["component_2x1x1"]
        native = group_fields["native_voxel"]
        matched_by_dt[dt_label] = {
            "component_2x1x1": group_results["component_2x1x1"],
            "native_voxel": group_results["native_voxel"],
            "component_vs_native_spatial_error": _spatial_comparison(component, native),
        }
        if previous_dt_label is not None:
            matched_by_dt[dt_label]["time_step_refinement_vs_previous"] = {
                "dt_fs_previous": float(previous_dt_label.removesuffix("fs")),
                "dt_fs_current": dt_fs,
                "component_2x1x1": _temporal_comparison(
                    component, previous_fields["component_2x1x1"]
                ),
                "native_voxel": _temporal_comparison(native, previous_fields["native_voxel"]),
            }
        previous_fields = group_fields
        previous_dt_label = dt_label
        del group_fields

    power_values = list(source_power_by_group.values())
    power_relative_difference = abs(power_values[0] - power_values[1]) / max(abs(power_values[1]), 1.0e-30)
    spatial_gate = all(
        result["component_vs_native_spatial_error"][key] <= 0.05
        for result in matched_by_dt.values()
        for key in (
            "electron_excess_energy_eV_A3_relative_L1",
            "lattice_enthalpy_eV_A3_relative_L1",
            "electron_to_lattice_energy_eV_A3_relative_L1",
        )
    )
    energy_gate = all(
        run["energy_residual_relative_to_source"] <= 1.0e-10
        for result in matched_by_dt.values()
        for run in (result["component_2x1x1"], result["native_voxel"])
    )
    timestep_refinement = list(matched_by_dt.values())[1].get("time_step_refinement_vs_previous", {})
    timestep_gate = bool(timestep_refinement) and all(
        timestep_refinement[model][key] <= 0.05
        for model in ("component_2x1x1", "native_voxel")
        for key in (
            "electron_excess_energy_eV_A3_relative_L1",
            "lattice_enthalpy_eV_A3_relative_L1",
            "electron_to_lattice_energy_eV_A3_relative_L1",
            "electron_temperature_K_peak_relative_change",
            "lattice_temperature_K_peak_relative_change",
            "ag_liquid_fraction_mean_absolute_change",
        )
    )
    source_gate = power_relative_difference <= 1.0e-12
    output = {
        "status": f"FULL_DOMAIN_QSGS_COMPONENT_TTM_MATCHED_{args.duration_ps:g}PS",
        "acceptance": {
            "source_power_match_passed": source_gate,
            "energy_ledger_passed": energy_gate,
            "spatial_L1_error_le_5_percent_passed": spatial_gate,
            "time_step_refinement_le_5_percent_passed": timestep_gate,
            "numerical_20fs_gate_passed": source_gate and energy_gate and spatial_gate and timestep_gate,
            "stable_melt_pool_observed": False,
            "note": "A passed numerical transport gate does not establish a melt pool or validate scaffold/interface heat transfer.",
        },
        "geometry_file": args.geometry.name,
        "electron_table": args.electron_table.name,
        "native_voxel_resolution_A": "approximately 3 Å source voxels",
        "tile_xy": args.tile_xy,
        "duration_ps": args.duration_ps,
        "component_coarse_voxels_xyz": list(COARSE_GROUP),
        "reference_coarse_voxels_xyz": list(NATIVE_GROUP),
        "heat_source": {
            "peak_flux_eV_ps_A2": 5.0,
            "gaussian_radius_A": 30.0,
            "gaussian_definition": "exp(-r^2/R^2)",
            "absorption_depth_A": 20.0,
            "source_mask": "native binary Ag voxels; exact source integrals then summed to component nodes",
            "source_center": "centre of full 3x3 periodic domain in x/y, depth from top z surface",
            "component_source_power_W": source_power_by_group["component_2x1x1"],
            "native_source_power_W": source_power_by_group["native_voxel"],
            "relative_power_difference": power_relative_difference,
        },
        "boundary_conditions": "periodic x/y including outer supercell seams; insulated z; boundary A; no scaffold heat storage or cross-phase transport",
        "time_steps_fs": list(dt_fs_values),
        "linear_preconditioner": args.preconditioner,
        "results_by_time_step": matched_by_dt,
        "interpretation": "full 3x3 QSGS domain numerical transport comparison only; candidate properties; not a physical interface model",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
