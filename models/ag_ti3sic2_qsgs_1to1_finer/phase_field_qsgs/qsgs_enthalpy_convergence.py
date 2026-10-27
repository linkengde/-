"""Serial mesh and timestep convergence screen for the native QSGS solver."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import numpy as np

from qsgs_enthalpy_solver import QSGSEnthalpyGrid


DATA_DIR = Path(__file__).resolve().parent


def _run_case(
    geometry_path: Path,
    boundary: str,
    duration_ps: float,
    refinement: int,
    maximum_dt_s: float,
) -> tuple[QSGSEnthalpyGrid, dict[str, float], np.ndarray]:
    model = QSGSEnthalpyGrid(
        None,
        boundary,
        native_geometry_npz=geometry_path,
        mesh_refinement=refinement,
    )
    stable_dt_s = model.stable_timestep_s()
    if maximum_dt_s > stable_dt_s * (1.0 + 1.0e-10):
        raise ValueError(
            f"Requested dt {maximum_dt_s * 1e15:.6g} fs exceeds refinement-{refinement} "
            f"stability limit {stable_dt_s * 1e15:.6g} fs"
        )
    start = time.monotonic()
    rows = model.run(
        duration_ps,
        source_on=True,
        sample_interval_ps=duration_ps,
        maximum_dt_s=maximum_dt_s,
    )
    elapsed = time.monotonic() - start
    total_enthalpy = float(model.enthalpy_density_J_m3.sum() * model.cell_volume_m3)
    ag_mask = model.ag_fraction > 0.0
    result = {
        "mesh_refinement": int(refinement),
        "grid_shape_xyz": list(model.temperature_shape),
        "spacing_A": model.spacing_A.tolist(),
        "stable_timestep_fs": stable_dt_s * 1.0e15,
        "maximum_timestep_fs": maximum_dt_s * 1.0e15,
        "step_count": int(math.ceil(duration_ps * 1.0e-12 / maximum_dt_s)),
        "elapsed_s": elapsed,
        "source_power_W": model.actual_source_power_W,
        "injected_energy_J": model.injected_energy_J,
        "final_enthalpy_J": total_enthalpy,
        "energy_balance_error_J": total_enthalpy - model.injected_energy_J,
        "max_ag_temperature_K": float(model.state.temperature_K[ag_mask].max()),
        "ag_weighted_mean_temperature_K": model._ag_weighted_mean_temperature(),
        "max_ag_liquid_fraction": float(model.state.ag_liquid_fraction.max()),
    }
    if rows:
        result["reported_dt_fs_at_sample"] = float(rows[-1]["dt_fs"])
    return model, result, model.state.temperature_K.copy()


def _weighted_relative_difference(
    reference: np.ndarray,
    comparison: np.ndarray,
    weights: np.ndarray,
) -> float:
    numerator = float(np.sum(weights * np.abs(reference - comparison)))
    denominator = float(np.sum(weights * np.abs(reference - 300.0)))
    return numerator / max(denominator, np.finfo(float).tiny)


def compare_native_mesh_and_timestep(
    duration_ps: float,
    boundary: str,
    base_dt_fs: float | None = None,
    output: Path | None = None,
) -> dict[str, object]:
    geometry_path = DATA_DIR.parent / "qsgs_morphology.npz"
    coarse_probe = QSGSEnthalpyGrid(None, boundary, native_geometry_npz=geometry_path)
    fine_probe = QSGSEnthalpyGrid(None, boundary, native_geometry_npz=geometry_path, mesh_refinement=2)
    minimum_stable_fs = min(coarse_probe.stable_timestep_s(), fine_probe.stable_timestep_s()) * 1.0e15
    chosen_dt_fs = min(0.8 * minimum_stable_fs, base_dt_fs) if base_dt_fs is not None else 0.8 * minimum_stable_fs
    if chosen_dt_fs <= 0.0:
        raise ValueError("base_dt_fs must be positive")

    coarse_dt, coarse_half_dt, fine_dt = chosen_dt_fs * 1.0e-15, chosen_dt_fs * 0.5e-15, chosen_dt_fs * 1.0e-15
    cases: dict[str, dict[str, float]] = {}
    coarse_model, cases["native_mesh_base_dt"], coarse_field = _run_case(
        geometry_path, boundary, duration_ps, 1, coarse_dt
    )
    _, cases["native_mesh_half_dt"], half_dt_field = _run_case(
        geometry_path, boundary, duration_ps, 1, coarse_half_dt
    )
    fine_model, cases["native_mesh_refined_base_dt"], fine_field = _run_case(
        geometry_path, boundary, duration_ps, 2, fine_dt
    )

    refined_common_grid = fine_field.reshape(
        coarse_model.temperature_shape[0],
        2,
        coarse_model.temperature_shape[1],
        2,
        coarse_model.temperature_shape[2],
        2,
    ).mean(axis=(1, 3, 5))
    weights = coarse_model.ag_fraction
    timestep_field_l1 = _weighted_relative_difference(coarse_field, half_dt_field, weights)
    mesh_field_l1 = _weighted_relative_difference(coarse_field, refined_common_grid, weights)
    timestep_peak_delta_K = float(np.max(np.abs(coarse_field[weights > 0] - half_dt_field[weights > 0])))
    mesh_peak_delta_K = float(np.max(np.abs(coarse_field[weights > 0] - refined_common_grid[weights > 0])))
    peak_rise_K = max(float(cases["native_mesh_base_dt"]["max_ag_temperature_K"]) - 300.0, 1.0e-30)
    tolerance = {
        "relative_ag_weighted_field_L1_timestep_halving_max": 1.0e-3,
        "relative_ag_weighted_field_L1_mesh_refinement_1_to_2_max": 1.0e-2,
        "maximum_peak_temperature_difference_as_peak_rise_max": 2.0e-2,
    }
    report: dict[str, object] = {
        "status": "NATIVE_QSGS_NUMERICAL_CONVERGENCE_SCREEN",
        "model": "single-temperature explicit finite-volume enthalpy conduction",
        "geometry_input": "existing qsgs_morphology.npz; each native voxel is uniformly subdivided for mesh refinement",
        "source_integration": "exact separable Gaussian/exponential integral over original binary Ag voxels",
        "boundary_case": boundary,
        "initial_temperature_K": 300.0,
        "duration_ps": duration_ps,
        "lateral_tile_xy": 1,
        "chosen_matched_dt_fs": chosen_dt_fs,
        "screening_tolerance": tolerance,
        "interpretation": "Working numerical gates for this exploratory start only; they are not universal criteria or paper-grade validation.",
        "runs": cases,
        "relative_ag_weighted_field_L1_timestep_halving": timestep_field_l1,
        "relative_ag_weighted_field_L1_mesh_refinement_1_to_2": mesh_field_l1,
        "max_ag_temperature_difference_timestep_K": timestep_peak_delta_K,
        "max_ag_temperature_difference_mesh_K": mesh_peak_delta_K,
        "max_temperature_difference_as_base_peak_rise_timestep": timestep_peak_delta_K / peak_rise_K,
        "max_temperature_difference_as_base_peak_rise_mesh": mesh_peak_delta_K / peak_rise_K,
        "screening_checks": {
            "timestep_field_l1_within_gate": timestep_field_l1 <= tolerance["relative_ag_weighted_field_L1_timestep_halving_max"],
            "mesh_field_l1_within_gate": mesh_field_l1 <= tolerance["relative_ag_weighted_field_L1_mesh_refinement_1_to_2_max"],
            "mesh_peak_difference_within_gate": mesh_peak_delta_K / peak_rise_K <= tolerance["maximum_peak_temperature_difference_as_peak_rise_max"],
        },
        "all_runs_energy_balance_error_J_below_1e-12_relative": all(
            abs(run["energy_balance_error_J"])
            <= 1.0e-12 * max(abs(run["injected_energy_J"]), 1.0e-300)
            for run in cases.values()
        ),
    }
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, indent=2) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--duration-ps", type=float, default=0.005)
    parser.add_argument("--boundary", choices=("A", "B"), default="A")
    parser.add_argument("--base-dt-fs", type=float, default=None)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    if args.duration_ps <= 0.0 or (args.base_dt_fs is not None and args.base_dt_fs <= 0.0):
        parser.error("duration and timestep must be positive")
    report = compare_native_mesh_and_timestep(args.duration_ps, args.boundary, args.base_dt_fs, args.output)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
