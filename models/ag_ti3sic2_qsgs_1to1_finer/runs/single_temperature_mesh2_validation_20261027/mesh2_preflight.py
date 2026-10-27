"""Benchmark native 1x versus 2x refinement without changing QSGS geometry."""

from __future__ import annotations

import gc
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np


RUN_DIR = Path(__file__).resolve().parent
PROJECT_DIR = RUN_DIR.parents[1]
PHASE_DIR = PROJECT_DIR / "phase_field_qsgs"
GEOMETRY = PROJECT_DIR / "qsgs_morphology.npz"
sys.path.insert(0, str(PHASE_DIR))
from qsgs_enthalpy_solver import QSGSEnthalpyGrid  # noqa: E402


def run_one(refinement: int, duration_ps: float) -> dict[str, object]:
    wall0 = time.perf_counter()
    model = QSGSEnthalpyGrid(
        None,
        "A",
        native_geometry_npz=GEOMETRY,
        tile_xy=1,
        mesh_refinement=refinement,
    )
    init_s = time.perf_counter() - wall0
    dt_fs = model.stable_timestep_s() * 1.0e15
    wall1 = time.perf_counter()
    rows = model.run(
        duration_ps,
        source_on=True,
        sample_interval_ps=duration_ps,
        track_liquid_connectivity=True,
    )
    run_s = time.perf_counter() - wall1
    final = rows[-1]
    temperature = np.asarray(model.state.temperature_K).copy()
    liquid = np.asarray(model.state.ag_liquid_fraction).copy()
    ag = np.asarray(model.ag_fraction, dtype=bool).copy()
    state = {
        "refinement": refinement,
        "shape_xyz": list(model.temperature_shape),
        "spacing_A": model.spacing_A.tolist(),
        "ag_volume_fraction": float(model.ag_fraction.mean()),
        "source_power_W": float(model.actual_source_power_W),
        "stable_dt_fs": dt_fs,
        "initialize_wall_s": init_s,
        "run_wall_s": run_s,
        "max_temperature_K": float(final["max_temperature_K"]),
        "Ag_liquid_volume_fraction": float(final["ag_liquid_volume_fraction"]),
        "connected_liquid_cluster_count": int(final.get("connected_liquid_cluster_count", 0)),
        "injected_energy_J": float(final["injected_energy_J"]),
        "total_energy_J": float(final["total_energy_J"]),
        "energy_balance_error_J": float(final["energy_balance_error_J"]),
        "max_rss_kib_process": int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
    }
    del model
    gc.collect()
    state["fields"] = (temperature, liquid, ag)
    return state


def main() -> None:
    with np.load(GEOMETRY, allow_pickle=False) as archive:
        base_shape = tuple(np.asarray(archive["skeleton"]).shape)
        geometry_sha256 = QSGSEnthalpyGrid._sha256_file(GEOMETRY)
    duration_ps = 0.001
    base = run_one(1, duration_ps)
    fine = run_one(2, duration_ps)
    temp_1x, liquid_1x, ag_1x = base.pop("fields")
    temp_2x, liquid_2x, ag_2x = fine.pop("fields")
    geometry_exact = bool(
        np.array_equal(ag_2x, np.repeat(np.repeat(np.repeat(ag_1x, 2, 0), 2, 1), 2, 2))
    )
    blocked_temperature = temp_2x.reshape(
        base_shape[0], 2, base_shape[1], 2, base_shape[2], 2
    ).mean(axis=(1, 3, 5))
    temperature_delta = blocked_temperature - temp_1x
    threshold_1x = (liquid_1x >= 0.5) & ag_1x
    threshold_2x = (liquid_2x >= 0.5) & ag_2x
    threshold_2x_volume = threshold_2x.reshape(
        base_shape[0], 2, base_shape[1], 2, base_shape[2], 2
    ).mean(axis=(1, 3, 5))
    threshold_union = int(np.count_nonzero(threshold_1x | (threshold_2x_volume >= 0.5)))
    report = {
        "status": "MESH2_SHORT_PREFLIGHT_ONLY",
        "geometry_sha256": geometry_sha256,
        "input_mask_shape_xyz": list(base_shape),
        "short_heating_duration_ps": duration_ps,
        "same_existing_volume_composition_only": True,
        "1x": base,
        "2x": fine,
        "geometry_is_exact_voxel_subdivision": geometry_exact,
        "source_power_relative_difference": fine["source_power_W"] / base["source_power_W"] - 1.0,
        "stable_timestep_ratio_2x_over_1x": fine["stable_dt_fs"] / base["stable_dt_fs"],
        "coarsened_2x_max_abs_temperature_difference_K": float(
            np.max(np.abs(temperature_delta[ag_1x]))
        ),
        "coarsened_2x_Ag_weighted_L1_temperature_difference_over_300K_rise": float(
            np.sum(ag_1x * np.abs(temperature_delta))
            / max(float(np.sum(ag_1x * np.abs(temp_1x - 300.0))), 1.0e-300)
        ),
        "threshold_liquid_mask_IoU_after_block_volume_average": (
            float(np.count_nonzero(threshold_1x & (threshold_2x_volume >= 0.5)) / threshold_union)
            if threshold_union else 1.0
        ),
        "runtime_ratio_2x_over_1x": fine["run_wall_s"] / max(base["run_wall_s"], 1e-12),
        "work_estimate_for_0p2ps_2x_seconds": fine["run_wall_s"] * (0.2 / duration_ps),
        "estimate_warning": (
            "The 1 fs benchmark measures setup/short-run throughput only. The 0.2 ps estimate is linear "
            "in elapsed time and does not include thermal-state-dependent overhead or filesystem variability. "
            "No melt-onset or mesh-convergence conclusion follows from this benchmark."
        ),
        "interpretation": (
            "This is a numerical and cost preflight on the existing QSGS Ag volume fraction; it does not "
            "alter or validate the experimental target composition."
        ),
    }
    (RUN_DIR / "mesh2_preflight.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
