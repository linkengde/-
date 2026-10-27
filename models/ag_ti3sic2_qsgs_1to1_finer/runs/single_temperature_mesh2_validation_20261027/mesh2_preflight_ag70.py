#!/usr/bin/env python3
"""Measure 1x/2x solver cost and source/energy agreement on the Ag 70 vol% mask."""

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
GEOMETRY = RUN_DIR / "qsgs_morphology_ag70vol.npz"
sys.path.insert(0, str(PHASE_DIR))
from qsgs_enthalpy_solver import QSGSEnthalpyGrid  # noqa: E402


def run_one(refinement: int, duration_ps: float) -> dict[str, object]:
    model = QSGSEnthalpyGrid(
        None,
        "A",
        native_geometry_npz=GEOMETRY,
        source_geometry_npz=GEOMETRY,
        tile_xy=1,
        mesh_refinement=refinement,
    )
    stable_dt_fs = model.stable_timestep_s() * 1.0e15
    start = time.perf_counter()
    rows = model.run(
        duration_ps,
        source_on=True,
        sample_interval_ps=duration_ps,
        track_liquid_connectivity=True,
        progress=False,
    )
    elapsed = time.perf_counter() - start
    final = rows[-1]
    result = {
        "refinement": refinement,
        "shape_xyz": list(model.temperature_shape),
        "spacing_A": model.spacing_A.tolist(),
        "ag_volume_fraction": float(model.ag_fraction.mean()),
        "source_power_W": float(model.actual_source_power_W),
        "source_method": model.source_integration_method,
        "stable_dt_fs": stable_dt_fs,
        "steps": int(round(duration_ps * 1000.0 / stable_dt_fs)),
        "run_wall_s": elapsed,
        "max_temperature_K": float(final["max_temperature_K"]),
        "ag_liquid_volume_fraction": float(final["ag_liquid_volume_fraction"]),
        "connected_liquid_cluster_count": int(final.get("connected_liquid_cluster_count", 0)),
        "injected_energy_J": float(final["injected_energy_J"]),
        "total_energy_J": float(final["total_energy_J"]),
        "energy_balance_error_J": float(final["energy_balance_error_J"]),
        "max_rss_kib_process": int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
        "fields": (
            np.asarray(model.state.temperature_K).copy(),
            np.asarray(model.state.ag_liquid_fraction).copy(),
            np.asarray(model.ag_fraction > 0.0).copy(),
        ),
    }
    del model, rows
    gc.collect()
    return result


def main() -> None:
    duration_ps = 0.001
    with np.load(GEOMETRY, allow_pickle=False) as archive:
        shape = tuple(np.asarray(archive["skeleton"]).shape)
    geometry_sha = QSGSEnthalpyGrid._sha256_file(GEOMETRY)
    base = run_one(1, duration_ps)
    fine = run_one(2, duration_ps)
    temp_1x, liquid_1x, ag_1x = base.pop("fields")
    temp_2x, liquid_2x, ag_2x = fine.pop("fields")
    coarsened = temp_2x.reshape(shape[0], 2, shape[1], 2, shape[2], 2).mean(axis=(1, 3, 5))
    liquid_volume_2x = liquid_2x.reshape(shape[0], 2, shape[1], 2, shape[2], 2).mean(axis=(1, 3, 5))
    temp_delta = coarsened - temp_1x
    liquid_union = (liquid_1x >= 0.5) | (liquid_volume_2x >= 0.5)
    liquid_union_count = int(np.count_nonzero(liquid_union))
    report = {
        "status": "AG70VOL_MESH2_COST_AND_CONSERVATION_PREFLIGHT",
        "geometry": str(GEOMETRY.relative_to(PROJECT_DIR)),
        "geometry_sha256": geometry_sha,
        "duration_ps": duration_ps,
        "composition": {"Ag": float(ag_1x.mean()), "Ti3SiC2": float(1.0 - ag_1x.mean())},
        "1x": base,
        "2x": fine,
        "source_power_relative_difference": fine["source_power_W"] / base["source_power_W"] - 1.0,
        "stable_timestep_ratio_2x_over_1x": fine["stable_dt_fs"] / base["stable_dt_fs"],
        "block_coarsened_temperature_max_abs_difference_K": float(np.max(np.abs(temp_delta[ag_1x]))),
        "block_coarsened_Ag_temperature_L1_over_rise": float(
            np.sum(ag_1x * np.abs(coarsened - temp_1x))
            / max(float(np.sum(ag_1x * np.abs(temp_1x - 300.0))), 1.0e-300)
        ),
        "liquid_mask_block_volume_IoU": (
            float(np.count_nonzero((liquid_1x >= 0.5) & (liquid_volume_2x >= 0.5)) / liquid_union_count)
            if liquid_union_count
            else None
        ),
        "liquid_mask_comparison_note": (
            "Both resolutions have no liquid voxels at this 1 fs preflight, so IoU is undefined."
            if not liquid_union_count
            else "Liquid masks are compared after block-averaging 2x cell occupancy into each 1x cell."
        ),
        "runtime_ratio_2x_over_1x": fine["run_wall_s"] / max(base["run_wall_s"], 1.0e-12),
        "linear_0p2ps_2x_work_estimate_seconds": fine["run_wall_s"] * 200.0,
        "estimate_warning": "A 1 fs short run estimates throughput only, not melt onset or mesh convergence; long-run thermal and tracking overhead may differ.",
        "limitations": [
            "The run begins uniformly at 300 K with Case A boundary conditions and candidate properties.",
            "The Ti3SiC2 skeleton has zero heat capacity and conductivity in Case A.",
            "Candidate property values remain unsourced first-pass inputs.",
        ],
    }
    out = RUN_DIR / "mesh2_preflight_ag70.json"
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
