#!/usr/bin/env python3
"""在50³窄相变区做一次FFT与AMG预条件器对照，不运行完整热循环。"""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
PHASE = PROJECT / "phase_change_validation"
sys.path.insert(0, str(PHASE))
sys.path.insert(0, str(HERE))

import numpy as np

from implicit_enthalpy_fvm import BoundarySet
from material_model import CONFIG
from optimized_solver import FastCandidateMaterials, NumbaEnthalpyFVM
from profiled_advance import advance_profiled
from run_3d_phase_change import load_seed_ag, make_ag_fraction, source_field


def run_one(preconditioner: str) -> dict:
    n = 50
    edge_m = float(CONFIG["domain_edge_m"]["ultrafine_500nm"])
    config = CONFIG["numerical_heating_test"]
    dx = edge_m / n
    ag = np.ascontiguousarray(make_ag_fraction("bicontinuous_qsgs_seed1105", n, load_seed_ag()))
    rate = float(config["source_rate_W_m3_at_10um"])
    source = np.ascontiguousarray(source_field(
        ag, dx, rate, float(config["source_region_radius_fraction_of_edge"])
    ))
    materials = FastCandidateMaterials()
    initial_temperature = materials.Tm_K - 0.4
    sim = NumbaEnthalpyFVM(
        ag, (dx, dx, dx), materials=materials, initial_temperature_K=initial_temperature,
        nonlinear_relative_tolerance=1.0e-8, linear_relative_tolerance=1.0e-8,
        maximum_nonlinear_iterations=48, maximum_linear_iterations=400,
        interface_resistance_m2_K_W=0.0,
    )
    dt_s = float(config["heat_duration_s_at_10um"]) / 100.0
    started = time.perf_counter()
    record, profile = advance_profiled(
        sim, dt_s, source, BoundarySet(), "ew_safe", fft_workers=4,
        preconditioner_method=preconditioner,
    )
    wall = time.perf_counter() - started
    iteration_rows = profile["newton_iterations"]
    return {
        "resolution": n,
        "preconditioner": preconditioner,
        "linear_policy": "ew_safe",
        "initial_temperature_K": initial_temperature,
        "dt_s": dt_s,
        "step_wall_seconds": wall,
        "amg_setup_seconds_sum": sum(float(x.get("amg_setup_seconds", 0.0)) for x in iteration_rows),
        "amg_levels": max((int(x["amg_levels"]) for x in iteration_rows if x.get("amg_levels") is not None), default=0),
        "newton_iterations": len(iteration_rows),
        "pcg_iterations": int(record["pcg_iterations_total"]),
        "final_nonlinear_relative_residual": float(record["nonlinear_relative_residual"]),
        "max_temperature_K": float(record["max_temperature_K"]),
        "max_ag_liquid_fraction": float(record["max_ag_liquid_fraction"]),
        "ag_bulk_liquid_fraction": float(record["ag_liquid_volume_fraction"]),
        "energy_residual_J": float(record["energy_residual_J"]),
        "normalized_energy_residual": abs(float(record["energy_residual_J"]))
        / max(abs(float(record["source_energy_step_J"])), abs(float(record["storage_change_step_J"])), 1e-30),
        "source_interpretation": "numerical localized volumetric pulse; not an arc boundary",
    }


def main():
    # 保持单进程内逐个运行，避免同时构建两个稀疏层级造成内存峰值。
    from numba import set_num_threads
    from run_benchmark import warm_numba

    set_num_threads(4)
    warmup_started = time.perf_counter()
    warm_numba()
    jit_warmup_seconds = time.perf_counter() - warmup_started
    rows = [run_one("fft"), run_one("amg")]
    reference = rows[0]
    accelerated = rows[1]
    accelerated["max_temperature_difference_vs_fft_K"] = abs(
        accelerated["max_temperature_K"] - reference["max_temperature_K"]
    )
    accelerated["bulk_liquid_fraction_relative_difference_vs_fft"] = abs(
        accelerated["ag_bulk_liquid_fraction"] - reference["ag_bulk_liquid_fraction"]
    ) / max(abs(reference["ag_bulk_liquid_fraction"]), 1e-30)
    if accelerated["max_temperature_difference_vs_fft_K"] > 1e-6:
        raise RuntimeError("AMG and FFT probe temperatures diverged")
    if accelerated["bulk_liquid_fraction_relative_difference_vs_fft"] > 1e-6:
        raise RuntimeError("AMG and FFT probe liquid fractions diverged")
    rows[0]["jit_warmup_seconds_excluded"] = jit_warmup_seconds
    rows[1]["jit_warmup_seconds_excluded"] = jit_warmup_seconds
    out = HERE / "results" / "AMG_PROBE.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rows, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    for row in rows:
        print(json.dumps(row, ensure_ascii=False, allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
