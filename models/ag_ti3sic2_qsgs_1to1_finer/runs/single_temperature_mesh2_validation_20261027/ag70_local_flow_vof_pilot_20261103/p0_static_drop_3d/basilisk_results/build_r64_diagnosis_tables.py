#!/usr/bin/env python3
"""Extract compact R64 AMR/pressure/force-balance tables from preserved logs."""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOGS = ROOT / "basilisk_cases" / "diagnostics"
OUT = ROOT / "basilisk_results"


def fields(line: str) -> dict[str, str]:
    out: dict[str, str] = {"record": line.split(maxsplit=1)[0]}
    for token in line.split()[1:]:
        if "=" in token:
            key, value = token.split("=", 1)
            out[key] = value
    return out


def read_records(name: str, kinds: tuple[str, ...]) -> list[dict[str, str]]:
    path = LOGS / name
    rows = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(kinds):
            rows.append(fields(line))
    return rows


def write_csv(name: str, rows: list[dict[str, str]], columns: list[str]) -> None:
    with (OUT / name).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


cases = [
    ("D2_freeze_geometric", "r64_exact_freeze.stderr.log", "exact", "initial geometric AMR frozen", "pass", "r64_exact_freeze.exit_code.txt"),
    ("D2_prefreeze_wavelet", "r64_exact_prefreeze.stderr.log", "exact", "one wavelet AMR pass during init, then frozen", "pressure diverged on first projection", "r64_exact_prefreeze.exit_code.txt"),
    ("D2_prefreeze_resetp", "r64_exact_prefreeze_resetp.stderr.log", "exact", "wavelet AMR during init; reset p=f-threshold Laplace field; then frozen", "same first-projection failure; pressure interpolation alone ruled out", "r64_exact_prefreeze_resetp.exit_code.txt"),
    ("D2_once_then_frozen", "r64_exact_once_then_frozen.stderr.log", "exact", "one wavelet AMR pass after first step, then frozen", "pressure diverged on second projection", "r64_exact_once_then_frozen.exit_code.txt"),
    ("D3_dynamic_exact", "r64_exact_dynamic.stderr.log", "exact", "wavelet AMR every step", "pressure diverged on second projection", "r64_exact_dynamic.exit_code.txt"),
    ("D3_dynamic_auto", "r64_auto_dynamic.stderr.log", "auto", "wavelet AMR every step", "two projections passed", "r64_auto_dynamic.exit_code.txt"),
]

fixed_rows: list[dict[str, str]] = [{
    "case_id": "D1_uniform_full_domain",
    "curvature": "exact",
    "mesh_policy": "uniform full 8R domain, target 512^3 cells",
    "initial_cells": "",
    "final_cells": "",
    "mesh_cell_sequence": "not initialized",
    "pressure_projections": "0",
    "first_projection_ratio": "",
    "second_projection_ratio": "",
    "max_velocity": "",
    "mass_drift_pct": "",
    "f_bounds": "",
    "nonfinite_cells": "",
    "status": "resource-limited; no solver result",
    "note": "30-GiB cgroup guard stopped the attempt before INIT; peak 32,243,068,928 bytes; cgroup OOM counters remained zero.",
}]

pressure_rows: list[dict[str, str]] = []
force_rows: list[dict[str, str]] = []

for case_id, log_name, curvature_mode, mesh_policy, conclusion, exit_name in cases:
    records = read_records(log_name, ("INIT ", "MESH ", "STEP ", "FORCE_BALANCE "))
    init = next((row for row in records if row["record"] == "INIT"), {})
    meshes = [row for row in records if row["record"] == "MESH"]
    steps = [row for row in records if row["record"] == "STEP"]
    forces = [row for row in records if row["record"] == "FORCE_BALANCE"]
    mesh_cells = [row.get("cells", "") for row in meshes if row.get("phase", "").startswith(("after", "before"))]
    initial_cells = init.get("cells", "")
    final_cells = steps[-1].get("cells", "") if steps else ""
    ratios = [row.get("pressure_residual_ratio", "") for row in steps]
    bad = [row for row in steps if int(row.get("pressure_iterations", "0") or 0) >= 100 or
           float(row.get("pressure_residual_ratio", "0") or 0) >= 1.]
    final = steps[-1] if steps else {}
    fixed_rows.append({
        "case_id": case_id,
        "curvature": curvature_mode,
        "mesh_policy": mesh_policy,
        "initial_cells": initial_cells,
        "final_cells": final_cells,
        "mesh_cell_sequence": "->".join([initial_cells, *mesh_cells]) if mesh_cells else initial_cells,
        "pressure_projections": str(len(steps)),
        "first_projection_ratio": ratios[0] if ratios else "",
        "second_projection_ratio": ratios[1] if len(ratios) > 1 else "",
        "max_velocity": max((float(row.get("vmax", "nan")) for row in steps), default=float("nan")),
        "mass_drift_pct": final.get("mass_drift_pct", ""),
        "f_bounds": f"[{final.get('fmin', '')},{final.get('fmax', '')}]" if final else "",
        "nonfinite_cells": final.get("nonfinite_cells", ""),
        "status": conclusion,
        "note": f"MPI-4; exit={(LOGS / exit_name).read_text().strip() if (LOGS / exit_name).exists() else 'recorded in log'}; see {log_name}.",
    })

    for step in steps:
        ratio = float(step.get("pressure_residual_ratio", "nan"))
        iters = int(step.get("pressure_iterations", "0"))
        pressure_rows.append({
            "case_id": case_id, "curvature": curvature_mode, "mesh_policy": mesh_policy,
            "step": step.get("i", ""), "projection_time": step.get("t", ""), "dt": step.get("dt", ""),
            "cells": step.get("cells", ""), "iterations": step.get("pressure_iterations", ""),
            "residual_before": step.get("pressure_resb", ""), "residual_after": step.get("pressure_resa", ""),
            "tolerance": step.get("pressure_tolerance", ""), "residual_to_tolerance": step.get("pressure_residual_ratio", ""),
            "rhs_sum": step.get("pressure_rhs_sum", ""), "pmin": step.get("pmin", ""), "pmax": step.get("pmax", ""),
            "vmax": step.get("vmax", ""), "urms_liquid": step.get("urms_liquid", ""),
            "fmin": step.get("fmin", ""), "fmax": step.get("fmax", ""),
            "mass_drift_pct": step.get("mass_drift_pct", ""), "nonfinite_cells": step.get("nonfinite_cells", ""),
            "projection_status": "diverged" if iters >= 100 or ratio >= 1. else "converged",
        })

    for force in forces:
        force_rows.append({
            "case_id": case_id, "curvature": curvature_mode, "mesh_policy": mesh_policy,
            "phase": force.get("phase", ""), "step": force.get("step", ""), "time": force.get("t", ""),
            "weighted_interface_area": force.get("weighted_interface_area", ""),
            "capillary_force_rms": force.get("capillary_rms", ""),
            "pressure_force_residual_rms": force.get("pressure_force_residual_rms", ""),
            "pressure_force_residual_max": force.get("pressure_force_residual_max", ""),
            "relative_residual_rms": force.get("relative_rms", ""),
        })

# Include the earlier two smaller-dt failure probes without rerunning them.
for factor in ("035", "020"):
    for step in read_records(f"exact_R64_factor{factor}.stderr.log", ("STEP ",)):
        iters = int(step.get("pressure_iterations", "0"))
        ratio = float(step.get("pressure_residual_ratio", "nan"))
        pressure_rows.append({
            "case_id": f"prior_exact_dt_factor_{factor}", "curvature": "exact", "mesh_policy": "dynamic AMR; prior timestep probe",
            "step": step.get("i", ""), "projection_time": step.get("t", ""), "dt": step.get("dt", ""),
            "cells": step.get("cells", ""), "iterations": step.get("pressure_iterations", ""),
            "residual_before": step.get("pressure_resb", ""), "residual_after": step.get("pressure_resa", ""),
            "tolerance": step.get("pressure_tolerance", ""), "residual_to_tolerance": step.get("pressure_residual_ratio", ""),
            "rhs_sum": step.get("pressure_rhs_sum", ""), "pmin": step.get("pmin", ""), "pmax": step.get("pmax", ""),
            "vmax": step.get("vmax", ""), "urms_liquid": step.get("urms_liquid", ""),
            "fmin": step.get("fmin", ""), "fmax": step.get("fmax", ""),
            "mass_drift_pct": step.get("mass_drift_pct", ""), "nonfinite_cells": step.get("nonfinite_cells", ""),
            "projection_status": "diverged" if iters >= 100 or ratio >= 1. else "converged",
        })

write_csv("r64_fixed_vs_amr.csv", fixed_rows, [
    "case_id", "curvature", "mesh_policy", "initial_cells", "final_cells", "mesh_cell_sequence",
    "pressure_projections", "first_projection_ratio", "second_projection_ratio", "max_velocity",
    "mass_drift_pct", "f_bounds", "nonfinite_cells", "status", "note",
])
write_csv("r64_pressure_projection.csv", pressure_rows, [
    "case_id", "curvature", "mesh_policy", "step", "projection_time", "dt", "cells", "iterations",
    "residual_before", "residual_after", "tolerance", "residual_to_tolerance", "rhs_sum", "pmin", "pmax",
    "vmax", "urms_liquid", "fmin", "fmax", "mass_drift_pct", "nonfinite_cells", "projection_status",
])
write_csv("r64_curvature_force_balance.csv", force_rows, [
    "case_id", "curvature", "mesh_policy", "phase", "step", "time", "weighted_interface_area",
    "capillary_force_rms", "pressure_force_residual_rms", "pressure_force_residual_max", "relative_residual_rms",
])

print(f"wrote fixed_vs_amr={len(fixed_rows)} pressure_projection={len(pressure_rows)} force_balance={len(force_rows)}")
