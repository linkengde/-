#!/usr/bin/env python3
"""Run the dimensionless Basilisk P0 static-drop verification matrix."""
from __future__ import annotations

import csv
import json
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "basilisk_results"
CASES = HERE / "basilisk_cases"
SOURCE = HERE / "basilisk_static_drop.c"
DEFAULT_SRC = Path("/tmp/basilisk-build/src")
BASILISK_SRC = Path(os.environ.get("BASILISK_SRC", str(DEFAULT_SRC))).resolve()
QCC = BASILISK_SRC / "qcc"
LEVEL_BY_RESOLUTION = {16: 7, 32: 8, 64: 9}
T_END = float(os.environ.get("P0_T_END", "0.05"))
TIMEOUT_SECONDS = int(os.environ.get("P0_CASE_TIMEOUT", "1800"))


def compile_solver() -> None:
    env = os.environ.copy()
    env["BASILISK"] = str(BASILISK_SRC)
    env["PATH"] = f"{BASILISK_SRC}:{env.get('PATH', '')}"
    proc = subprocess.run(
        [str(QCC), "-grid=octree", "-O2", SOURCE.name, "-o", "basilisk_static_drop", "-lm"],
        cwd=HERE,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    (RESULTS / "compile.log").write_text(proc.stdout, encoding="utf-8")
    if proc.returncode:
        raise RuntimeError(f"Basilisk solver compilation failed; see {RESULTS / 'compile.log'}")


def parse_result(stderr: str) -> dict[str, str]:
    matches = [line for line in stderr.splitlines() if line.startswith("RESULT ")]
    if len(matches) != 1:
        raise RuntimeError(f"expected one RESULT line, got {len(matches)}: {stderr[-1000:]}")
    row: dict[str, str] = {}
    for token in shlex.split(matches[0])[1:]:
        key, value = token.split("=", 1)
        row[key] = value
    return row


def write_csv(path: Path, rows: list[dict[str, str]], fields: list[str] | None = None) -> None:
    if fields is None:
        fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def run_case(name: str, resolution: int, kind: str, curvature: str, dt_factor: float,
             end_time: float) -> dict[str, str]:
    level = LEVEL_BY_RESOLUTION[resolution]
    out_path = CASES / f"{name}.stdout.log"
    err_path = CASES / f"{name}.stderr.log"
    cmd = ["./basilisk_static_drop", str(level), kind, curvature, str(dt_factor), str(end_time)]
    proc = subprocess.run(
        cmd,
        cwd=HERE,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=TIMEOUT_SECONDS,
        check=False,
    )
    out_path.write_text(proc.stdout, encoding="utf-8")
    err_path.write_text(proc.stderr, encoding="utf-8")
    if proc.returncode:
        raise RuntimeError(f"{name} failed ({proc.returncode}); see {err_path}")
    row = parse_result(proc.stderr)
    row["case"] = name
    row["dt_factor"] = str(dt_factor)
    row["requested_end_time"] = str(end_time)
    row["run_mode"] = "geometry_only" if end_time == 0. else "dynamic"
    print(
        f"{name} [{row['run_mode']}]: t={row['t']} p_err={row['pressure_error_pct']}% "
        f"Umax={row['vmax']} theta={row['theta_deg']} mass={row['mass_drift_pct']}%",
        flush=True,
    )
    return row


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    CASES.mkdir(parents=True, exist_ok=True)
    if not QCC.is_file():
        raise RuntimeError(f"qcc not found at {QCC}; set BASILISK_SRC")
    compile_solver()

    plan: list[tuple[str, int, str, str, float, float]] = []
    for resolution in (16, 32):
        plan.extend([
            (f"sphere_exact_R{resolution}_dt0p5", resolution, "sphere", "exact", 0.5, T_END),
            (f"sphere_vof_R{resolution}_dt0p5", resolution, "sphere", "auto", 0.5, T_END),
            (f"sessile72_vof_R{resolution}_dt0p5", resolution, "sessile", "auto", 0.5, T_END),
        ])
    plan.extend([
        ("sphere_vof_R32_dt0p25", 32, "sphere", "auto", 0.25, T_END),
        ("sessile72_vof_R32_dt0p25", 32, "sessile", "auto", 0.25, T_END),
        ("sphere_vof_R64_geometry", 64, "sphere", "auto", 0.5, 0.),
        ("sessile72_vof_R64_geometry", 64, "sessile", "auto", 0.5, 0.),
    ])

    rows: list[dict[str, str]] = []
    for name, resolution, kind, curvature, dt_factor, end_time in plan:
        rows.append(run_case(name, resolution, kind, curvature, dt_factor, end_time))

    fields = [
        "case", "run_mode", "kind", "curvature", "level", "R_over_dx", "t", "cells", "steps", "dt_factor", "dt",
        "vmax", "Ca", "volume", "initial_volume", "initial_volume_error_pct", "mass_drift_pct",
        "p_jump", "p_expected_nominal", "pressure_error_pct", "p_expected_volume",
        "pressure_error_volume_pct", "theta_deg", "theta_error_deg", "fitted_radius",
        "kappa_mean", "kappa_std", "kappa_min", "kappa_max", "kappa_error_pct", "requested_end_time",
    ]
    write_csv(RESULTS / "vof_static_drop_summary.csv", [r for r in rows if r["kind"] == "sphere"], fields)
    write_csv(RESULTS / "vof_contact_angle_summary.csv", [r for r in rows if r["kind"] == "sessile72"], fields)
    write_csv(RESULTS / "vof_mass_conservation.csv", rows,
              ["case", "kind", "curvature", "R_over_dx", "t", "volume", "initial_volume",
               "initial_volume_error_pct", "mass_drift_pct"])

    conv = []
    for r in rows:
        if r["dt_factor"] != "0.5":
            continue
        conv.append({k: r.get(k, "") for k in (
            "case", "run_mode", "kind", "curvature", "R_over_dx", "pressure_error_pct", "vmax", "Ca",
            "mass_drift_pct", "theta_error_deg", "kappa_error_pct", "kappa_std", "dt", "steps")})
    write_csv(RESULTS / "vof_grid_convergence.csv", conv)

    sensitivity = []
    for family in ("sphere_vof", "sessile72_vof"):
        base = next(r for r in rows if r["case"] == f"{family}_R32_dt0p5")
        fine = next(r for r in rows if r["case"] == f"{family}_R32_dt0p25")
        for metric in ("pressure_error_pct", "vmax", "mass_drift_pct", "theta_error_deg", "kappa_error_pct"):
            sensitivity.append({
                "case_family": family,
                "metric": metric,
                "dt_factor_0p5": base.get(metric, ""),
                "dt_factor_0p25": fine.get(metric, ""),
            })
    write_csv(RESULTS / "vof_timestep_sensitivity.csv", sensitivity)
    (RESULTS / "basilisk_run_configuration.json").write_text(json.dumps({
        "solver": "Basilisk qcc, octree grid, centered incompressible Navier-Stokes + VOF + height-function curvature",
        "curvature_modes": {"exact": "prescribed kappa=2/R for isolated force/pressure balance", "auto": "VOF height-function curvature"},
        "dimensions": 3,
        "rho_liquid": 1.0,
        "rho_gas": 1.0,
        "dynamic_viscosity": 1.0,
        "surface_tension": 1.0,
        "radius": 1.0,
        "gravity": [0.0, 0.0, 0.0],
        "sessile_contact_angle_deg": 72.0,
        "contact_boundary": "back, corresponding to z=Z0; h.t and h.r",
        "dynamic_time_end": T_END,
        "geometry_only_time": 0.0,
        "time_step_factors": [0.5, 0.25],
        "dynamic_resolutions_R_over_dx": [16, 32],
        "geometry_only_resolutions_R_over_dx": [64],
        "time_step_sensitivity_resolution_R_over_dx": 32,
        "R64_dynamic_status": "not_completed; see basilisk_cases/diagnostics/R64_DYNAMIC_ATTEMPTS.md",
        "physical_material_properties_used": False,
    }, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (RuntimeError, subprocess.TimeoutExpired) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
