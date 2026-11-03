#!/usr/bin/env python3
"""Run and merge the R/Δx=64 VOF-curvature sphere and 72° sessile-drop cases."""
from __future__ import annotations

import csv
import json
import os
import shlex
import subprocess
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "basilisk_static_drop.c"
RESULTS = HERE / "basilisk_results"
CASES = HERE / "basilisk_cases"
BASILISK_SRC = Path(os.environ["BASILISK_SRC"]).resolve()
MPI_NP = int(os.environ.get("P0_MPI_NP", "4"))
MPIEXEC = shlex.split(os.environ.get("MPIEXEC", "mpirun"))
BINARY = Path(tempfile.gettempdir()) / "basilisk_static_drop_p0_mpi"


def run(command: list[str], *, env: dict[str, str], timeout: int = 1200,
        stdout_path: Path | None = None, stderr_path: Path | None = None) -> subprocess.CompletedProcess[str]:
    proc = subprocess.run(command, cwd=HERE, env=env, text=True, capture_output=True,
                          timeout=timeout, check=False)
    if stdout_path:
        stdout_path.write_text(proc.stdout, encoding="utf-8")
    if stderr_path:
        stderr_path.write_text(proc.stderr, encoding="utf-8")
    if proc.returncode:
        raise RuntimeError(f"command failed ({proc.returncode}): {' '.join(command)}")
    return proc


def parse_result(stderr: str) -> dict[str, str]:
    lines = [line for line in stderr.splitlines() if line.startswith("RESULT ")]
    if len(lines) != 1:
        raise RuntimeError(f"expected one RESULT line, got {len(lines)}")
    row: dict[str, str] = {}
    for item in shlex.split(lines[0])[1:]:
        key, value = item.split("=", 1)
        row[key] = value
    return row


def upsert_csv(path: Path, row: dict[str, str], fields: list[str] | None = None) -> None:
    if path.exists():
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            columns = reader.fieldnames or fields
            rows = list(reader)
    else:
        columns, rows = fields, []
    if not columns:
        raise RuntimeError(f"no CSV header defined for {path}")
    rows = [old for old in rows if old.get("case") != row.get("case")]
    rows.append({key: row.get(key, "") for key in columns})
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    if MPI_NP < 2:
        raise RuntimeError("P0_MPI_NP must be at least 2")
    RESULTS.mkdir(exist_ok=True)
    CASES.mkdir(exist_ok=True)
    env = os.environ.copy()
    env["BASILISK"] = str(BASILISK_SRC)
    env["PATH"] = f"{BASILISK_SRC}:{env.get('PATH', '')}"
    env.setdefault("CC99", "mpicc -std=c99 -pedantic -D_GNU_SOURCE=1 -Wno-unused-result")
    qcc = BASILISK_SRC / "qcc"
    compile_command = [str(qcc), "-grid=octree", "-D_MPI=1", "-O2",
                       SOURCE.name, "-o", str(BINARY), "-lm"]
    compile_proc = subprocess.run(compile_command, cwd=HERE, env=env, text=True,
                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                  timeout=300, check=False)
    (RESULTS / "compile_mpi.log").write_text(compile_proc.stdout, encoding="utf-8")
    if compile_proc.returncode:
        raise RuntimeError(f"MPI qcc build failed; see {RESULTS / 'compile_mpi.log'}")

    cases = [
        ("sphere_vof_R64_dt0p5", ["9", "sphere", "auto", "0.5", "0.05"]),
        ("sessile72_vof_R64_dt0p5", ["9", "sessile", "auto", "0.5", "0.05"]),
    ]
    sphere_fields = ["case", "run_mode", "kind", "curvature", "level", "R_over_dx", "t", "cells", "steps",
                     "dt_factor", "dt", "vmax", "Ca", "volume", "initial_volume", "initial_volume_error_pct",
                     "mass_drift_pct", "p_jump", "p_expected_nominal", "pressure_error_pct", "p_expected_volume",
                     "pressure_error_volume_pct", "theta_deg", "theta_error_deg", "fitted_radius", "kappa_mean",
                     "kappa_std", "kappa_min", "kappa_max", "kappa_error_pct", "requested_end_time"]
    for name, args in cases:
        out = CASES / f"{name}.stdout.log"
        err = CASES / f"{name}.stderr.log"
        proc = run([*MPIEXEC, "-np", str(MPI_NP), str(BINARY), *args], env=env,
                   stdout_path=out, stderr_path=err)
        row = parse_result(proc.stderr)
        row.update(case=name, run_mode="dynamic", dt_factor="0.5", requested_end_time="0.05")
        target = "vof_static_drop_summary.csv" if row["kind"] == "sphere" else "vof_contact_angle_summary.csv"
        upsert_csv(RESULTS / target, row, sphere_fields)
        mass_fields = ["case", "kind", "curvature", "R_over_dx", "t", "volume", "initial_volume",
                       "initial_volume_error_pct", "mass_drift_pct"]
        upsert_csv(RESULTS / "vof_mass_conservation.csv", row, mass_fields)
        grid_fields = ["case", "run_mode", "kind", "curvature", "R_over_dx", "pressure_error_pct", "vmax", "Ca",
                       "mass_drift_pct", "theta_error_deg", "kappa_error_pct", "kappa_std", "dt", "steps"]
        upsert_csv(RESULTS / "vof_grid_convergence.csv", row, grid_fields)
        print(f"{name}: pressure error {row['pressure_error_pct']}%, Umax={row['vmax']}, "
              f"mass drift={row['mass_drift_pct']}%, angle={row['theta_deg']}°", flush=True)

    config_path = RESULTS / "basilisk_run_configuration.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    config["dynamic_resolutions_R_over_dx"] = [16, 32, 64]
    config["dynamic_curvature_modes_by_resolution"] = {"16": ["exact", "vof"], "32": ["exact", "vof"], "64": ["vof"]}
    config["dynamic_contact_angle_resolutions_R_over_dx"] = [16, 32, 64]
    config["geometry_only_resolutions_R_over_dx"] = [64]
    config["R64_exact_curvature_status"] = "failed; pressure projection diverged; see diagnostics/R64_DYNAMIC_ATTEMPTS.md"
    config["R64_execution"] = {"mode": "MPI", "runtime": "Open MPI", "ranks": MPI_NP}
    config["R32_MPI_serial_consistency_check"] = "recorded under basilisk_cases/diagnostics/"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
