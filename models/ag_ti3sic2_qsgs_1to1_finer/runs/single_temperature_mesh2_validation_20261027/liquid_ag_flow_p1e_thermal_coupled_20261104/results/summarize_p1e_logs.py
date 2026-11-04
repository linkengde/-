#!/usr/bin/env python3
"""Extract final diagnostics and maximum logged mass drift from P1E runs."""
import csv
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASES = [
    ("R8_fixed_gradient", "thermal_R8_Pe001_fixed"),
    ("R8_AMR_gradient", "thermal_R8_Pe001_amr"),
    ("R8_AMR_zero_gradient", "thermal_R8_Pe001_zero_gradient"),
    ("R16_AMR_gradient", "thermal_R16_Pe001_amr_half_t0"),
]
LOG_COLUMNS = [
    "t_over_t0", "mass_drift_frac", "xcm", "mean_u_over_Uref",
    "interval_speed_over_Uref", "umax", "interface_T_rms", "dt",
    "cells", "mgp_i", "mgp_resa", "PeT", "gradT",
]
OUTPUT = HERE / "p1e_thermal_coupled_summary.csv"
FIELDS = [
    "case", "cells_per_radius", "PeT", "gradT", "t_over_t0",
    "interval_speed_over_Uref", "theory_speed_over_Uref", "speed_error_pct",
    "max_abs_mass_drift_pct", "final_abs_mass_drift_pct", "max_velocity",
    "interface_T_rms", "final_grid_cells", "exit_code",
]

rows = []
for case, stem in CASES:
    out_path = HERE / f"{stem}.stdout.log"
    err_path = HERE / f"{stem}.stderr.log"
    code_path = HERE / f"{stem}.exit_code.txt"
    stdout = out_path.read_text()
    stderr = err_path.read_text()
    report_line = next((line for line in stderr.splitlines() if line.startswith("REPORT ")), None)
    if report_line is None:
        raise SystemExit(f"missing REPORT line in {err_path}")
    report = dict(re.findall(r"([A-Za-z0-9_]+)=([^\s]+)", report_line))
    samples = []
    for line in stdout.splitlines():
        parts = line.split()
        if len(parts) != len(LOG_COLUMNS):
            continue
        try:
            samples.append(dict(zip(LOG_COLUMNS, map(float, parts))))
        except ValueError:
            continue
    if not samples:
        raise SystemExit(f"no numeric samples in {out_path}")
    drifts = [abs(sample["mass_drift_frac"]) for sample in samples]
    grad = float(report["gradT"])
    rows.append({
        "case": case,
        "cells_per_radius": int(report["cells_per_radius"]),
        "PeT": report["PeT"],
        "gradT": report["gradT"],
        "t_over_t0": report["t_over_t0"],
        "interval_speed_over_Uref": report["measured_speed_over_Uref"],
        "theory_speed_over_Uref": report["theory_over_Uref"] if grad else "",
        "speed_error_pct": report["rel_error_pct"] if grad else "",
        "max_abs_mass_drift_pct": f"{100*max(drifts):.9g}",
        "final_abs_mass_drift_pct": f"{100*drifts[-1]:.9g}",
        "max_velocity": report["max_u"],
        "interface_T_rms": report["interface_T_rms"],
        "final_grid_cells": int(report["cells"]),
        "exit_code": int(code_path.read_text().strip()),
    })

with OUTPUT.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
print(f"wrote {OUTPUT}")
