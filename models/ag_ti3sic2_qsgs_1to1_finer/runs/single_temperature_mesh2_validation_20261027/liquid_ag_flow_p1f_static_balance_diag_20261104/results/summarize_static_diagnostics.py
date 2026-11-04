#!/usr/bin/env python3
"""Summarize static-drop force-balance and grid-update diagnostics."""
import csv
import os
from pathlib import Path

ROOT = Path(os.environ.get("STATIC_DIAG_RESULTS_ROOT", Path(__file__).resolve().parent))
CASES = [
    "r16_amr_dt025_zero",
    "r16_fixed_dt025_zero",
    "r16_amr_frozen_dt025_zero",
    "r16_amr_dynamic_dt0125_zero",
    "r32_amr_dynamic_dt025_zero",
    "r32_amr_frozen_dt025_zero",
    "r16_amr_dynamic_dt025_plus",
    "r16_amr_dynamic_dt025_minus",
    "r16_amr_dynamic_tol1e-5",
    "r32_amr_dynamic_tol1e-5",
    "r16_amr_interface_only_dt025_zero",
]
FIELDS = [
    "case", "rows", "t_end_t0", "exit_code", "pressure_tolerance", "max_liquid_u", "time_max_liquid_u_t0",
    "max_liquid_u_rms", "time_max_liquid_u_rms_t0", "comparison_window_t0",
    "max_liquid_u_in_window", "max_liquid_u_rms_in_window", "max_mass_drift_pct",
    "pressure_jump_mean_in_window", "pressure_jump_min_in_window", "pressure_jump_max_in_window",
    "meridional_curvature_mean_in_window", "curvature_error_from_1_max_in_window_pct",
    "max_mgp_tolerance_ratio", "max_mgpf_tolerance_ratio", "samples_u_gt_0p01",
]

def maximum(rows, field):
    row = max(rows, key=lambda r: float(r[field]))
    return float(row[field]), float(row["t_over_t0"])

out = []
for name in CASES:
    folder = ROOT/name
    with (folder/"xcm_history.csv").open(newline="") as f:
        rows = list(csv.DictReader(f))
    code = int((folder/"exit_code.txt").read_text().strip())
    tolerance = 1e-5 if "tol1e-5" in name else 1e-4
    lo, hi = (0.065, 0.080) if name.startswith("r32_") else (0.173, 0.176)
    near = [r for r in rows if lo <= float(r["t_over_t0"]) <= hi]
    umax, ut = maximum(rows, "liquid_u_max")
    urms, urmst = maximum(rows, "liquid_u_rms")
    mgp_norm = [float(r["mgp_resa"])*float(r["dt"])**2/tolerance for r in rows]
    mgpf_norm = [float(r["mgpf_resa"])*(float(r["dt"])/2)**2/tolerance for r in rows]
    kapmean = sum(float(r["meridional_curvature_mean"]) for r in near)/len(near)
    out.append({
        "case": name,
        "rows": len(rows),
        "t_end_t0": rows[-1]["t_over_t0"],
        "exit_code": code,
        "pressure_tolerance": tolerance,
        "max_liquid_u": umax,
        "time_max_liquid_u_t0": ut,
        "max_liquid_u_rms": urms,
        "time_max_liquid_u_rms_t0": urmst,
        "comparison_window_t0": f"{lo:.3f}-{hi:.3f}",
        "max_liquid_u_in_window": max(float(r["liquid_u_max"]) for r in near),
        "max_liquid_u_rms_in_window": max(float(r["liquid_u_rms"]) for r in near),
        "max_mass_drift_pct": 100*max(abs(float(r["mass_drift_frac"])) for r in rows),
        "pressure_jump_mean_in_window": sum(float(r["pressure_jump"]) for r in near)/len(near),
        "pressure_jump_min_in_window": min(float(r["pressure_jump"]) for r in near),
        "pressure_jump_max_in_window": max(float(r["pressure_jump"]) for r in near),
        "meridional_curvature_mean_in_window": kapmean,
        "curvature_error_from_1_max_in_window_pct": 100*max(abs(float(r["meridional_curvature_mean"])-1.) for r in near),
        "max_mgp_tolerance_ratio": max(mgp_norm),
        "max_mgpf_tolerance_ratio": max(mgpf_norm),
        "samples_u_gt_0p01": sum(float(r["liquid_u_max"]) > .01 for r in rows),
    })

with (ROOT/"static_balance_summary.csv").open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(out)
for r in out:
    print(f"{r['case']}: exit={r['exit_code']} t/t0={float(r['t_end_t0']):.4f} "
          f"U_liq,max={float(r['max_liquid_u']):.5g} at {float(r['time_max_liquid_u_t0']):.4f} "
          f"U_liq,rms,max={float(r['max_liquid_u_rms']):.5g}; "
          f"Umax in {r['comparison_window_t0']}t0={float(r['max_liquid_u_in_window']):.5g}")
