#!/usr/bin/env python3
"""Fit center-of-mass speed over one common physical-time window for P1F-1."""
import csv
import math
import os
from pathlib import Path

HERE = Path(os.environ.get("P1F_RESULTS_ROOT", Path(__file__).resolve().parent))
CASES = [
    ("DT-1", 1.0, "timestep/dt_1x", 1),
    ("DT-2", 0.5, "timestep/dt_0p5x", 1),
    ("DT-3", 0.25, "timestep/dt_0p25x", 1),
    ("zero-gradient", 0.25, "timestep/zero_gradient_dt_0p25x", 0),
]
FIT_LOW, FIT_HIGH = 0.30, 0.50
U_THEORY = -0.13333333333333333
ZERO_LIMIT = 0.01*abs(U_THEORY)
FIELDS = [
    "case", "dt_factor", "gradient", "PeT", "DT_cap", "fit_window_t0",
    "fit_speed_over_Uref", "fit_R2", "theory_speed_over_Uref",
    "speed_error_pct", "max_mass_drift_pct", "max_velocity",
    "mean_interface_T", "mean_interface_T_rms_departure",
    "mean_interface_gradT_tangent_rms", "sigma_interface_min",
    "sigma_interface_max", "mean_dt_in_fit_window", "max_dt_in_fit_window",
    "max_mgp_resa_in_fit_window", "finite_history", "exit_code",
]

def fit_speed(rows):
    points = [(float(r["t_over_t0"]), float(r["xcm"])) for r in rows
              if FIT_LOW <= float(r["t_over_t0"]) <= FIT_HIGH]
    if len(points) < 3:
        raise ValueError(f"need >=3 points in fit window, found {len(points)}")
    mt = sum(t for t, _ in points)/len(points)
    mx = sum(x for _, x in points)/len(points)
    denom = sum((t-mt)**2 for t, _ in points)
    slope = sum((t-mt)*(x-mx) for t, x in points)/denom
    ss_res = sum((x - (mx+slope*(t-mt)))**2 for t, x in points)
    ss_tot = sum((x-mx)**2 for _, x in points)
    r2 = 1. - ss_res/ss_tot if ss_tot > 0 else 1.
    # Uref*t0 = 1 in this benchmark, so dxcm/d(t/t0) equals U/Uref.
    return slope, r2, points

out_rows = []
for label, factor, rel, gradient in CASES:
    folder = HERE/rel
    code = int((folder/"exit_code.txt").read_text().strip())
    with (folder/"xcm_history.csv").open(newline="") as f:
        rows = list(csv.DictReader(f))
    finite_history = all(math.isfinite(float(v)) for r in rows for v in r.values())
    speed, r2, fit_points = fit_speed(rows)
    fit_rows = [r for r in rows if FIT_LOW <= float(r["t_over_t0"]) <= FIT_HIGH]
    max_drift_pct = 100.*max(abs(float(r["mass_drift_frac"])) for r in rows)
    dts = [float(r["dt"]) for r in fit_rows]
    sigma_min = min(float(r["sigma_interface_min"]) for r in fit_rows)
    sigma_max = max(float(r["sigma_interface_max"]) for r in fit_rows)
    row = {
        "case": label,
        "dt_factor": factor,
        "gradient": gradient,
        "PeT": 0.01,
        "DT_cap": f"{0.00506658613*factor:.12g}",
        "fit_window_t0": f"{FIT_LOW:.2f}-{FIT_HIGH:.2f}",
        "fit_speed_over_Uref": f"{speed:.12g}",
        "fit_R2": f"{r2:.12g}",
        "theory_speed_over_Uref": f"{U_THEORY:.12g}" if gradient else "",
        "speed_error_pct": f"{100.*abs(speed/U_THEORY-1.):.9g}" if gradient else "",
        "max_mass_drift_pct": f"{max_drift_pct:.9g}",
        "max_velocity": f"{max(float(r['max_velocity']) for r in rows):.12g}",
        "mean_interface_T": f"{sum(float(r['interface_T_mean']) for r in fit_rows)/len(fit_rows):.12g}",
        "mean_interface_T_rms_departure": f"{sum(float(r['interface_T_rms_departure']) for r in fit_rows)/len(fit_rows):.12g}",
        "mean_interface_gradT_tangent_rms": f"{sum(float(r['interface_gradT_tangent_rms']) for r in fit_rows)/len(fit_rows):.12g}",
        "sigma_interface_min": f"{sigma_min:.12g}",
        "sigma_interface_max": f"{sigma_max:.12g}",
        "mean_dt_in_fit_window": f"{sum(dts)/len(dts):.12g}",
        "max_dt_in_fit_window": f"{max(dts):.12g}",
        "max_mgp_resa_in_fit_window": f"{max(float(r['mgp_resa']) for r in fit_rows):.12g}",
        "finite_history": finite_history,
        "exit_code": code,
    }
    out_rows.append(row)

by_case = {r["case"]: r for r in out_rows}
half = float(by_case["DT-2"]["fit_speed_over_Uref"])
quarter = float(by_case["DT-3"]["fit_speed_over_Uref"])
delta_pct = 100.*abs(half-quarter)/abs(quarter)
speed_ok = delta_pct <= 1.
mass_ok = all(float(r["max_mass_drift_pct"]) < 0.1 for r in out_rows)
zero_ok = abs(float(by_case["zero-gradient"]["fit_speed_over_Uref"])) <= ZERO_LIMIT
exit_ok = all(int(r["exit_code"]) == 0 and r["finite_history"] for r in out_rows)
passed = speed_ok and mass_ok and zero_ok and exit_ok

with (HERE/"p1f_timestep_summary.csv").open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(out_rows)
(HERE/"p1f_timestep_decision.txt").write_text(
    f"{'PASS' if passed else 'FAIL'}\n"
    f"DT-2_vs_DT-3_speed_difference_pct={delta_pct:.9g}\n"
    f"criteria: difference<=1%, max mass drift<0.1%, zero-gradient |fit speed|<={ZERO_LIMIT:.9g} Uref, all exit codes 0\n"
    f"checks: speed={speed_ok} mass={mass_ok} zero_gradient={zero_ok} exit_codes={exit_ok}\n"
)
print(f"DT-2 vs DT-3 fit speed difference: {delta_pct:.6g}%")
print(f"P1F-1 decision: {'PASS' if passed else 'FAIL'}")
