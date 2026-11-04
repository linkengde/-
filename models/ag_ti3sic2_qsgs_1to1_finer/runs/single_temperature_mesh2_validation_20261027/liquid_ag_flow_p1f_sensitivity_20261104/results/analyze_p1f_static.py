#!/usr/bin/env python3
"""Summarize zero-gradient static-drop pressure initialization diagnostics."""
import csv
import os
from pathlib import Path

HERE = Path(os.environ.get("P1F_RESULTS_ROOT", Path(__file__).resolve().parent))
CASES = [
    ("zero pressure initialization", "timestep/zero_gradient_dt_0p25x", "0"),
    ("positive Laplace jump", "timestep/zero_gradient_laplace_pressure", "+2 Gamma0/R"),
    ("negative Laplace jump", "timestep/zero_gradient_laplace_pressure_negative", "-2 Gamma0/R"),
]
LOW, HIGH = 0.30, 0.50

def read_rows(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))

def fit(points):
    ts = [p[0] for p in points]
    xs = [p[1] for p in points]
    mt, mx = sum(ts)/len(ts), sum(xs)/len(xs)
    den = sum((t-mt)**2 for t in ts)
    slope = sum((t-mt)*(x-mx) for t, x in points)/den
    ssr = sum((x-(mx+slope*(t-mt)))**2 for t, x in points)
    sst = sum((x-mx)**2 for x in xs)
    return slope, (1.-ssr/sst if sst else 1.)

out = []
for label, rel, pressure in CASES:
    folder = HERE/rel
    rows = read_rows(folder/"xcm_history.csv")
    code = int((folder/"exit_code.txt").read_text().strip())
    points = [(float(r["t_over_t0"]), float(r["xcm"])) for r in rows
              if LOW <= float(r["t_over_t0"]) <= HIGH]
    speed, r2 = fit(points)
    late = [r for r in rows if float(r["t_over_t0"]) >= LOW]
    peak = max(rows, key=lambda r: float(r["max_velocity"]))
    peak_late = max(late, key=lambda r: float(r["max_velocity"]))
    first, last = rows[0], rows[-1]
    out.append({
        "case": label,
        "initial_pressure_inside": pressure,
        "fit_window_t0": f"{LOW:.2f}-{HIGH:.2f}",
        "centroid_speed_over_Uref": f"{speed:.12g}",
        "fit_R2": f"{r2:.12g}",
        "centroid_shift_full_run": f"{float(last['xcm'])-float(first['xcm']):.12g}",
        "max_velocity_all": peak["max_velocity"],
        "time_of_max_velocity_t0": peak["t_over_t0"],
        "max_velocity_after_0p30t0": peak_late["max_velocity"],
        "time_of_late_max_velocity_t0": peak_late["t_over_t0"],
        "max_mass_drift_pct": f"{100.*max(abs(float(r['mass_drift_frac'])) for r in rows):.12g}",
        "max_mgp_resa_fit_window": f"{max(float(r['mgp_resa']) for r in rows if LOW <= float(r['t_over_t0']) <= HIGH):.12g}",
        "exit_code": code,
    })

path = HERE/"p1f_static_balance_summary.csv"
with path.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0]), lineterminator="\n")
    w.writeheader()
    w.writerows(out)
print(f"Static-drop diagnostics written to {path}")
