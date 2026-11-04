#!/usr/bin/env python3
"""Fit the frozen-AMR Marangoni control over the P1F comparison window."""
import csv
import os
from pathlib import Path

case = Path(os.environ.get(
    "FROZEN_CONTROL_CASE",
    Path(__file__).resolve().parent/"r16_amr_frozen_dt025_marangoni",
))
with (case/"xcm_history.csv").open(newline="") as f:
    rows = list(csv.DictReader(f))
points = [(float(r["t_over_t0"]), float(r["xcm"])) for r in rows
          if 0.30 <= float(r["t_over_t0"]) <= 0.50]
if len(points) < 3:
    raise SystemExit("frozen dynamic control lacks the complete 0.30–0.50 t0 fit interval")
ts, xs = [p[0] for p in points], [p[1] for p in points]
mt, mx = sum(ts)/len(ts), sum(xs)/len(xs)
speed = sum((t-mt)*(x-mx) for t,x in points)/sum((t-mt)**2 for t in ts)
ssr = sum((x-(mx+speed*(t-mt)))**2 for t,x in points)
sst = sum((x-mx)**2 for x in xs)
r2 = 1.-ssr/sst if sst else 1.
theory = -0.13333333333333333
out = {
    "case": case.name,
    "exit_code": int((case/"exit_code.txt").read_text().strip()),
    "fit_window_t0": "0.30-0.50",
    "fit_speed_over_Uref": speed,
    "theory_speed_over_Uref": theory,
    "speed_error_pct": 100.*abs(speed/theory-1.),
    "fit_R2": r2,
    "max_mass_drift_pct": 100.*max(abs(float(r["mass_drift_frac"])) for r in rows),
    "max_liquid_speed": max(float(r["liquid_u_max"]) for r in rows),
    "max_liquid_rms_speed": max(float(r["liquid_u_rms"]) for r in rows),
    "max_cells": max(int(r["cells"]) for r in rows),
    "t_end_t0": float(rows[-1]["t_over_t0"]),
}
path = case/"frozen_dynamic_summary.csv"
with path.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out), lineterminator="\n")
    w.writeheader()
    w.writerow(out)
print(out)
