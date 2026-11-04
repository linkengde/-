"""Fit common-window velocity and summarize P1F-2 thermal sensitivity."""
import csv
import os
from pathlib import Path

HERE = Path(os.environ.get("P1F_RESULTS_ROOT", Path(__file__).resolve().parent))
PES = [0.005, 0.01, 0.02, 0.05]
FIT_LOW, FIT_HIGH = 0.30, 0.50
FIELDS = [
    "PeT", "thermal_diffusivity", "dt_cap", "fit_window_t0",
    "fit_speed_over_Uref", "speed_change_vs_Pe0p01_pct", "fit_R2",
    "max_mass_drift_pct", "mean_interface_T", "mean_interface_T_rms_departure",
    "mean_interface_gradT_tangent_rms", "sigma_interface_min",
    "sigma_interface_max", "mean_dt_in_fit_window", "max_dt_in_fit_window",
    "max_velocity", "max_mgp_resa_in_fit_window", "exit_code",
]

def analyze(pe):
    name = f"PeT_{pe:g}".replace(".", "p")
    folder = HERE/"diffusivity"/name
    with (folder/"xcm_history.csv").open(newline="") as f:
        rows = list(csv.DictReader(f))
    code = int((folder/"exit_code.txt").read_text().strip())
    fit = [r for r in rows if FIT_LOW <= float(r["t_over_t0"]) <= FIT_HIGH]
    if len(fit) < 3:
        raise ValueError(f"{name}: not enough data in fixed fit window")
    ts = [float(r["t_over_t0"]) for r in fit]
    xs = [float(r["xcm"]) for r in fit]
    mt, mx = sum(ts)/len(ts), sum(xs)/len(xs)
    den = sum((t-mt)**2 for t in ts)
    speed = sum((t-mt)*(x-mx) for t,x in zip(ts,xs))/den
    ssr = sum((x-(mx+speed*(t-mt)))**2 for t,x in zip(ts,xs))
    sst = sum((x-mx)**2 for x in xs)
    r2 = 1.-ssr/sst if sst else 1.
    dts = [float(r["dt"]) for r in fit]
    return {
        "PeT": pe,
        "thermal_diffusivity": 0.066/pe,
        "dt_cap": 0.00506658613,
        "fit_window_t0": "0.30-0.50",
        "fit_speed_over_Uref": speed,
        "fit_speed_raw": speed,
        "fit_R2": r2,
        "max_mass_drift_pct": 100.*max(abs(float(r["mass_drift_frac"])) for r in rows),
        "mean_interface_T": sum(float(r["interface_T_mean"]) for r in fit)/len(fit),
        "mean_interface_T_rms_departure": sum(float(r["interface_T_rms_departure"]) for r in fit)/len(fit),
        "mean_interface_gradT_tangent_rms": sum(float(r["interface_gradT_tangent_rms"]) for r in fit)/len(fit),
        "sigma_interface_min": min(float(r["sigma_interface_min"]) for r in fit),
        "sigma_interface_max": max(float(r["sigma_interface_max"]) for r in fit),
        "mean_dt_in_fit_window": sum(dts)/len(dts),
        "max_dt_in_fit_window": max(dts),
        "max_velocity": max(float(r["max_velocity"]) for r in rows),
        "max_mgp_resa_in_fit_window": max(float(r["mgp_resa"]) for r in fit),
        "exit_code": code,
    }

rows = [analyze(pe) for pe in PES]
reference = next(r["fit_speed_raw"] for r in rows if r["PeT"] == 0.01)
for row in rows:
    row["speed_change_vs_Pe0p01_pct"] = 100.*(row["fit_speed_raw"]-reference)/abs(reference)
with (HERE/"p1f_diffusivity_summary.csv").open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\n")
    writer.writeheader()
    writer.writerows([{k:v for k,v in r.items() if k in FIELDS} for r in rows])
if any(r["exit_code"] != 0 or r["max_mass_drift_pct"] >= 0.1 for r in rows):
    raise SystemExit("P1F-2 has failed runs or mass drift >= 0.1%; inspect the logs")
print("P1F-2 summaries written; thermal sensitivity remains a dimensionless study.")
