#!/usr/bin/env python3
"""Collect final logged diagnostics from the six P1H height/boundary smoke runs."""
import csv
import pathlib
import sys


def read_rows(path):
    if not path.exists():
        return []
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def read_config(path):
    if not path.exists():
        return {}
    result = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            result[key] = value
    return result


def number(row, key, default=""):
    value = row.get(key, "")
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def main(root):
    records = []
    for nz in (64, 76, 86):
        for boundary in ("closed", "open"):
            folder = root / f"nz{nz}_{boundary}"
            config = read_config(folder / "run_config.env")
            snapshot_time = float(config.get("max_time_ps", "0.16"))
            history = read_rows(folder / "p1h_history.csv")
            pressure = read_rows(folder / "p1h_arc_pressure_history.csv")
            if not history:
                continue
            last = history[-1]
            last_pressure = pressure[-1] if pressure else {}
            records.append({
                "grid": f"64x64x{nz}",
                "gas_cap_cells": nz - 53,
                "gas_cap_nm": (nz - 53) * 19.31887901542944 / 64,
                "top_boundary": boundary,
                "last_logged_time_ps": number(last, "time_ps"),
                "snapshot_time_ps": snapshot_time,
                "tmax_K": number(last, "Tmax_K"),
                "tmax_Ag_K": number(last, "Tmax_Ag_K"),
                "tmax_skeleton_K": number(last, "Tmax_skeleton_K"),
                "liquid_Ag_volume_nm3": number(last, "liquid_Ag_volume_nm3"),
                "Ag_volume_drift_pct_last": number(last, "Ag_volume_drift_pct"),
                "Ag_volume_drift_pct_max": max(number(r, "Ag_volume_drift_pct", 0.) for r in history),
                "mean_surface_depression_nm": number(last, "mean_surface_depression_nm"),
                "source_center_depression_nm": number(last, "source_center_depression_nm"),
                "max_velocity_m_s": number(last, "max_velocity_m_s"),
                "energy_residual_J_last": number(last, "energy_residual_J"),
                "max_abs_energy_residual_J": max(abs(number(r, "energy_residual_J", 0.)) for r in history),
                "max_effective_interface_pressure_Pa": number(last_pressure, "max_effective_interface_Pa"),
                "integrated_pressure_force_z_N": number(last_pressure, "integrated_pressure_force_z_N"),
                "max_time_snapshot_present": int((folder / f"p1h_marangoni_off_max_time_t{snapshot_time:010.6f}_ps.bin.gz").exists()),
                "final_restart_present": int(any(folder.glob(f"p1h_restart_i*_t{snapshot_time:010.6f}_ps.dump.gz"))),
                "run_exit_reason": "MAX_TIME_STOP" if "MAX_TIME_STOP" in (folder / "run.stderr.log").read_text(errors="replace") else "INCOMPLETE",
            })
    out = root / "rectangular_height_sensitivity_summary.csv"
    if records:
        with out.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(records[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(records)
    print(f"wrote {out} ({len(records)} cases)")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: summarize_p1h_rectangular_height_sensitivity.py OUTPUT_ROOT")
    main(pathlib.Path(sys.argv[1]))
