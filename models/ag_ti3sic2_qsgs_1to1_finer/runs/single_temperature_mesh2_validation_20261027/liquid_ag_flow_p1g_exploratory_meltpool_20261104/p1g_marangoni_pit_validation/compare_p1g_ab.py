#!/usr/bin/env python3
"""Create reproducible A/B tables and a compact P1G comparison plot."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_csv(path: Path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def nearest(rows, t, key="time_ps", tol=1.e-6):
    if not rows:
        return None
    row = min(rows, key=lambda r: abs(float(r[key]) - t))
    return row if abs(float(row[key]) - t) <= tol else None


def value(row, key):
    return "" if row is None else row.get(key, "")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("root", type=Path, help="P1G run root containing results/")
    p.add_argument("--output-dir", type=Path, default=None)
    a = p.parse_args()
    out = a.output_dir or (a.root / "results")
    out.mkdir(parents=True, exist_ok=True)
    case_dirs = {
        "A": a.root / "results/formal_case_A_marangoni_off_2100K_cut_cool_20261105",
        "B": a.root / "results/formal_case_B_marangoni_on_2100K_cut_cool_20261105",
    }
    hist = {k: read_csv(v / "p1g_history.csv") for k, v in case_dirs.items()}
    pits = {k: read_csv(v / "P1G_PIT_GEOMETRY.csv") for k, v in case_dirs.items()}

    columns = [
        "time_ps", "A_heat_on", "B_heat_on",
        "A_liquid_Ag_nm3", "B_liquid_Ag_nm3", "B_minus_A_liquid_Ag_nm3",
        "A_Tmax_Ag_K", "B_Tmax_Ag_K", "A_Tmax_skeleton_K", "B_Tmax_skeleton_K",
        "A_max_velocity_m_s", "B_max_velocity_m_s",
        "A_mean_surface_depression_nm", "B_mean_surface_depression_nm",
        "A_source_center_depression_nm", "B_source_center_depression_nm",
        "A_Ag_volume_drift_pct", "B_Ag_volume_drift_pct",
        "A_input_minus_bottom_mW", "B_input_minus_bottom_mW",
        "A_energy_residual_J", "B_energy_residual_J",
        "A_pit_depth_max_nm", "B_pit_depth_max_nm", "A_pit_area_nm2", "B_pit_area_nm2",
        "A_rim_volume_nm3", "B_rim_volume_nm3",
        "A_max_surface_rise_nm", "B_max_surface_rise_nm",
    ]
    out_csv = out / "P1G_A_B_COMPARISON.csv"
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        for ra in hist["A"]:
            t = float(ra["time_ps"])
            rb = nearest(hist["B"], t)
            pa, pb = nearest(pits["A"], t), nearest(pits["B"], t)
            la, lb = float(value(ra, "liquid_Ag_volume_nm3")), float(value(rb, "liquid_Ag_volume_nm3"))
            row = {
                "time_ps": f"{t:.9g}", "A_heat_on": value(ra, "heat_on_after_step"), "B_heat_on": value(rb, "heat_on_after_step"),
                "A_liquid_Ag_nm3": f"{la:.12g}", "B_liquid_Ag_nm3": f"{lb:.12g}", "B_minus_A_liquid_Ag_nm3": f"{lb-la:.12g}",
                "A_Tmax_Ag_K": value(ra, "Tmax_Ag_K"), "B_Tmax_Ag_K": value(rb, "Tmax_Ag_K"),
                "A_Tmax_skeleton_K": value(ra, "Tmax_skeleton_K"), "B_Tmax_skeleton_K": value(rb, "Tmax_skeleton_K"),
                "A_max_velocity_m_s": value(ra, "max_velocity_m_s"), "B_max_velocity_m_s": value(rb, "max_velocity_m_s"),
                "A_mean_surface_depression_nm": value(ra, "mean_surface_depression_nm"), "B_mean_surface_depression_nm": value(rb, "mean_surface_depression_nm"),
                "A_source_center_depression_nm": value(ra, "source_center_depression_nm"), "B_source_center_depression_nm": value(rb, "source_center_depression_nm"),
                "A_Ag_volume_drift_pct": value(ra, "Ag_volume_drift_pct"), "B_Ag_volume_drift_pct": value(rb, "Ag_volume_drift_pct"),
                "A_input_minus_bottom_mW": value(ra, "input_minus_bottom_mW"), "B_input_minus_bottom_mW": value(rb, "input_minus_bottom_mW"),
                "A_energy_residual_J": value(ra, "energy_residual_J"), "B_energy_residual_J": value(rb, "energy_residual_J"),
                "A_pit_depth_max_nm": value(pa, "pit_depth_max_nm"), "B_pit_depth_max_nm": value(pb, "pit_depth_max_nm"),
                "A_pit_area_nm2": value(pa, "pit_area_nm2"), "B_pit_area_nm2": value(pb, "pit_area_nm2"),
                "A_rim_volume_nm3": value(pa, "rim_volume_nm3"), "B_rim_volume_nm3": value(pb, "rim_volume_nm3"),
                "A_max_surface_rise_nm": value(pa, "max_surface_rise_nm"), "B_max_surface_rise_nm": value(pb, "max_surface_rise_nm"),
            }
            writer.writerow(row)

    fig, axs = plt.subplots(2, 2, figsize=(11, 7.2), dpi=180, sharex=True)
    colors = {"A": "#4169a1", "B": "#e28a24"}
    for k, label in [("A", "Marangoni off"), ("B", "Marangoni on")]:
        tt = np.array([float(r["time_ps"]) for r in hist[k]])
        axs[0, 0].plot(tt, [float(r["liquid_Ag_volume_nm3"]) for r in hist[k]], color=colors[k], label=label, lw=2)
        axs[0, 1].plot(tt, [float(r["max_velocity_m_s"]) for r in hist[k]], color=colors[k], label=label, lw=2)
        axs[1, 0].plot(tt, [-float(r["source_center_depression_nm"]) for r in hist[k]], color=colors[k], label=label, lw=2)
        axs[1, 1].plot(tt, [float(r["Tmax_skeleton_K"]) for r in hist[k]], color=colors[k], label=label, lw=2)
    axs[0, 0].set_ylabel("Liquid Ag (nm³)")
    axs[0, 1].set_ylabel("Maximum speed (m/s)")
    axs[1, 0].set_ylabel("Source-region surface elevation change (nm)\n(positive = uplift; negative = depression)")
    axs[1, 1].set_ylabel("Ti₃SiC₂ maximum temperature (K)")
    axs[1, 0].set_xlabel("Time (ps)")
    axs[1, 1].set_xlabel("Time (ps)")
    for ax in axs.flat:
        ax.axvline(0.624017255, color="0.4", ls="--", lw=0.9)
        ax.grid(alpha=0.22)
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle("P1G exploratory A/B run — matched 300 K starts; 64³ frozen grid\nDashed line: heat off at scaffold temperature threshold", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out_png = out / "P1G_A_B_COMPARISON.png"
    fig.savefig(out_png, bbox_inches="tight")
    print(f"wrote {out_csv} ({out_csv.stat().st_size} bytes)")
    print(f"wrote {out_png} ({out_png.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
