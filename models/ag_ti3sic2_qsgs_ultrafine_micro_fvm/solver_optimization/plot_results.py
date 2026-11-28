#!/usr/bin/env python3
"""Render summary plots from committed timing and threshold CSV files."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
FIGURES = HERE / "figures"
FIGURES.mkdir(exist_ok=True)


def read_rows(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def performance_plot():
    rows = read_rows(RESULTS / "OPTIMIZATION_COMPARISON.csv")
    n50_ids = [
        "bicontinuous_qsgs_seed1105_n50_baseline_fixed_dt0.25_t1",
        "bicontinuous_qsgs_seed1105_n50_numba_fixed_dt0.25_t1",
        "bicontinuous_qsgs_seed1105_n50_numba_fixed_dt0.25_t2",
        "bicontinuous_qsgs_seed1105_n50_numba_fixed_dt0.25_t4",
        "bicontinuous_qsgs_seed1105_n50_numba_ew_safe_dt0.25_t4",
    ]
    by_id = {row["run_id"]: row for row in rows}
    labels = ["Reference", "Numba 1T", "Numba 2T", "Numba 4T", "EW-safe 4T"]
    values = [float(by_id[key]["solver_step_wall_seconds"]) for key in n50_ids]
    fig, ax = plt.subplots(figsize=(8.6, 4.8), constrained_layout=True)
    bars = ax.bar(labels, values, color=["#76869b", "#4595c9", "#3483b5", "#24709e", "#d18f36"])
    ax.set_ylabel("Summed solver time per accepted step (s)")
    ax.set_title("50³ QGS thermal-cycle solver time")
    ax.grid(axis="y", alpha=0.25)
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, value, f"{value:.1f}", ha="center", va="bottom")
    fig.savefig(FIGURES / "performance_50.png", dpi=180)
    plt.close(fig)

    ids = ["qgs_n100_numba_ew_t4_full", "qgs_n100_numba_ew_safe_t4_full"]
    values = [float(by_id[key]["wall_seconds"]) / 60.0 for key in ids]
    retries = [int(by_id[key]["retries"]) for key in ids]
    fig, ax = plt.subplots(figsize=(6.8, 4.8), constrained_layout=True)
    bars = ax.bar(["EW", "EW-safe"], values, color=["#5688b4", "#d18f36"])
    ax.set_ylabel("Full-cycle wall time (min)")
    ax.set_title("100³ QGS transient cycle")
    ax.set_ylim(0.0, max(values) * 1.22)
    ax.grid(axis="y", alpha=0.25)
    for bar, value, retry in zip(bars, values, retries):
        ax.text(bar.get_x() + bar.get_width() / 2, value,
                f"{value:.1f} min\n{retry} retries", ha="center", va="bottom")
    fig.savefig(FIGURES / "performance_100.png", dpi=180)
    plt.close(fig)


def threshold_plot():
    path = RESULTS / "checkpoints" / "qgs_n100_numba_ew_safe_t4_full" / "liquid_threshold_sensitivity_heated.csv"
    rows = read_rows(path)
    by_ag = {}
    for row in rows:
        by_ag.setdefault(float(row["Ag_fraction_threshold"]), []).append(row)
    fig, ax = plt.subplots(figsize=(7.5, 4.8), constrained_layout=True)
    for ag_threshold, items in sorted(by_ag.items()):
        items.sort(key=lambda row: float(row["liquid_threshold"]))
        ax.plot([float(x["liquid_threshold"]) for x in items],
                [float(x["thresholded_geometric_volume_m3"]) / 1.0e-19 for x in items],
                marker="o", label=f"Ag threshold {ag_threshold:g}")
    ax.set_xlabel("Liquid-fraction threshold")
    ax.set_ylabel("Thresholded geometric volume (10⁻¹⁹ m³)")
    ax.set_title("100³ heated endpoint: liquid-volume threshold sensitivity")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(FIGURES / "liquid_threshold_sensitivity.png", dpi=180)
    plt.close(fig)


def main():
    performance_plot()
    threshold_plot()


if __name__ == "__main__":
    main()
