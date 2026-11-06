#!/usr/bin/env python3
"""Plot P1I temperature and liquid-pool histories from archived CSV files."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


def load_rows(path: Path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--history", required=True, type=Path)
    parser.add_argument("--clusters", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    history = load_rows(args.history)
    clusters = load_rows(args.clusters)
    time = [float(row["time_ps"]) for row in history]
    skeleton = [float(row["Tmax_skeleton_K"]) for row in history]
    silver = [float(row["Tmax_Ag_K"]) for row in history]
    liquid = [float(row["liquid_Ag_volume_nm3"]) for row in history]
    ctime = [float(row["time_ps"]) for row in clusters]
    largest = [float(row["largest_liquid_volume_nm3"]) for row in clusters]

    hold_start = 0.62344771
    heat_off = 0.800950919
    fig, (ax_t, ax_v) = plt.subplots(2, 1, figsize=(10.5, 7.5), sharex=True,
                                    constrained_layout=True)
    fig.suptitle("P1I: 3× arc-pressure / 3× Marangoni exploratory run", fontsize=14)
    for ax in (ax_t, ax_v):
        ax.axvspan(hold_start, heat_off, color="#f2aa45", alpha=0.18,
                   label="2100 K hold before safety cutoff")
        ax.axvline(heat_off, color="#a23b3b", linestyle="--", linewidth=1.2,
                   label="heat-off at 2300 K predictive guard")
        ax.grid(True, alpha=0.25)

    ax_t.plot(time, skeleton, color="#233b4a", linewidth=2, label="Ti₃SiC₂ max T")
    ax_t.plot(time, silver, color="#e57f22", linewidth=1.7, label="Ag max T")
    ax_t.axhline(2100, color="#e4a11b", linestyle=":", linewidth=1.4,
                 label="2100 K hold trigger")
    ax_t.axhline(2300, color="#a23b3b", linestyle="-.", linewidth=1.4,
                 label="2300 K property ceiling")
    ax_t.set_ylabel("Maximum temperature (K)")
    ax_t.legend(loc="upper right", fontsize=8, ncol=2)

    ax_v.plot(time, liquid, color="#de6b28", linewidth=2, label="total liquid Ag")
    ax_v.scatter(ctime, largest, s=17, color="#f2c900", edgecolor="#886f00",
                 linewidth=0.35, label="largest connected liquid cluster")
    ax_v.set_ylabel("Liquid Ag volume (nm³)")
    ax_v.set_xlabel("Simulation time (ps)")
    ax_v.legend(loc="upper right", fontsize=8)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=170, facecolor="white")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
