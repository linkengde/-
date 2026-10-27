#!/usr/bin/env python3
"""Plot matched Ag70 1x/2x peak and liquid-fraction history."""

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metrics", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    with args.metrics.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    colors = {"1": "#1769aa", "2": "#d55e00"}
    for refinement in ("1", "2"):
        selected = [row for row in rows if row["refinement"] == refinement]
        if not selected:
            continue
        time = [float(row["time_ps"]) for row in selected]
        peak = [float(row["peak_temperature_K"]) for row in selected]
        liquid = [float(row["ag_liquid_volume_fraction"] or 0.0) for row in selected]
        axes[0].plot(time, peak, label=f"{refinement}×", color=colors[refinement], linewidth=1.7)
        axes[1].plot(time, liquid, label=f"{refinement}×", color=colors[refinement], linewidth=1.7)
    axes[0].axhline(1235.08, color="black", linestyle="--", linewidth=1, label="Ag Tm candidate")
    axes[0].set(xlabel="time (ps)", ylabel="maximum temperature (K)", title="Ag70 Case A heating")
    axes[0].legend(frameon=False)
    axes[0].grid(alpha=.2)
    axes[1].set(xlabel="time (ps)", ylabel="Ag liquid volume fraction", title="No connected liquid observed if curve is zero")
    axes[1].set_ylim(bottom=0)
    axes[1].legend(frameon=False)
    axes[1].grid(alpha=.2)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)
    plt.close(fig)
    print(args.output)


if __name__ == "__main__":
    main()
