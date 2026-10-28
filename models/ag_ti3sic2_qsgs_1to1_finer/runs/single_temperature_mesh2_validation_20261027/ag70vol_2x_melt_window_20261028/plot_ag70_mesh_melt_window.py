#!/usr/bin/env python3
"""Plot existing 1x data against all completed 2x continuation segments."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt


HERE = Path(__file__).resolve().parent
ONE_X = HERE.parent / "ag70vol_complete_20261028" / "MELT_FREEZE_METRICS.csv"
TWO_X_SEGMENTS = HERE / "heating_0p36_to_pool" / "segments"
OUTPUT = HERE / "ag70_1x_vs_2x_melt_window.png"
TM_K = 1235.08


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def series(rows: list[dict[str, str]], key: str) -> tuple[list[float], list[float]]:
    points: list[tuple[float, float]] = []
    for row in rows:
        if row.get("stage") != "heating":
            continue
        try:
            points.append((float(row["time_ps"]), float(row[key])))
        except (KeyError, TypeError, ValueError):
            continue
    points.sort()
    return [point[0] for point in points], [point[1] for point in points]


def main() -> None:
    one_rows = read_rows(ONE_X)
    two_rows: list[dict[str, str]] = []
    if TWO_X_SEGMENTS.exists():
        for path in sorted(TWO_X_SEGMENTS.glob("heating_*.csv")):
            segment_rows = read_rows(path)
            for row in segment_rows:
                row["stage"] = "heating"
            two_rows.extend(segment_rows)
    if not one_rows:
        raise SystemExit(f"Missing 1x reference metrics: {ONE_X}")

    fig, (temperature_ax, liquid_ax) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    for rows, label, color in ((one_rows, "1×", "tab:blue"), (two_rows, "2×", "tab:orange")):
        if not rows:
            continue
        time, temperature = series(rows, "peak_temperature_K")
        if time:
            temperature_ax.plot(time, temperature, label=label, color=color, linewidth=1.5)
        time, liquid = series(rows, "ag_liquid_volume_fraction")
        if time:
            liquid_ax.plot(time, liquid, label=label, color=color, linewidth=1.5)

    temperature_ax.axhline(TM_K, color="black", linestyle="--", linewidth=1, label="Ag melting reference")
    temperature_ax.set_ylabel("Peak temperature (K)")
    temperature_ax.grid(True, alpha=0.25)
    temperature_ax.legend(loc="best")
    liquid_ax.set_yscale("symlog", linthresh=1.0e-10)
    liquid_ax.set_ylabel("Ag liquid volume fraction")
    liquid_ax.set_xlabel("Simulation time (ps)")
    liquid_ax.grid(True, alpha=0.25)
    liquid_ax.legend(loc="best")
    fig.suptitle("Ag 70:30 vol% Case A heating: 1× reference and 2× continuation")
    fig.tight_layout()
    fig.savefig(OUTPUT, dpi=180)
    print(OUTPUT)
    print(f"1x heating rows: {sum(row.get('stage') == 'heating' for row in one_rows)}")
    print(f"2x heating rows: {sum(row.get('stage') == 'heating' for row in two_rows)}")


if __name__ == "__main__":
    main()
