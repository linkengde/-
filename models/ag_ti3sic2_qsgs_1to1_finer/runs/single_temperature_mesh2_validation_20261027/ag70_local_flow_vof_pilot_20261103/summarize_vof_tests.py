#!/usr/bin/env python3
"""Summarize saved static-drop VOF validation runs and make a comparison plot."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
CASES = {
    "80": "binary_cell_center",
    "160": "binary_cell_center",
    "iso80": "binary_cell_center",
    "viscous80": "binary_cell_center",
    "isoSub8_80": "subcell_8x8_area_fraction",
    "foamSub8_80": "subcell_8x8_area_fraction",
    "foamSub8_160": "subcell_8x8_area_fraction",
    "foamSub8_nu1_80": "subcell_8x8_area_fraction",
    "foamSub8_nu1_t2_80": "subcell_8x8_area_fraction",
    "foamSub8_nu1_t10_80": "subcell_8x8_area_fraction",
}


def display_label(tag: str, data: dict) -> str:
    n = data["grid_shape_xy"][0]
    solver = "iso" if "interIsoFoam" in data["solver"] else "inter"
    init = "sub8" if tag.startswith(("isoSub8", "foamSub8")) else "binary"
    nu = data["kinematic_viscosity_both_phases_dimensionless"]
    end = data["final_time"]
    if nu > 0.1:
        return f"nu=1, t={end:g} {init}"
    return f"{solver} {n} {init}"


def main() -> None:
    rows = []
    labels = []
    for tag, initialization in CASES.items():
        data = json.loads((ROOT / f"static_drop_metrics_{tag}.json").read_text(encoding="utf-8"))
        rows.append({
            "tag": tag,
            "solver": data["solver"].replace("OpenFOAM ", "").replace(" 2406", ""),
            "grid_xy": "x".join(map(str, data["grid_shape_xy"])),
            "nu_dimensionless": data["kinematic_viscosity_both_phases_dimensionless"],
            "final_time": data["final_time"],
            "initialization": initialization,
            "relative_area_drift": data["relative_liquid_area_drift"],
            "expected_laplace_jump": data["expected_2d_laplace_pressure_jump"],
            "measured_bulk_pressure_jump": data["measured_bulk_pressure_jump"],
            "relative_pressure_jump_error": data["relative_pressure_jump_error"],
            "final_max_speed": data["final_max_speed"],
            "aspect_ratio": data["final_drop_covariance_aspect_ratio"],
            "continuity_error": data["cumulative_continuity_error_last"],
        })
        labels.append(display_label(tag, data))

    out = ROOT / "mesh_sensitivity.csv"
    with out.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    errors = np.asarray([100 * row["relative_pressure_jump_error"] for row in rows])
    speeds = np.asarray([row["final_max_speed"] for row in rows])
    x = np.arange(len(rows))
    fig, axes = plt.subplots(2, 1, figsize=(11, 7.5), sharex=True, constrained_layout=True)
    colors = ["#245b78" if row["nu_dimensionless"] <= 0.1 else "#a65d03" for row in rows]
    axes[0].bar(x, errors, color=colors)
    axes[0].axhline(0, color="#222222", linewidth=0.8)
    axes[0].set_ylabel("pressure-jump error (%)")
    axes[0].set_title("Static circular-drop VOF verification (dimensionless; not Ag flow)")
    axes[0].grid(axis="y", alpha=0.25)
    axes[1].bar(x, speeds, color=colors)
    axes[1].set_yscale("log")
    axes[1].set_ylabel("final maximum speed (dimensionless, log)")
    axes[1].set_xticks(x, labels, rotation=30, ha="right")
    axes[1].grid(axis="y", alpha=0.25)
    fig.savefig(ROOT / "static_drop_validation_comparison.png", dpi=180, facecolor="white")
    plt.close(fig)
    print(f"wrote {out.name} and static_drop_validation_comparison.png ({len(rows)} runs)")


if __name__ == "__main__":
    main()
