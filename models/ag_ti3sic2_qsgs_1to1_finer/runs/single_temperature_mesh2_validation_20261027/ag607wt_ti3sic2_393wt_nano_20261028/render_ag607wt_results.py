#!/usr/bin/env python3
"""Render matched temperature fields and melt-pool histories from completed runs."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", "/tmp/ag607wt-matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/ag607wt-cache")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
GEOMETRY = HERE / "nano_geometry_ag607wt_ti3sic2_393wt.npz"


def read_case(refinement: int, output: Path) -> dict[str, Any]:
    summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
    milestones = summary.get("milestones") or {}
    relative_checkpoint = milestones.get("stable_pool_checkpoint") or summary["last_checkpoint"]
    checkpoint = output / relative_checkpoint
    with np.load(checkpoint, allow_pickle=False) as archive:
        temperature = np.asarray(archive["temperature_K"], dtype=np.float64).copy()
        liquid = np.asarray(archive["ag_liquid_fraction"], dtype=np.float64).copy()
        metadata = json.loads(str(archive["metadata_json"].item()))
    if refinement == 2:
        temperature = temperature.reshape(100, 2, 100, 2, 50, 2).mean(axis=(1, 3, 5))
        liquid = liquid.reshape(100, 2, 100, 2, 50, 2).mean(axis=(1, 3, 5))
    return {
        "refinement": refinement,
        "summary": summary,
        "checkpoint": checkpoint,
        "time_ps": float(metadata["time_ps"]),
        "temperature": temperature,
        "liquid": liquid,
        "tracker": metadata.get("tracker_metadata") or {},
    }


def read_history(output: Path) -> list[dict[str, Any]]:
    path = output / "segments.jsonl"
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line]


def render(output_1x: Path, output_2x: Path, destination: Path) -> tuple[Path, Path]:
    with np.load(GEOMETRY, allow_pickle=False) as archive:
        skeleton = np.asarray(archive["skeleton"], dtype=bool)
        size_A = np.asarray(archive["size_A"], dtype=float)
    spacing_A = size_A / np.asarray(skeleton.shape, dtype=float)
    cases = [read_case(1, output_1x), read_case(2, output_2x)]
    vmin = 300.0
    vmax = max(float(case["temperature"].max()) for case in cases)
    domain_xy = [0.0, float(size_A[0]), 0.0, float(size_A[1])]
    domain_xz = [0.0, float(size_A[0]), 0.0, float(size_A[2])]

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    image = None
    for col, case in enumerate(cases):
        temperature = case["temperature"]
        liquid = case["liquid"]
        active_z = np.sum(liquid, axis=(0, 1))
        if float(active_z.max()) > 0.0:
            z_index = int(np.argmax(active_z))
        else:
            z_index = int(np.unravel_index(np.argmax(temperature), temperature.shape)[2])
        tracker_records = case["tracker"].get("history", [])
        main = [
            record for record in tracker_records
            if record.get("is_main_cluster") and abs(float(record["time_ps"]) - case["time_ps"]) <= 1.0e-8
        ]
        if main:
            main_record = max(main, key=lambda row: int(row["voxel_count"]))
            y_index = int(float(main_record["centroid_A"][1]) / spacing_A[1]) % 100
        else:
            y_index = int(np.unravel_index(np.argmax(temperature), temperature.shape)[1])
        xy_ax = axes[0, col]
        image = xy_ax.imshow(
            temperature[:, :, z_index].T,
            origin="lower",
            extent=domain_xy,
            cmap="inferno",
            vmin=vmin,
            vmax=max(vmax, vmin + 1.0),
            aspect="equal",
        )
        if np.any(liquid[:, :, z_index] >= 0.5):
            xy_ax.contour(
                (np.arange(100) + 0.5) * spacing_A[0],
                (np.arange(100) + 0.5) * spacing_A[1],
                (liquid[:, :, z_index] >= 0.5).T,
                levels=[0.5],
                colors=["cyan"],
                linewidths=0.8,
            )
        xy_ax.set_title(f"{case['refinement']}× XY slice z={z_index * spacing_A[2]:.1f} Å | t={case['time_ps']:.3f} ps")
        xy_ax.set_xlabel("x (Å)")
        xy_ax.set_ylabel("y (Å)")

        xz_ax = axes[1, col]
        xz_ax.imshow(
            temperature[:, y_index, :].T,
            origin="lower",
            extent=domain_xz,
            cmap="inferno",
            vmin=vmin,
            vmax=max(vmax, vmin + 1.0),
            aspect="auto",
        )
        if np.any(liquid[:, y_index, :] >= 0.5):
            xz_ax.contour(
                (np.arange(100) + 0.5) * spacing_A[0],
                (np.arange(50) + 0.5) * spacing_A[2],
                (liquid[:, y_index, :] >= 0.5).T,
                levels=[0.5],
                colors=["cyan"],
                linewidths=0.8,
            )
        xz_ax.set_title(f"{case['refinement']}× XZ slice y={y_index * spacing_A[1]:.1f} Å")
        xz_ax.set_xlabel("x (Å)")
        xz_ax.set_ylabel("z (Å)")
    if image is not None:
        fig.colorbar(image, ax=axes.ravel().tolist(), label="Temperature (K)", shrink=0.82)
    fig.suptitle("Ag 60.7 wt% / Ti₃SiC₂ 39.3 wt% | 2× fields block-averaged to 1×")
    destination.mkdir(parents=True, exist_ok=True)
    fields_path = destination / "temperature_fields_and_liquid_pool.png"
    fig.savefig(fields_path, dpi=180)
    plt.close(fig)

    histories = [read_history(output_1x), read_history(output_2x)]
    summaries = [case["summary"] for case in cases]
    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True, constrained_layout=True)
    colors = ["#2667a5", "#d47f1f"]
    labels = ["1×", "2×"]
    metrics = [
        ("peak_temperature_K", "Peak temperature (K)"),
        ("ag_liquid_volume_fraction", "Ag liquid volume fraction"),
        ("largest_cluster_volume_A3", "Largest connected liquid cluster (Å³)"),
    ]
    for axis, (key, title) in zip(axes, metrics, strict=True):
        for idx, records in enumerate(histories):
            records = [record for record in records if key in record]
            if not records:
                continue
            times = [float(record["end_time_ps"]) for record in records]
            values = [float(record[key]) for record in records]
            axis.plot(times, values, color=colors[idx], label=labels[idx], linewidth=1.4)
            heat_off = summaries[idx].get("heat_off_time_ps")
            if heat_off is not None:
                axis.axvline(float(heat_off), color=colors[idx], linestyle="--", alpha=0.5)
        axis.set_ylabel(title)
        axis.grid(True, alpha=0.25)
        if axis.get_legend_handles_labels()[0]:
            axis.legend()
    axes[-1].set_xlabel("Simulation time (ps)")
    axes[0].set_title("Heating, stable-pool shutoff, and insulated-domain cooling histories")
    history_path = destination / "melt_solidification_histories.png"
    fig.savefig(history_path, dpi=180)
    plt.close(fig)
    return fields_path, history_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-1x", type=Path, default=HERE / "runs" / "full_1x")
    parser.add_argument("--output-2x", type=Path, default=HERE / "runs" / "full_2x")
    parser.add_argument("--destination", type=Path, default=HERE / "runs")
    args = parser.parse_args()
    outputs = render(args.output_1x.resolve(), args.output_2x.resolve(), args.destination.resolve())
    print(json.dumps({"status": "rendered", "outputs": [str(path) for path in outputs]}, indent=2))


if __name__ == "__main__":
    main()
