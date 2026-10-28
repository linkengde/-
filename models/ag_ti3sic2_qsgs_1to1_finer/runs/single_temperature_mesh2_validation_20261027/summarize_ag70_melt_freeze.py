#!/usr/bin/env python3
"""Build a combined heating/cooling metrics table and evolution figure."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


FIELDS = [
    "stage",
    "time_ps",
    "peak_temperature_K",
    "mean_ag_temperature_K",
    "ag_liquid_volume_fraction",
    "max_ag_liquid_fraction",
    "connected_liquid_cluster_count",
    "liquid_voxel_count",
    "largest_cluster_voxels",
    "largest_cluster_volume_A3",
    "largest_cluster_centroid_A",
    "largest_cluster_depth_from_top_A",
    "largest_cluster_lateral_size_x_A",
    "largest_cluster_lateral_size_y_A",
    "largest_track_id",
    "largest_track_duration_ps",
    "total_energy_J",
    "injected_energy_J",
    "energy_balance_error_J",
    "relative_energy_balance_error",
]


def normalize(row: dict, stage: str) -> dict:
    return {
        "stage": stage,
        "time_ps": row.get("end_time_ps", row.get("time_ps", "")),
        "peak_temperature_K": row.get("peak_temperature_K", row.get("max_temperature_K", "")),
        "mean_ag_temperature_K": row.get("mean_ag_temperature_K", ""),
        "ag_liquid_volume_fraction": row.get("ag_liquid_volume_fraction", ""),
        "max_ag_liquid_fraction": row.get("max_ag_liquid_fraction", ""),
        "connected_liquid_cluster_count": row.get("connected_cluster_count", row.get("connected_liquid_cluster_count", "")),
        "liquid_voxel_count": row.get("liquid_voxel_count", ""),
        "largest_cluster_voxels": row.get("largest_cluster_voxels", row.get("largest_liquid_cluster_voxels", "")),
        "largest_cluster_volume_A3": row.get("largest_cluster_volume_A3", row.get("largest_liquid_cluster_volume_A3", "")),
        "largest_cluster_centroid_A": row.get("largest_cluster_centroid_A", row.get("largest_liquid_cluster_centroid_A", "")),
        "largest_cluster_depth_from_top_A": row.get("largest_cluster_depth_from_top_A", row.get("depth_from_top_A", "")),
        "largest_cluster_lateral_size_x_A": row.get("largest_cluster_lateral_size_x_A", row.get("lateral_size_x_A", "")),
        "largest_cluster_lateral_size_y_A": row.get("largest_cluster_lateral_size_y_A", row.get("lateral_size_y_A", "")),
        "largest_track_id": row.get("largest_track_id", row.get("largest_liquid_track_id", "")),
        "largest_track_duration_ps": row.get("largest_track_duration_ps", row.get("largest_liquid_track_duration_ps", row.get("largest_track_continuous_duration_ps", ""))),
        "total_energy_J": row.get("total_energy_J", row.get("total_enthalpy_J", "")),
        "injected_energy_J": row.get("injected_energy_J", ""),
        "energy_balance_error_J": row.get("energy_balance_error_J", ""),
        "relative_energy_balance_error": row.get("relative_energy_balance_error", ""),
    }


def read_jsonl(path: Path) -> list[dict]:
    with path.open() as handle:
        return [json.loads(line) for line in handle if line.strip()]


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def initial_checkpoint_row(checkpoint: Path, geometry: Path) -> dict:
    with np.load(checkpoint, allow_pickle=False) as data:
        temperature = np.asarray(data["temperature_K"])
        liquid = np.asarray(data["ag_liquid_fraction"])
        metadata = json.loads(str(data["metadata_json"].item()))
    with np.load(geometry, allow_pickle=False) as data:
        ag = ~np.asarray(data["skeleton"], dtype=bool)
    refinement = temperature.shape[0] // ag.shape[0]
    if temperature.shape != tuple(np.asarray(ag.shape) * refinement):
        raise ValueError("Initial checkpoint shape is incompatible with the QSGS geometry")
    if refinement > 1:
        ag = np.repeat(np.repeat(np.repeat(ag, refinement, axis=0), refinement, axis=1), refinement, axis=2)
    values = {
        "time_ps": metadata["time_ps"],
        "peak_temperature_K": float(temperature.max()),
        "mean_ag_temperature_K": float(temperature[ag].mean()),
        "ag_liquid_volume_fraction": float(np.sum(liquid[ag]) / ag.sum()),
        "max_ag_liquid_fraction": float(liquid[ag].max()),
        "connected_cluster_count": 0,
        "liquid_voxel_count": 0,
        "largest_cluster_voxels": 0,
        "largest_cluster_volume_A3": 0.0,
        "largest_liquid_cluster_centroid_A": "",
        "largest_liquid_track_id": "",
        "largest_liquid_track_duration_ps": 0.0,
        "total_energy_J": metadata["total_enthalpy_J"],
        "injected_energy_J": metadata["injected_energy_J"],
        "energy_balance_error_J": metadata["energy_balance_error_J"],
        "relative_energy_balance_error": metadata["relative_energy_balance_error"],
    }
    return normalize(values, "heating")


def enrich_from_tracks(rows: list[dict], track_path: Path) -> None:
    if not track_path.exists():
        return
    by_time: dict[float, list[dict]] = {}
    for track in read_csv(track_path):
        by_time.setdefault(round(float(track["time_ps"]), 9), []).append(track)
    for row in rows:
        candidates = by_time.get(round(float(row["time_ps"]), 9), [])
        if not candidates:
            continue
        track = max(candidates, key=lambda item: int(item.get("voxel_count") or 0))
        row["largest_cluster_voxels"] = track.get("voxel_count", row["largest_cluster_voxels"])
        row["largest_cluster_volume_A3"] = track.get("volume_A3", row["largest_cluster_volume_A3"])
        row["largest_cluster_centroid_A"] = track.get("centroid_A", "")
        row["largest_cluster_depth_from_top_A"] = track.get("depth_from_top_A", "")
        row["largest_cluster_lateral_size_x_A"] = track.get("lateral_size_x_A", "")
        row["largest_cluster_lateral_size_y_A"] = track.get("lateral_size_y_A", "")
        row["largest_track_id"] = track.get("track_id", "")
        row["largest_track_duration_ps"] = track.get("continuous_duration_ps", "")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("heating_dir", type=Path)
    parser.add_argument("cooling_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--initial-checkpoint", type=Path)
    parser.add_argument("--geometry", type=Path)
    args = parser.parse_args()

    heating = [normalize(row, "heating") for row in read_jsonl(args.heating_dir / "segments.jsonl")]
    enrich_from_tracks(heating, args.heating_dir / "heating.clusters.csv")
    if args.initial_checkpoint:
        if not args.geometry:
            parser.error("--geometry is required with --initial-checkpoint")
        heating.insert(0, initial_checkpoint_row(args.initial_checkpoint, args.geometry))
    cooling = [normalize(row, "cooling") for row in read_csv(args.cooling_dir / "cooling.csv")]
    enrich_from_tracks(cooling, args.cooling_dir / "cooling.clusters.csv")
    rows = sorted(heating + cooling, key=lambda row: float(row["time_ps"]))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    metrics_path = args.output_dir / "MELT_FREEZE_METRICS.csv"
    with metrics_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    heat_rows = sorted(heating, key=lambda row: float(row["time_ps"]))
    cool_rows = sorted(cooling, key=lambda row: float(row["time_ps"]))
    times = [float(row["time_ps"]) for row in rows]
    peaks = [float(row["peak_temperature_K"]) for row in rows]
    liquid = [float(row["ag_liquid_volume_fraction"] or 0.0) for row in rows]
    cluster_volume = [float(row["largest_cluster_volume_A3"] or 0.0) for row in rows]

    first_liquid = next(float(row["time_ps"]) for row in heating if float(row["ag_liquid_volume_fraction"] or 0.0) > 0)
    first_pool = next(float(row["time_ps"]) for row in heating if int(row["largest_cluster_voxels"] or 0) >= 2)
    heat_off = max(float(row["time_ps"]) for row in heating)
    all_solid = next(float(row["time_ps"]) for row in cooling if float(row["ag_liquid_volume_fraction"] or 0.0) == 0)
    fig, axes = plt.subplots(3, 2, figsize=(14, 10), sharex="col", constrained_layout=True)
    series = [
        (peaks, "Peak temperature (K)", "#b23a48"),
        (liquid, "Ag liquid volume fraction", "#1769aa"),
        (cluster_volume, "Largest liquid cluster (Å³)", "#00876c"),
    ]
    for row_index, (values, ylabel, color) in enumerate(series):
        for col_index, axis in enumerate(axes[row_index]):
            axis.plot(times, values, color=color, linewidth=1.6)
            axis.set_ylabel(ylabel)
            axis.grid(alpha=0.2)
            axis.axvline(first_liquid, color="#7b2cbf", linestyle=":", linewidth=1)
            axis.axvline(first_pool, color="#00876c", linestyle=":", linewidth=1)
            axis.axvline(heat_off, color="#e07a00", linestyle="--", linewidth=1)
            axis.axvline(all_solid, color="#555555", linestyle="-.", linewidth=1)
            if row_index == 0:
                axis.axhline(1235.08, color="#333333", linestyle="--", linewidth=1, label="Ag melting reference")
                if col_index == 0:
                    axis.legend(frameon=False, loc="upper left")
        axes[row_index, 0].set_xlim(min(times), max(times))
        axes[row_index, 1].set_xlim(4.08, max(times))
    axes[0, 1].text(first_liquid + 0.003, axes[0, 1].get_ylim()[0] + 0.94 * (axes[0, 1].get_ylim()[1] - axes[0, 1].get_ylim()[0]), "first liquid", va="top", fontsize=8, color="#7b2cbf")
    axes[0, 1].text(first_pool + 0.003, axes[0, 1].get_ylim()[0] + 0.84 * (axes[0, 1].get_ylim()[1] - axes[0, 1].get_ylim()[0]), "multi-voxel", va="top", fontsize=8, color="#00876c")
    axes[0, 1].text(heat_off + 0.003, axes[0, 1].get_ylim()[0] + 0.72 * (axes[0, 1].get_ylim()[1] - axes[0, 1].get_ylim()[0]), "source off", va="top", fontsize=8, color="#e07a00")
    axes[0, 1].text(all_solid + 0.003, axes[0, 1].get_ylim()[0] + 0.60 * (axes[0, 1].get_ylim()[1] - axes[0, 1].get_ylim()[0]), "all solid", va="top", fontsize=8, color="#555555")
    for axis in axes[-1]:
        axis.set_xlabel("Simulation time (ps)")
    axes[-1, 0].set_xticks([1, 2, 3, 4])
    axes[-1, 1].set_xticks([4.1, 4.2, 4.3, 4.4, 4.5])
    axes[0, 0].set_title("Full heating and cooling record")
    axes[0, 1].set_title("Melt-pool window (expanded)")
    figure_path = args.output_dir / "ag70_melt_freeze_evolution.png"
    fig.suptitle("Ag70 vol% Case A: heating and insulated source-off cooling")
    fig.savefig(figure_path, dpi=180)
    plt.close(fig)
    print(json.dumps({"metrics_csv": str(metrics_path), "evolution_figure": str(figure_path), "heating_rows": len(heating), "cooling_rows": len(cooling)}, indent=2))


if __name__ == "__main__":
    main()
