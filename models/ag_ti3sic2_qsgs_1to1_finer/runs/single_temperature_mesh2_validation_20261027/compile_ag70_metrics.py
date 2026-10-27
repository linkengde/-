#!/usr/bin/env python3
"""Compile native 1x/2x sampled heating records into the delivery CSV."""

from __future__ import annotations

import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = HERE / "ag70vol_melt_window"
FIELDS = [
    "refinement", "time_ps", "peak_temperature_K", "mean_ag_temperature_K",
    "ag_liquid_volume_fraction", "max_ag_liquid_fraction",
    "connected_liquid_cluster_count_at_threshold_0p5",
    "largest_cluster_voxels_at_threshold_0p5", "largest_cluster_volume_A3",
    "largest_cluster_track_id", "largest_cluster_track_duration_ps",
    "largest_cluster_centroid_A", "largest_cluster_depth_from_top_A",
    "largest_cluster_lateral_size_x_A", "largest_cluster_lateral_size_y_A",
    "ag_latent_energy_J", "total_enthalpy_J", "injected_energy_J",
    "energy_balance_error_J", "checkpoint_path",
]


def number(row: dict[str, str], key: str) -> str:
    return row.get(key, "") or ""


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def get_rows(refinement: int) -> list[dict[str, str]]:
    source = BASE / f"{refinement}x_0p2ps" / "heating.csv"
    run_dir = HERE / f"ag70vol_resume_{refinement}x"
    raw = read_csv(source)
    for path in sorted((run_dir / "segments").glob("heating_*.csv")):
        raw.extend(read_csv(path))
    checkpoint_times: dict[float, str] = {}
    base_checkpoint = BASE / f"{refinement}x_0p2ps" / "checkpoints" / "caseA_ag70_0p2ps_0p2ps.npz"
    if base_checkpoint.exists():
        checkpoint_times[0.2] = str(base_checkpoint.relative_to(HERE))
    for candidate in (run_dir / "checkpoints").glob(f"caseA_ag70_{refinement}x_heat_*.npz"):
        time_text = candidate.stem.rsplit("_", 1)[-1].removesuffix("ps").replace("p", ".")
        try:
            checkpoint_times[round(float(time_text), 9)] = str(candidate.relative_to(HERE))
        except ValueError:
            continue
    for candidate in (run_dir / "checkpoints").glob(f"caseA_ag70_{refinement}x_final_*.npz"):
        time_text = candidate.stem.rsplit("_", 1)[-1].removesuffix("ps").replace("p", ".")
        try:
            checkpoint_times[round(float(time_text), 9)] = str(candidate.relative_to(HERE))
        except ValueError:
            continue
    dedup: dict[float, dict[str, str]] = {}
    for row in raw:
        try:
            time_ps = float(row["time_ps"])
        except (KeyError, ValueError):
            continue
        dedup[round(time_ps, 9)] = row
    if not dedup:
        return []
    result: list[dict[str, str]] = []
    for time_ps, row in sorted(dedup.items()):
        checkpoint = checkpoint_times.get(round(time_ps, 9), "")
        result.append({
            "refinement": str(refinement),
            "time_ps": f"{time_ps:.9f}",
            "peak_temperature_K": number(row, "max_temperature_K"),
            "mean_ag_temperature_K": number(row, "mean_ag_temperature_K"),
            "ag_liquid_volume_fraction": number(row, "ag_liquid_volume_fraction"),
            "max_ag_liquid_fraction": number(row, "max_ag_liquid_fraction"),
            "connected_liquid_cluster_count_at_threshold_0p5": number(row, "connected_liquid_cluster_count"),
            "largest_cluster_voxels_at_threshold_0p5": number(row, "largest_liquid_cluster_voxels"),
            "largest_cluster_volume_A3": number(row, "largest_liquid_cluster_volume_A3"),
            "largest_cluster_track_id": number(row, "largest_liquid_track_id"),
            "largest_cluster_track_duration_ps": number(row, "largest_liquid_track_duration_ps"),
            "largest_cluster_centroid_A": number(row, "largest_liquid_cluster_centroid_A"),
            "largest_cluster_depth_from_top_A": number(row, "largest_liquid_cluster_depth_from_top_A"),
            "largest_cluster_lateral_size_x_A": number(row, "largest_liquid_cluster_lateral_size_x_A"),
            "largest_cluster_lateral_size_y_A": number(row, "largest_liquid_cluster_lateral_size_y_A"),
            "ag_latent_energy_J": number(row, "ag_latent_energy_J"),
            "total_enthalpy_J": number(row, "total_energy_J"),
            "injected_energy_J": number(row, "injected_energy_J"),
            "energy_balance_error_J": number(row, "energy_balance_error_J"),
            "checkpoint_path": checkpoint,
        })
    return result


def main() -> None:
    rows = get_rows(1) + get_rows(2)
    rows.sort(key=lambda row: (int(row["refinement"]), float(row["time_ps"])))
    out = HERE / "MELT_FREEZE_METRICS.csv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"output": str(out), "rows": len(rows), "rows_1x": len(get_rows(1)), "rows_2x": len(get_rows(2))}, indent=2))


if __name__ == "__main__":
    main()
