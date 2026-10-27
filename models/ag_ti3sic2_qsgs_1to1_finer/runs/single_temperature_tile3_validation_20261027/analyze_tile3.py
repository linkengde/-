"""Compare the original single cell with the central cell of the 3x3 run."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import label


RUN_DIR = Path(__file__).resolve().parent
PROJECT_DIR = RUN_DIR.parents[1]
BASE_RUN = PROJECT_DIR / "runs/single_temperature_restart_track_20261027"
GEOMETRY = PROJECT_DIR / "qsgs_morphology.npz"
HEAT_DIR = RUN_DIR / "heat"
SAMPLE_TIMES_PS = (0.2, 0.3, 0.4, 0.5, 0.6)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def nearest_row(rows: list[dict[str, str]], time_ps: float) -> dict[str, str]:
    row = min(rows, key=lambda item: abs(float(item["time_ps"]) - time_ps))
    if abs(float(row["time_ps"]) - time_ps) > 1.0e-8:
        raise ValueError(f"No sample near {time_ps} ps")
    return row


def nearest_cluster(rows: list[dict[str, str]], time_ps: float) -> dict[str, str]:
    candidates = [row for row in rows if row["is_main_cluster"].lower() == "true"]
    return nearest_row(candidates, time_ps)


def largest_component(mask: np.ndarray, spacing_A: np.ndarray) -> dict[str, object]:
    structure = np.zeros((3, 3, 3), dtype=bool)
    structure[1, 1, 1] = True
    for axis in range(3):
        low = [1, 1, 1]
        high = [1, 1, 1]
        low[axis] = 0
        high[axis] = 2
        structure[tuple(low)] = True
        structure[tuple(high)] = True
    labels, count = label(mask, structure=structure)
    if count == 0:
        return {"voxel_count": 0, "volume_A3": 0.0, "centroid_A": None,
                "depth_from_top_A": None, "lateral_size_x_A": 0.0, "lateral_size_y_A": 0.0}
    counts = np.bincount(labels.ravel(), minlength=count + 1)
    counts[0] = 0
    component_id = int(np.argmax(counts))
    xyz = np.where(labels == component_id)
    voxel_count = int(counts[component_id])
    centroid = [float((axis_indices.mean() + 0.5) * spacing_A[axis])
                for axis, axis_indices in enumerate(xyz)]
    return {
        "voxel_count": voxel_count,
        "volume_A3": float(voxel_count * np.prod(spacing_A)),
        "centroid_A": centroid,
        "depth_from_top_A": float(mask.shape[2] * spacing_A[2] - centroid[2]),
        "lateral_size_x_A": float((xyz[0].max() - xyz[0].min() + 1) * spacing_A[0]),
        "lateral_size_y_A": float((xyz[1].max() - xyz[1].min() + 1) * spacing_A[1]),
    }


def main() -> None:
    with np.load(GEOMETRY, allow_pickle=False) as archive:
        ag_mask_base = ~np.asarray(archive["skeleton"], dtype=bool)
        domain_size_A = np.asarray(archive["size_A"], dtype=float)
    geometry_sha256 = hashlib.sha256(GEOMETRY.read_bytes()).hexdigest()
    expected_geometry_sha256 = "7e0e1509b414597ce4e2a1c9fc387ed9656caf507f1186e2d28de79c3c1bfc43"
    if geometry_sha256 != expected_geometry_sha256:
        raise RuntimeError("QSGS geometry hash differs from the verified original input")
    ag_mask = np.tile(ag_mask_base, (3, 3, 1))
    spacing_A = domain_size_A / np.asarray(ag_mask_base.shape)
    source_center_local_A = np.asarray([449.4525, 452.7862269241275]) - domain_size_A[:2]
    xy = (np.arange(ag_mask_base.shape[0]) + 0.5) * spacing_A[0]
    yy = (np.arange(ag_mask_base.shape[1]) + 0.5) * spacing_A[1]
    lateral_radius_A = np.sqrt(
        (xy[:, None] - source_center_local_A[0]) ** 2
        + (yy[None, :] - source_center_local_A[1]) ** 2
    )
    central = (slice(ag_mask_base.shape[0], 2 * ag_mask_base.shape[0]),
               slice(ag_mask_base.shape[1], 2 * ag_mask_base.shape[1]), slice(None))

    base_rows = read_csv(BASE_RUN / "caseA_heat0p6.csv")
    tiled_rows = read_csv(HEAT_DIR / "caseA_tile3_heat0p6.csv")
    base_clusters = read_csv(BASE_RUN / "caseA_heat0p6.clusters.csv")
    tiled_clusters = read_csv(HEAT_DIR / "caseA_tile3_heat0p6.clusters.csv")

    comparisons: list[dict[str, object]] = []
    for time_ps in SAMPLE_TIMES_PS:
        base = nearest_row(base_rows, time_ps)
        tiled = nearest_row(tiled_rows, time_ps)
        base_main = nearest_cluster(base_clusters, time_ps)
        tiled_main = nearest_cluster(tiled_clusters, time_ps)
        checkpoint = HEAT_DIR / "checkpoints" / f"caseA_tile3_heat_{str(time_ps).replace('.', 'p')}ps.npz"
        with np.load(checkpoint, allow_pickle=False) as archive:
            temperature = np.asarray(archive["temperature_K"])
            liquid = np.asarray(archive["ag_liquid_fraction"])
            track_labels = np.asarray(archive["tracker_track_labels"])
            metadata = json.loads(str(archive["metadata_json"].item()))
        base_time_tag = f"{time_ps:.3f}".replace(".", "p")
        base_checkpoint = BASE_RUN / "checkpoints" / f"caseA_heat_{base_time_tag}ps.npz"
        with np.load(base_checkpoint, allow_pickle=False) as archive:
            base_temperature_field = np.asarray(archive["temperature_K"])
            base_liquid_field = np.asarray(archive["ag_liquid_fraction"])
        if abs(float(metadata["time_ps"]) - time_ps) > 1.0e-9:
            raise ValueError(f"Checkpoint time does not match requested sample: {checkpoint}")

        ag_central = ag_mask[central]
        liquid_central = liquid[central]
        temperature_central = temperature[central]
        temperature_difference = temperature_central - base_temperature_field
        liquid_difference = liquid_central - base_liquid_field
        abs_temperature_difference = np.where(ag_central, np.abs(temperature_difference), -np.inf)
        max_difference_index = tuple(int(value) for value in np.unravel_index(
            np.argmax(abs_temperature_difference), abs_temperature_difference.shape
        ))
        max_difference_coordinates_A = [
            float((max_difference_index[axis] + 0.5) * spacing_A[axis]) for axis in range(3)
        ]
        temperature_rise_norm = float(np.sum(ag_central * np.abs(base_temperature_field - 300.0)))
        liquid_l1_norm = float(np.sum(ag_central * base_liquid_field))
        central_threshold_mask = (liquid_central >= 0.5) & ag_central
        base_threshold_mask = (base_liquid_field >= 0.5) & ag_central
        threshold_union = int(np.count_nonzero(central_threshold_mask | base_threshold_mask))
        radius_50 = np.broadcast_to((lateral_radius_A <= 50.0)[:, :, None], ag_central.shape) & ag_central
        radius_100 = np.broadcast_to((lateral_radius_A <= 100.0)[:, :, None], ag_central.shape) & ag_central
        labels_central = track_labels[central]
        tile3_main_track_id = int(tiled_main["track_id"])
        central_main_mask = (labels_central == tile3_main_track_id) & (liquid_central >= 0.5) & ag_central
        central_liquid_mask = (liquid_central >= 0.5) & ag_central
        central_main = largest_component(central_main_mask, spacing_A)
        central_largest = largest_component(central_liquid_mask, spacing_A)
        central_main_voxels = int(central_main["voxel_count"])
        global_main_voxels = int(tiled_main["voxel_count"])
        base_centroid_A = json.loads(base_main["centroid_A"])
        local_centroid_A = central_main["centroid_A"]
        centroid_delta_A = (
            None if local_centroid_A is None else
            [float(local_centroid_A[i] - base_centroid_A[i]) for i in range(3)]
        )
        comparisons.append({
            "time_ps": time_ps,
            "base_cell": {
                "peak_temperature_K": float(base["max_temperature_K"]),
                "Ag_liquid_volume_fraction": float(base["ag_liquid_volume_fraction"]),
                "whole_cell_main_cluster_voxels": int(base_main["voxel_count"]),
                "whole_cell_main_cluster_volume_A3": float(base_main["volume_A3"]),
                "whole_cell_main_cluster_centroid_A": base_centroid_A,
                "whole_cell_main_cluster_depth_from_top_A": float(base_main["depth_from_top_A"]),
                "whole_cell_main_cluster_lateral_size_x_A": float(base_main["lateral_size_x_A"]),
                "whole_cell_main_cluster_lateral_size_y_A": float(base_main["lateral_size_y_A"]),
                "main_track_id_within_run": int(base_main["track_id"]),
                "main_track_continuous_duration_ps": float(base_main["continuous_duration_ps"]),
                "energy_balance_error_J": float(base["energy_balance_error_J"]),
            },
            "tile3_whole_domain": {
                "peak_temperature_K": float(tiled["max_temperature_K"]),
                "Ag_liquid_volume_fraction": float(tiled["ag_liquid_volume_fraction"]),
                "connected_liquid_cluster_count": int(tiled["connected_liquid_cluster_count"]),
                "liquid_voxel_count": int(tiled["liquid_voxel_count"]),
                "main_cluster_voxels": global_main_voxels,
                "main_cluster_volume_A3": float(tiled_main["volume_A3"]),
                "main_cluster_centroid_A": json.loads(tiled_main["centroid_A"]),
                "main_cluster_track_id_within_run": tile3_main_track_id,
                "main_track_continuous_duration_ps": float(tiled_main["continuous_duration_ps"]),
                "energy_balance_error_J": float(tiled["energy_balance_error_J"]),
                "injected_energy_J": float(tiled["injected_energy_J"]),
            },
            "tile3_central_cell": {
                "peak_Ag_temperature_K": float(np.max(temperature_central[ag_central])),
                "Ag_liquid_volume_fraction": float(np.sum(liquid_central * ag_central) / ag_central.sum()),
                "main_global_track_voxels_inside_cell": central_main_voxels,
                "main_global_track_volume_inside_cell_A3": central_main["volume_A3"],
                "main_global_track_centroid_local_A": central_main["centroid_A"],
                "main_global_track_depth_from_top_A": central_main["depth_from_top_A"],
                "main_global_track_lateral_size_x_A": central_main["lateral_size_x_A"],
                "main_global_track_lateral_size_y_A": central_main["lateral_size_y_A"],
                "global_main_track_fraction_inside_cell": (
                    central_main_voxels / global_main_voxels if global_main_voxels else 0.0
                ),
                "largest_local_cluster_voxels": central_largest["voxel_count"],
                "largest_local_cluster_volume_A3": central_largest["volume_A3"],
                "largest_local_cluster_centroid_local_A": central_largest["centroid_A"],
                "largest_local_cluster_depth_from_top_A": central_largest["depth_from_top_A"],
                "largest_local_cluster_lateral_size_x_A": central_largest["lateral_size_x_A"],
                "largest_local_cluster_lateral_size_y_A": central_largest["lateral_size_y_A"],
            },
            "central_vs_base_difference": {
                "peak_temperature_K": float(np.max(temperature_central[ag_central])) - float(base["max_temperature_K"]),
                "Ag_liquid_volume_fraction": float(np.sum(liquid_central * ag_central) / ag_central.sum()) - float(base["ag_liquid_volume_fraction"]),
                "main_track_cluster_voxels": central_main_voxels - int(base_main["voxel_count"]),
                "main_track_cluster_volume_A3": float(central_main["volume_A3"]) - float(base_main["volume_A3"]),
                "main_track_cluster_volume_relative_difference": (
                    float(central_main["volume_A3"]) / float(base_main["volume_A3"]) - 1.0
                    if float(base_main["volume_A3"]) else None
                ),
                "main_track_cluster_centroid_delta_A": centroid_delta_A,
                "main_track_cluster_depth_delta_A": (
                    None if central_main["depth_from_top_A"] is None else
                    float(central_main["depth_from_top_A"]) - float(base_main["depth_from_top_A"])
                ),
                "main_track_duration_delta_ps": float(tiled_main["continuous_duration_ps"]) - float(base_main["continuous_duration_ps"]),
                "central_field_max_abs_temperature_difference_K": float(np.max(abs_temperature_difference[ag_central])),
                "central_field_max_abs_temperature_difference_location_index_xyz": list(max_difference_index),
                "central_field_max_abs_temperature_difference_location_A": max_difference_coordinates_A,
                "central_field_max_abs_temperature_difference_lateral_radius_from_source_A": float(
                    lateral_radius_A[max_difference_index[0], max_difference_index[1]]
                ),
                "central_field_Ag_weighted_L1_temperature_difference_over_300K_rise": float(
                    np.sum(ag_central * np.abs(temperature_difference)) / max(temperature_rise_norm, 1e-300)
                ),
                "central_field_max_abs_temperature_difference_within_50A_K": float(
                    np.max(np.abs(temperature_difference[radius_50]))
                ),
                "central_field_max_abs_temperature_difference_within_100A_K": float(
                    np.max(np.abs(temperature_difference[radius_100]))
                ),
                "central_field_max_abs_liquid_fraction_difference": float(
                    np.max(np.abs(liquid_difference[ag_central]))
                ),
                "central_field_relative_L1_liquid_fraction_difference": float(
                    np.sum(ag_central * np.abs(liquid_difference)) / max(liquid_l1_norm, 1e-300)
                ),
                "central_field_f_l_ge_0p5_mask_IoU": (
                    float(np.count_nonzero(central_threshold_mask & base_threshold_mask) / threshold_union)
                    if threshold_union else 1.0
                ),
            },
        })

    result = {
        "status": "3X3_FINITE_SIZE_COMPARISON",
        "geometry_sha256": geometry_sha256,
        "threshold_ag_liquid_fraction": 0.5,
        "cluster_connectivity": "six face neighbours; periodic x/y for whole domain; central-cell crop reported without periodic seams",
        "identity_warning": "Track IDs are local to each run and are not compared as physical identities across runs.",
        "samples": comparisons,
    }
    (RUN_DIR / "tile3_comparison.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
