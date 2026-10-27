#!/usr/bin/env python3
"""Compare 1x and block-coarsened 2x Ag 70 vol% heating checkpoints."""

from __future__ import annotations

import json
import hashlib
from pathlib import Path

import numpy as np


RUN_DIR = Path(__file__).resolve().parent
PROJECT_DIR = RUN_DIR.parents[1]
PHASE_DIR = PROJECT_DIR / "phase_field_qsgs"
WINDOW_DIR = RUN_DIR / "ag70vol_melt_window"
GEOMETRY = RUN_DIR / "qsgs_morphology_ag70vol.npz"
CHECKPOINT_NAME = "caseA_ag70_0p2ps_0p2ps.npz"
SUMMARY_1X = WINDOW_DIR / "1x_0p2ps" / "summary.json"
SUMMARY_2X = WINDOW_DIR / "2x_0p2ps" / "summary.json"
CHECKPOINT_1X = WINDOW_DIR / "1x_0p2ps" / "checkpoints" / CHECKPOINT_NAME
CHECKPOINT_2X = WINDOW_DIR / "2x_0p2ps" / "checkpoints" / CHECKPOINT_NAME


def load_fields(path: Path) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    with np.load(path, allow_pickle=False) as archive:
        metadata = json.loads(str(archive["metadata_json"].item()))
        temperature = np.asarray(archive["temperature_K"], dtype=float)
        liquid = np.asarray(archive["ag_liquid_fraction"], dtype=float)
    return temperature, liquid, metadata


def peak_coordinate(field: np.ndarray, spacing_A: np.ndarray) -> list[float]:
    index = np.asarray(np.unravel_index(int(np.argmax(field)), field.shape), dtype=float)
    return ((index + 0.5) * spacing_A).tolist()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    summary_1x = json.loads(SUMMARY_1X.read_text(encoding="utf-8"))
    summary_2x = json.loads(SUMMARY_2X.read_text(encoding="utf-8"))
    properties = json.loads((PHASE_DIR / "thermophysical_candidates.json").read_text(encoding="utf-8"))
    melting_temperature_K = float(properties["ag"]["melting_temperature_K"])
    with np.load(GEOMETRY, allow_pickle=False) as archive:
        skeleton = np.asarray(archive["skeleton"], dtype=bool)
        size_A = np.asarray(archive["size_A"], dtype=float)
    temperature_1x, liquid_1x, metadata_1x = load_fields(CHECKPOINT_1X)
    temperature_2x, liquid_2x, metadata_2x = load_fields(CHECKPOINT_2X)
    if temperature_1x.shape != skeleton.shape:
        raise ValueError("The 1x checkpoint does not match the target geometry")
    expected_shape_2x = tuple(2 * value for value in skeleton.shape)
    if temperature_2x.shape != expected_shape_2x:
        raise ValueError("The 2x checkpoint is not an exact 2x spatial refinement")
    geometry_hash = sha256(GEOMETRY)
    if geometry_hash != summary_1x["geometry_sha256"] or geometry_hash != summary_2x["geometry_sha256"]:
        raise ValueError("Run summaries do not match the selected Ag 70 vol% geometry")

    block = 2
    coarsened_temperature = temperature_2x.reshape(
        skeleton.shape[0], block,
        skeleton.shape[1], block,
        skeleton.shape[2], block,
    ).mean(axis=(1, 3, 5))
    coarsened_liquid = liquid_2x.reshape(
        skeleton.shape[0], block,
        skeleton.shape[1], block,
        skeleton.shape[2], block,
    ).mean(axis=(1, 3, 5))
    ag_mask = ~skeleton
    delta = coarsened_temperature - temperature_1x
    common_ag_rise = float(np.sum(ag_mask * np.abs(temperature_1x - 300.0)))
    ag_count = int(np.count_nonzero(ag_mask))
    ag_mae = float(np.sum(ag_mask * np.abs(delta)) / ag_count)
    ag_rmse = float(np.sqrt(np.sum(ag_mask * delta**2) / ag_count))
    ag_mean_delta = float(np.sum(ag_mask * delta) / ag_count)
    ag_l1_over_rise = float(np.sum(ag_mask * np.abs(delta)) / max(common_ag_rise, 1.0e-300))
    max_abs_delta = float(np.max(np.abs(delta[ag_mask])))
    liquid_1x_mask = (liquid_1x >= 0.5) & ag_mask
    liquid_2x_mask = (coarsened_liquid >= 0.5) & ag_mask
    liquid_union = int(np.count_nonzero(liquid_1x_mask | liquid_2x_mask))
    spacing_1x = size_A / np.asarray(skeleton.shape)
    peak_xyz_1x = peak_coordinate(temperature_1x, spacing_1x)
    peak_xyz_2x = peak_coordinate(coarsened_temperature, spacing_1x)
    spacing_2x = size_A / np.asarray(temperature_2x.shape)
    peak_xyz_2x_fine = peak_coordinate(temperature_2x, spacing_2x)
    source_center_xy = size_A[:2] / 2.0

    def lateral_distance_from_source(coordinate: list[float]) -> float:
        difference = np.asarray(coordinate[:2], dtype=float) - source_center_xy
        difference -= np.round(difference / size_A[:2]) * size_A[:2]
        return float(np.linalg.norm(difference))

    report = {
        "status": "AG70VOL_2X_HEATING_WINDOW_COMPARISON_NOT_CONVERGENCE_PROOF",
        "geometry": str(GEOMETRY.relative_to(PROJECT_DIR)),
        "geometry_sha256": summary_1x["geometry_sha256"],
        "composition_target": "Ag:Ti3SiC2 = 70:30 by volume",
        "time_ps": 0.2,
        "comparison_method": "2x T and liquid-fraction fields block-averaged over 2x2x2 refined voxels and compared on the common Ag voxels of the 1x QSGS mask",
        "geometry_is_same_for_both_resolutions": summary_1x["geometry_sha256"] == summary_2x["geometry_sha256"],
        "grid_shape_xyz": {"1x": list(temperature_1x.shape), "2x": list(temperature_2x.shape)},
        "source_power_W": {"1x": summary_1x["source_power_W"], "2x": summary_2x["source_power_W"]},
        "relative_source_power_difference": summary_2x["source_power_W"] / summary_1x["source_power_W"] - 1.0,
        "peak_temperature_K": {
            "1x": summary_1x["max_temperature_K"],
            "2x": summary_2x["max_temperature_K"],
            "Ag_melting_temperature_candidate_K": melting_temperature_K,
            "1x_fraction_of_candidate_melting_temperature": summary_1x["max_temperature_K"] / melting_temperature_K,
            "2x_fraction_of_candidate_melting_temperature": summary_2x["max_temperature_K"] / melting_temperature_K,
            "2x_minus_1x_K": summary_2x["max_temperature_K"] - summary_1x["max_temperature_K"],
            "relative_to_1x_temperature_rise": (
                summary_2x["max_temperature_K"] - summary_1x["max_temperature_K"]
            ) / (summary_1x["max_temperature_K"] - 300.0),
        },
        "maximum_temperature_location_A": {
            "1x": peak_xyz_1x,
            "2x_fine_grid": peak_xyz_2x_fine,
            "2x_block_averaged": peak_xyz_2x,
        },
        "source_hotspot_alignment": {
            "source_center_xy_A": source_center_xy.tolist(),
            "1x_peak_lateral_distance_A": lateral_distance_from_source(peak_xyz_1x),
            "2x_fine_peak_lateral_distance_A": lateral_distance_from_source(peak_xyz_2x_fine),
            "2x_block_averaged_peak_lateral_distance_A": lateral_distance_from_source(peak_xyz_2x),
        },
        "Ag_weighted_temperature_field_difference": {
            "mean_signed_2x_minus_1x_K": ag_mean_delta,
            "mean_absolute_K": ag_mae,
            "root_mean_square_K": ag_rmse,
            "maximum_absolute_K": max_abs_delta,
            "L1_error_over_1x_Ag_temperature_rise": ag_l1_over_rise,
        },
        "liquid_fraction": {
            "1x_Ag_volume_fraction": summary_1x["ag_liquid_volume_fraction"],
            "2x_Ag_volume_fraction": summary_2x["ag_liquid_volume_fraction"],
            "connected_clusters_1x": summary_1x["connected_liquid_cluster_count_final"],
            "connected_clusters_2x": summary_2x["connected_liquid_cluster_count_final"],
            "threshold_mask_IoU": (
                float(np.count_nonzero(liquid_1x_mask & liquid_2x_mask) / liquid_union)
                if liquid_union else None
            ),
            "mask_comparison_note": "Both resolutions are fully solid at 0.2 ps; IoU is undefined because the union is empty." if not liquid_union else "IoU compares 1x liquid mask with the block-averaged 2x liquid occupancy mask.",
        },
        "energy_ledger": {
            "1x": {
                "injected_energy_J": summary_1x["injected_energy_J"],
                "total_enthalpy_J": summary_1x["total_enthalpy_J"],
                "balance_error_J": summary_1x["energy_balance_error_J"],
                "relative_balance_error": summary_1x["relative_energy_balance_error"],
            },
            "2x": {
                "injected_energy_J": summary_2x["injected_energy_J"],
                "total_enthalpy_J": summary_2x["total_enthalpy_J"],
                "balance_error_J": summary_2x["energy_balance_error_J"],
                "relative_balance_error": summary_2x["relative_energy_balance_error"],
            },
            "injected_energy_relative_difference": summary_2x["injected_energy_J"] / summary_1x["injected_energy_J"] - 1.0,
        },
        "checkpoint_integrity": {
            "1x_identity_sha256": metadata_1x["identity_sha256"],
            "2x_identity_sha256": metadata_2x["identity_sha256"],
            "1x_internal_energy_error_matches_summary": bool(np.isclose(
                summary_1x["energy_balance_error_J"], metadata_1x["energy_balance_error_J"], rtol=1.0e-10, atol=1.0e-30)),
            "2x_internal_energy_error_matches_summary": bool(np.isclose(
                summary_2x["energy_balance_error_J"], metadata_2x["energy_balance_error_J"], rtol=1.0e-10, atol=1.0e-30)),
        },
        "interpretation": (
            "The selected 70:30 vol% geometry has no Ag liquid in either the 1x or 2x Case A heating run through 0.2 ps. "
            "The 2x peak is lower than the 1x peak, while a single pair of resolutions does not establish full mesh convergence. "
            "This 0.2 ps window therefore cannot validate molten-pool size or lifetime for the selected composition."
        ),
        "limitations": [
            "Candidate thermal properties are not yet fully traced to archived source tables.",
            "Case A treats Ti3SiC2 as having zero heat capacity and zero conductivity; it is an idealized adiabatic scaffold bound.",
            "External domain faces are insulated and no environmental cooling or vaporization is modeled.",
        ],
    }
    path = RUN_DIR / "mesh2_ag70_comparison.json"
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
