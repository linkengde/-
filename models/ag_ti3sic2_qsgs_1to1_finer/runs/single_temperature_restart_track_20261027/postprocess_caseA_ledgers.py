"""Add derived latent-energy fields and index restart checkpoints for Case A."""

from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

import numpy as np


RUN_DIR = Path(__file__).resolve().parent
PROJECT_DIR = RUN_DIR.parents[1]
sys.path.insert(0, str(PROJECT_DIR / "phase_field_qsgs"))
from thermophysical_model import AgEnthalpy  # noqa: E402


def _read_geometry() -> tuple[np.ndarray, np.ndarray, float, float]:
    geometry_path = PROJECT_DIR / "qsgs_morphology.npz"
    with np.load(geometry_path, allow_pickle=False) as archive:
        ag_mask = ~np.asarray(archive["skeleton"], dtype=bool)
        size_A = np.asarray(archive["size_A"], dtype=float)
    cell_volume_m3 = float(np.prod(size_A / np.asarray(ag_mask.shape)) * 1.0e-30)
    ag = AgEnthalpy()
    return ag_mask, size_A, cell_volume_m3, ag.reference_density_kg_m3 * ag.latent_heat_J_kg


AG_MASK, DOMAIN_SIZE_A, CELL_VOLUME_M3, AG_LATENT_CAPACITY_PER_M3 = _read_geometry()
AG_CELL_COUNT = int(AG_MASK.sum())


def _field_ledger(liquid_fraction: np.ndarray, total_energy_J: float, injected_energy_J: float) -> dict[str, float]:
    latent_energy = float(
        np.sum(AG_MASK * AG_LATENT_CAPACITY_PER_M3 * liquid_fraction) * CELL_VOLUME_M3
    )
    error = float(total_energy_J - injected_energy_J)
    return {
        "ag_latent_energy_J": latent_energy,
        "sensible_enthalpy_J": float(total_energy_J - latent_energy),
        "energy_balance_relative_error": error / injected_energy_J if injected_energy_J else 0.0,
    }


def _atomic_npz(path: Path, arrays: dict[str, np.ndarray]) -> None:
    temporary = path.with_name(path.name + ".tmp")
    try:
        with temporary.open("wb") as handle:
            np.savez_compressed(handle, **arrays)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def enrich_checkpoint(path: Path) -> dict[str, object]:
    with np.load(path, allow_pickle=False) as archive:
        arrays = {name: np.asarray(archive[name]).copy() for name in archive.files}
    metadata = json.loads(str(arrays["metadata_json"].item()))
    identity = metadata["identity"]
    tile_xy = int(identity.get("tile_xy", 1))
    ag_mask = np.tile(AG_MASK, (tile_xy, tile_xy, 1))
    size_A = DOMAIN_SIZE_A * np.asarray([tile_xy, tile_xy, 1.0])
    temperature = arrays["temperature_K"]
    liquid = arrays["ag_liquid_fraction"]
    enthalpy = arrays["enthalpy_density_J_m3"]
    cell_volume = float(np.prod(size_A / np.asarray(ag_mask.shape)) * 1.0e-30)
    total_energy = float(enthalpy.sum() * cell_volume)
    injected = float(metadata["injected_energy_J"])
    latent = float(np.sum(ag_mask * AG_LATENT_CAPACITY_PER_M3 * liquid) * cell_volume)
    error = total_energy - injected
    metadata.update(
        {
            "total_enthalpy_J": total_energy,
            "ag_latent_energy_J": latent,
            "sensible_enthalpy_J": total_energy - latent,
            "energy_balance_error_J": error,
            "relative_energy_balance_error": error / injected if injected else 0.0,
            "derived_ledger_note": "Latent energy is the Ag reference-density fusion latent heat multiplied by the saved liquid fraction field; sensible enthalpy is total enthalpy minus this value.",
        }
    )
    arrays["metadata_json"] = np.asarray(json.dumps(metadata))
    _atomic_npz(path, arrays)

    ag_count = int(ag_mask.sum())
    liquid_volume_fraction = float(np.sum(ag_mask * liquid) / ag_count) if ag_count else 0.0
    tracker_metadata = metadata.get("tracker_metadata") or {}
    history = tracker_metadata.get("history", [])
    current_time = float(metadata["time_ps"])
    main_rows = [
        row for row in history
        if bool(row.get("is_main_cluster")) and abs(float(row["time_ps"]) - current_time) <= 1.0e-9
    ]
    main = main_rows[-1] if main_rows else {}
    tracker_labels = arrays["tracker_track_labels"]
    active_sizes = np.bincount(tracker_labels.ravel().astype(np.int64)) if tracker_labels.size else np.zeros(1)
    if active_sizes.size:
        active_sizes[0] = 0
    return {
        "checkpoint": str(path.relative_to(RUN_DIR)),
        "time_ps": current_time,
        "max_temperature_K": float(temperature.max()),
        "max_temperature_seen_K": float(metadata["max_temperature_seen_K"]),
        "mean_ag_temperature_K": float(np.sum(ag_mask * temperature) / ag_count) if ag_count else 300.0,
        "ag_liquid_volume_fraction": liquid_volume_fraction,
        "max_ag_liquid_fraction": float(liquid.max()),
        "connected_cluster_count": int(np.count_nonzero(active_sizes)),
        "largest_cluster_voxels": int(active_sizes.max()) if active_sizes.size else 0,
        "largest_main_track_id": main.get("track_id"),
        "largest_main_track_duration_ps": main.get("continuous_duration_ps"),
        "largest_main_cluster_volume_A3": main.get("volume_A3"),
        "largest_main_cluster_centroid_A": json.dumps(main.get("centroid_A")) if main else "",
        "largest_main_cluster_depth_from_top_A": main.get("depth_from_top_A"),
        "largest_main_cluster_lateral_size_x_A": main.get("lateral_size_x_A"),
        "largest_main_cluster_lateral_size_y_A": main.get("lateral_size_y_A"),
        "injected_energy_J": injected,
        "total_enthalpy_J": total_energy,
        "ag_latent_energy_J": latent,
        "sensible_enthalpy_J": total_energy - latent,
        "energy_balance_error_J": error,
        "relative_energy_balance_error": error / injected if injected else 0.0,
    }


def enrich_timeseries(path: Path) -> None:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
        fields = list(rows[0]) if rows else []
    if not rows:
        return
    extras = ("ag_latent_energy_J", "sensible_enthalpy_J", "energy_balance_relative_error")
    fields = [field for field in fields if field not in extras] + list(extras)
    for row in rows:
        liquid_fraction = float(row["ag_liquid_volume_fraction"])
        total_energy = float(row["total_energy_J"])
        injected = float(row["injected_energy_J"])
        latent_capacity = float(AG_CELL_COUNT * CELL_VOLUME_M3 * AG_LATENT_CAPACITY_PER_M3)
        latent = liquid_fraction * latent_capacity
        error = float(row["energy_balance_error_J"])
        row.update(
            {
                "ag_latent_energy_J": latent,
                "sensible_enthalpy_J": total_energy - latent,
                "energy_balance_relative_error": error / injected if injected else 0.0,
            }
        )
    temporary = path.with_name(path.name + ".tmp")
    try:
        with temporary.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _update_candidate_summaries() -> None:
    heat_json = RUN_DIR / "caseA_heat0p6.json"
    heat_csv = RUN_DIR / "caseA_heat0p6.csv"
    heat_clusters = RUN_DIR / "caseA_heat0p6.clusters.csv"
    summary = json.loads(heat_json.read_text())
    records = list(csv.DictReader(heat_clusters.open(newline="")))
    main_records = [row for row in records if row["is_main_cluster"].lower() == "true"]
    longest_main = max((float(row["continuous_duration_ps"]) for row in main_records), default=0.0)
    screen = summary["tracked_liquid_persistence_screen"]
    screen["longest_track_while_it_was_the_largest_component_ps"] = longest_main
    screen["provisional_0p3_ps_gate_met"] = longest_main >= 0.3
    screen["provisional_0p5_ps_enhanced_gate_met"] = longest_main >= 0.5
    screen["interpretation"] = (
        "Heating-only screen for the largest component's persistent identity. IDs continue by voxel overlap: "
        "the largest-overlap child inherits an ID on split and the largest-overlap parent supplies the ID on merge. "
        "The 0.3/0.5 ps thresholds are provisional screens, not universal thermal-stability criteria."
    )
    summary["liquid_tracking"] = {
        "threshold_ag_liquid_fraction": summary["liquid_connectivity_threshold"],
        "track_count": len({int(row["track_id"]) for row in records}),
        "merge_component_event_count": sum("merge" in row["event"] for row in records),
        "split_parent_event_count": len(
            {
                (float(row["time_ps"]), int(parent_id))
                for row in records if "split" in row["event"]
                for parent_id in json.loads(row["parent_track_ids"])
            }
        ),
        "cluster_history_csv": heat_clusters.name,
    }
    summary["energy_ledger_definition"] = {
        "total_enthalpy_J": "Integral of the stored volumetric enthalpy field relative to uniform 300 K.",
        "ag_latent_energy_J": "Ag reference-density fusion latent heat multiplied by the stored Ag liquid fraction field.",
        "sensible_enthalpy_J": "Total enthalpy minus the Ag latent-energy term.",
        "relative_energy_balance_error": "(total enthalpy - cumulative injected energy) / cumulative injected energy.",
    }
    heat_summary_rows = list(csv.DictReader(heat_csv.open(newline="")))
    for branch, cutoff in (("cooling_from_0p5", 0.5), ("cooling_from_0p6", 0.6)):
        json_path = RUN_DIR / branch / ("caseA_off0p5_cool0p5.json" if cutoff == 0.5 else "caseA_off0p6_cool0p6.json")
        if not json_path.exists():
            continue
        branch_summary = json.loads(json_path.read_text())
        heating_state = min(heat_summary_rows, key=lambda row: abs(float(row["time_ps"]) - cutoff))
        heat_track_age = float(heating_state["largest_liquid_track_duration_ps"])
        cluster_history = list(
            csv.DictReader((json_path.parent / branch_summary["liquid_tracking"]["cluster_history_csv"]).open(newline=""))
        )
        total_track_age = max(
            (float(row["continuous_duration_ps"]) for row in cluster_history), default=0.0
        )
        branch_summary["heat_off_candidate_screen"] = {
            "candidate_heat_off_time_ps": cutoff,
            "largest_track_id_at_heat_off": int(heating_state["largest_liquid_track_id"]),
            "same_track_age_at_heat_off_ps": heat_track_age,
            "provisional_0p3_ps_heating_gate_met": heat_track_age >= 0.3,
            "provisional_0p5_ps_heating_gate_met": heat_track_age >= 0.5,
            "total_longest_track_age_through_cooling_ps": total_track_age,
            "interpretation": "The provisional persistence gates apply at the heat-off instant; longer total lifetime includes post-cutoff solidification and is not a steady-pool stability claim.",
        }
        branch_screen = branch_summary["tracked_liquid_persistence_screen"]
        branch_screen["longest_total_track_age_through_cooling_ps"] = total_track_age
        branch_screen["longest_single_track_sampled_duration_ps"] = heat_track_age
        branch_screen["provisional_0p3_ps_gate_met"] = heat_track_age >= 0.3
        branch_screen["provisional_0p5_ps_enhanced_gate_met"] = heat_track_age >= 0.5
        branch_screen["interpretation"] = (
            "Persistence gates are evaluated at heat-off and use the largest persistent track at that instant. "
            "Total track age through cooling is separately recorded and does not mean a steady molten pool."
        )
        branch_summary["energy_ledger_definition"] = summary["energy_ledger_definition"]
        json_path.write_text(json.dumps(branch_summary, indent=2) + "\n")
    heat_json.write_text(json.dumps(summary, indent=2) + "\n")
    half_json = RUN_DIR / "time_step_halfdt/caseA_heat0p5_maxdt0p02fs.json"
    if half_json.exists():
        half_summary = json.loads(half_json.read_text())
        half_summary["energy_ledger_definition"] = summary["energy_ledger_definition"]
        half_json.write_text(json.dumps(half_summary, indent=2) + "\n")
    freeze_json = RUN_DIR / "cooling_from_0p5/freeze_extension/caseA_off0p5_freeze0p86.json"
    if freeze_json.exists():
        freeze_summary = json.loads(freeze_json.read_text())
        freeze_summary["energy_ledger_definition"] = summary["energy_ledger_definition"]
        freeze_json.write_text(json.dumps(freeze_summary, indent=2) + "\n")


def main() -> None:
    csv_files = [
        RUN_DIR / "caseA_heat0p6.csv",
        RUN_DIR / "time_step_halfdt/caseA_heat0p5_maxdt0p02fs.csv",
        RUN_DIR / "cooling_from_0p5/caseA_off0p5_cool0p5.csv",
        RUN_DIR / "cooling_from_0p5/freeze_extension/caseA_off0p5_freeze0p86.csv",
        RUN_DIR / "cooling_from_0p6/caseA_off0p6_cool0p6.csv",
    ]
    for path in csv_files:
        if path.exists():
            enrich_timeseries(path)
    checkpoints = sorted(RUN_DIR.rglob("*.npz"))
    index = [enrich_checkpoint(path) for path in checkpoints]
    index_path = RUN_DIR / "checkpoint_index.csv"
    with index_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(index[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(index)
    _update_candidate_summaries()
    print(f"updated {len(csv_files)} energy timeseries and {len(index)} restart checkpoints")


if __name__ == "__main__":
    main()
