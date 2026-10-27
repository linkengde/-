"""Summarize finite-size heating and cooling results for the 3x3 Case A domain."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np


RUN_DIR = Path(__file__).resolve().parent
PROJECT_DIR = RUN_DIR.parents[1]
with np.load(PROJECT_DIR / "qsgs_morphology.npz", allow_pickle=False) as _geometry:
    AG_MASK_BASE = ~np.asarray(_geometry["skeleton"], dtype=bool)
AG_MASK_TILE3 = np.tile(AG_MASK_BASE, (3, 3, 1))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def summarize(path: Path) -> dict[str, object]:
    rows = read_csv(path)
    if not rows:
        raise RuntimeError(f"No rows in completed run: {path}")
    liquid_rows = [row for row in rows if float(row["ag_liquid_volume_fraction"]) > 0.0]
    solid_rows = [row for row in rows if float(row["ag_liquid_volume_fraction"]) == 0.0]
    connected_rows = [
        row for row in rows if int(row.get("connected_liquid_cluster_count", "0")) > 0
    ]
    max_abs_energy_error = max(abs(float(row["energy_balance_error_J"])) for row in rows)
    max_relative_error = max(
        abs(float(row["energy_balance_error_J"]))
        / max(abs(float(row["injected_energy_J"])), 1e-300)
        for row in rows
    )
    final = rows[-1]
    return {
        "csv": str(path.relative_to(RUN_DIR)),
        "sample_count": len(rows),
        "start_time_ps": float(rows[0]["time_ps"]),
        "end_time_ps": float(final["time_ps"]),
        "first_nonzero_liquid_fraction_sample_ps": (
            float(liquid_rows[0]["time_ps"]) if liquid_rows else None
        ),
        "last_liquid_sample_ps": float(liquid_rows[-1]["time_ps"]) if liquid_rows else None,
        "first_fully_solid_sample_ps": float(solid_rows[0]["time_ps"]) if solid_rows else None,
        "first_connected_liquid_cluster_sample_ps": (
            float(connected_rows[0]["time_ps"]) if connected_rows else None
        ),
        "first_no_f_l_ge_0p5_cluster_sample_ps": (
            float(next(row["time_ps"] for row in rows
                       if int(row.get("connected_liquid_cluster_count", "0")) == 0))
            if any(int(row.get("connected_liquid_cluster_count", "0")) == 0 for row in rows)
            else None
        ),
        "final_peak_temperature_K": float(final["max_temperature_K"]),
        "final_Ag_liquid_volume_fraction": float(final["ag_liquid_volume_fraction"]),
        "final_connected_liquid_cluster_count": int(
            final.get("connected_liquid_cluster_count", "0")
        ),
        "injected_energy_J": float(final["injected_energy_J"]),
        "final_total_enthalpy_J": float(final.get("total_energy_J", "nan")),
        "max_abs_energy_balance_error_J": max_abs_energy_error,
        "max_relative_energy_balance_error": max_relative_error,
    }


def checkpoint_summary(path: Path) -> dict[str, object]:
    with np.load(path, allow_pickle=False) as archive:
        metadata = json.loads(str(archive["metadata_json"].item()))
        temperature = np.asarray(archive["temperature_K"])
        liquid = np.asarray(archive["ag_liquid_fraction"])
        enthalpy = np.asarray(archive["enthalpy_density_J_m3"])
    shape = temperature.shape
    if shape != (300, 300, 50):
        raise RuntimeError(f"Unexpected tile shape in checkpoint: {path}: {shape}")
    identity = metadata["identity"]
    if (identity["tile_xy"] != 3 or identity["boundary_case"] != "A"
            or identity["geometry_sha256"] != "7e0e1509b414597ce4e2a1c9fc387ed9656caf507f1186e2d28de79c3c1bfc43"):
        raise RuntimeError(f"Checkpoint model identity differs from verified 3x3 Case A: {path}")
    cell_volume_m3 = float(np.prod(np.asarray(identity["spacing_A"], dtype=float) * 1.0e-10))
    enthalpy_from_field_J = float(enthalpy.sum() * cell_volume_m3)
    energy_error_from_field_J = enthalpy_from_field_J - float(metadata["injected_energy_J"])
    if abs(energy_error_from_field_J - float(metadata["energy_balance_error_J"])) > 1e-28:
        raise RuntimeError(f"Checkpoint energy ledger does not match its fields: {path}")
    roi = (slice(100, 200), slice(100, 200), slice(None))
    ag = AG_MASK_TILE3[roi]
    central_temperature = temperature[roi]
    central_liquid = liquid[roi]
    return {
        "checkpoint": str(path.relative_to(RUN_DIR)),
        "time_ps": float(metadata["time_ps"]),
        "boundary_case": identity["boundary_case"],
        "tile_xy": identity["tile_xy"],
        "whole_domain_peak_Ag_temperature_K": float(np.max(temperature[AG_MASK_TILE3])),
        "central_cell_peak_Ag_temperature_K": float(np.max(central_temperature[ag])),
        "central_cell_Ag_liquid_fraction": float(np.sum(central_liquid * ag) / ag.sum()),
        "central_cell_liquid_Ag_volume_fraction_ge_0p5": float(
            np.count_nonzero((central_liquid >= 0.5) & ag) / ag.sum()
        ),
        "whole_domain_Ag_liquid_fraction": float(
            np.sum(liquid * AG_MASK_TILE3) / AG_MASK_TILE3.sum()
        ),
        "total_enthalpy_J_from_field": enthalpy_from_field_J,
        "energy_balance_error_J_from_field": energy_error_from_field_J,
        "injected_energy_J": float(metadata["injected_energy_J"]),
        "energy_balance_error_J_metadata": float(metadata["energy_balance_error_J"]),
        "latent_energy_J": float(metadata["ag_latent_energy_J"]),
    }


def cluster_snapshots(path: Path, cutoff_ps: float) -> list[dict[str, object]]:
    if not path.exists():
        return []
    rows = read_csv(path)
    times = [cutoff_ps, cutoff_ps + 0.1, cutoff_ps + 0.2, cutoff_ps + 0.3,
             cutoff_ps + 0.4, cutoff_ps + 0.5]
    output = []
    for target in times:
        nearest = min(rows, key=lambda row: abs(float(row["time_ps"]) - target)) if rows else None
        if nearest is None or abs(float(nearest["time_ps"]) - target) > 0.0051:
            continue
        same_time = [
            row for row in rows
            if abs(float(row["time_ps"]) - float(nearest["time_ps"])) < 1e-10
        ]
        main = next((row for row in same_time if row["is_main_cluster"].lower() == "true"), None)
        output.append({
            "time_ps": float(nearest["time_ps"]),
            "offset_after_heat_off_ps": float(nearest["time_ps"]) - cutoff_ps,
            "main_cluster": None if main is None else {
                "track_id": int(main["track_id"]),
                "event": main["event"],
                "continuous_duration_ps": float(main["continuous_duration_ps"]),
                "voxel_count": int(main["voxel_count"]),
                "volume_A3": float(main["volume_A3"]),
                "centroid_A": json.loads(main["centroid_A"]),
                "depth_from_top_A": float(main["depth_from_top_A"]),
                "lateral_size_x_A": float(main["lateral_size_x_A"]),
                "lateral_size_y_A": float(main["lateral_size_y_A"]),
            },
            "connected_liquid_cluster_count_at_sample": len(same_time),
        })
    return output


def main() -> None:
    heating = json.loads((RUN_DIR / "tile3_comparison.json").read_text())
    cool5 = RUN_DIR / "cooling_from_0p5/caseA_tile3_off0p5_cool0p5.csv"
    cool6 = RUN_DIR / "cooling_from_0p6/caseA_tile3_off0p6_cool0p6.csv"
    extension_dir = RUN_DIR / "cooling_from_0p5/freeze_extension"
    checkpoints = sorted(
        list((RUN_DIR / "heat/checkpoints").glob("*.npz"))
        + list((RUN_DIR / "cooling_from_0p5/checkpoints").glob("*.npz"))
        + list((RUN_DIR / "cooling_from_0p6/checkpoints").glob("*.npz"))
        + list((RUN_DIR / "cooling_from_0p5/freeze_extension/checkpoints").glob("*.npz"))
        + list((RUN_DIR / "cooling_from_0p5/freeze_extension").glob("*.checkpoint.npz"))
        + list((RUN_DIR / "cooling_from_0p5").glob("*.checkpoint.npz"))
        + list((RUN_DIR / "cooling_from_0p6").glob("*.checkpoint.npz"))
    )
    # Exclude transient atomic-write files if any appear during a running job.
    checkpoints = [p for p in checkpoints if not p.name.endswith(".tmp.npz")]
    checkpoint_rows = [checkpoint_summary(path) for path in checkpoints]
    index_fields = [
        "checkpoint", "time_ps", "boundary_case", "tile_xy",
        "whole_domain_peak_Ag_temperature_K", "central_cell_peak_Ag_temperature_K",
        "central_cell_Ag_liquid_fraction", "central_cell_liquid_Ag_volume_fraction_ge_0p5",
        "whole_domain_Ag_liquid_fraction", "total_enthalpy_J_from_field",
        "injected_energy_J", "energy_balance_error_J_from_field", "latent_energy_J",
    ]
    with (RUN_DIR / "tile3_checkpoint_index.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=index_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row[field] for field in index_fields} for row in checkpoint_rows)
    output = {
        "status": "3X3_FINITE_SIZE_HEAT_AND_COOLING_SUMMARY",
        "heating_comparison": heating,
        "cooling_from_0p5ps": summarize(cool5),
        "cooling_from_0p6ps": summarize(cool6),
        "cooling_main_cluster_samples_from_0p5ps": cluster_snapshots(
            RUN_DIR / "cooling_from_0p5/caseA_tile3_off0p5_cool0p5.clusters.csv", 0.5
        ),
        "cooling_main_cluster_samples_from_0p6ps": cluster_snapshots(
            RUN_DIR / "cooling_from_0p6/caseA_tile3_off0p6_cool0p6.clusters.csv", 0.6
        ),
        "freeze_extension_csv": (
            str((extension_dir / "caseA_tile3_off0p5_to0p86.csv").relative_to(RUN_DIR))
            if (extension_dir / "caseA_tile3_off0p5_to0p86.csv").exists() else None
        ),
        "checkpoint_index_csv": "tile3_checkpoint_index.csv",
        "checkpoint_snapshots": checkpoint_rows,
        "interpretation": (
            "Global tile and central-cell quantities are reported separately. The former uses the full "
            "nine-cell domain; central-cell liquid volume fractions use only the central native cell. "
            "Case A treats Ti3SiC2 as zero-capacity/zero-conductivity insulation, and adiabatic cooling "
            "redistributes energy internally rather than losing heat to an environment."
        ),
    }
    output_path = RUN_DIR / "tile3_cooling_summary.json"
    output_path.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
