#!/usr/bin/env python3
"""Compare same-geometry Ag70 Case A/B checkpoints near 14.51 ps."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
RUN_DIR = HERE.parent.parent
PROJECT_DIR = RUN_DIR.parents[1]
PHASE_DIR = PROJECT_DIR / "phase_field_qsgs"
GEOMETRY = RUN_DIR / "qsgs_morphology_ag70vol.npz"
CASE_A = HERE / "checkpoints/caseA_ag70_1x_heat_match_14p51ps.npz"
CASE_B = HERE / "checkpoints/caseB_ag70_1x_heat_final_14p506ps.npz"
OUTPUT = HERE / "caseA_caseB_matched_14p51ps.csv"


def main() -> None:
    import sys

    sys.path.insert(0, str(PHASE_DIR))
    from qsgs_liquid_tracking import liquid_cluster_labels

    with np.load(GEOMETRY, allow_pickle=False) as data:
        scaffold = np.asarray(data["skeleton"], dtype=bool)
        size_A = np.asarray(data["size_A"], dtype=float)
    ag = ~scaffold
    voxel_nm3 = float(np.prod(size_A / np.asarray(ag.shape)) / 1000.0)

    rows = []
    for case, path in (("A", CASE_A), ("B", CASE_B)):
        with np.load(path, allow_pickle=False) as data:
            temperature = np.asarray(data["temperature_K"], dtype=float)
            liquid = np.asarray(data["ag_liquid_fraction"], dtype=float)
            metadata = json.loads(str(data["metadata_json"].item()))
        mask = ag & (liquid >= 0.5)
        labels, cluster_count = liquid_cluster_labels(mask)
        counts = np.bincount(labels.ravel(), minlength=cluster_count + 1)
        counts[0] = 0
        largest_voxels = int(counts.max()) if counts.size > 1 else 0
        largest_cluster_nm3 = largest_voxels * voxel_nm3
        row = {
            "case": case,
            "boundary_case": metadata["identity"]["boundary_case"],
            "time_ps": float(metadata["time_ps"]),
            "grid": "x".join(str(v) for v in metadata["identity"]["grid_shape_xyz"]),
            "geometry_sha256": metadata["identity"]["geometry_sha256"],
            "source_actual_power_W": float(metadata["actual_source_power_W"]),
            "source_peak_flux_W_m2": float(metadata["identity"]["source_peak_flux_W_m2"]),
            "source_gaussian_radius_A": float(metadata["identity"]["source_gaussian_radius_A"]),
            "source_absorption_depth_A": float(metadata["identity"]["source_absorption_depth_A"]),
            "cumulative_source_energy_J": float(metadata["injected_energy_J"]),
            "total_enthalpy_J": float(metadata["total_enthalpy_J"]),
            "relative_energy_error": float(metadata["relative_energy_balance_error"]),
            "ag_peak_temperature_K": float(temperature[ag].max()),
            "ag_mean_temperature_K": float(temperature[ag].mean()),
            "scaffold_peak_temperature_K": float(temperature[scaffold].max()),
            "scaffold_mean_temperature_K": float(temperature[scaffold].mean()),
            "ag_liquid_fraction_of_ag": float(np.sum(liquid[ag]) / ag.sum()),
            "liquid_volume_fraction_of_domain": float(np.sum(liquid[ag]) / liquid.size),
            "liquid_volume_all_fractions_nm3": float(np.sum(liquid[ag]) * voxel_nm3),
            "connected_cluster_count_f_l_ge_0p5": int(cluster_count),
            "largest_cluster_voxels_f_l_ge_0p5": largest_voxels,
            "largest_cluster_volume_nm3_f_l_ge_0p5": largest_cluster_nm3,
            "scaffold_stored_enthalpy_J": float(metadata.get("scaffold_stored_enthalpy_J", 0.0)),
            "checkpoint": path.name,
        }
        # The Case B checkpoint stores phase-resolved energy in the Case B summary.
        if case == "B":
            summary = json.loads((HERE / "caseB_heat_summary.json").read_text())
            row["scaffold_stored_enthalpy_J"] = float(summary["final_metrics"]["ti3sic2_stored_enthalpy_J"])
        else:
            row["scaffold_stored_enthalpy_J"] = 0.0
        rows.append(row)

    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
