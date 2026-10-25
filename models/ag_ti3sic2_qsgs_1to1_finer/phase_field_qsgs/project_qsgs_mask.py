#!/usr/bin/env python3
"""Conservatively project the saved QSGS mask to 5 Å and 10 Å grids."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from qsgs_projection import project_phase_mask


HERE = Path(__file__).resolve().parent
DEFAULT_INPUT = HERE.parent / "qsgs_morphology.npz"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=HERE)
    parser.add_argument("--spacings", type=float, nargs="+", default=[5.0, 10.0])
    args = parser.parse_args()

    with np.load(args.input, allow_pickle=False) as saved:
        skeleton_mask = saved["skeleton"].astype(bool)
        size_A = saved["size_A"].astype(float)
    if len(args.spacings) == 0:
        raise ValueError("Provide at least one requested grid spacing")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    source_voxel_volume = float(np.prod(size_A / np.asarray(skeleton_mask.shape)))
    total_volume = float(np.prod(size_A))
    source_skeleton_volume = float(skeleton_mask.sum() * source_voxel_volume)
    source_ag_volume = total_volume - source_skeleton_volume
    records = []

    for requested_spacing_A in args.spacings:
        result = project_phase_mask(skeleton_mask, size_A, requested_spacing_A)
        tag = f"{requested_spacing_A:g}A"
        output_path = args.output_dir / f"qsgs_phase_fractions_{tag}.npz"
        np.savez_compressed(
            output_path,
            ag_fraction=result.ag_fraction,
            skeleton_fraction=result.skeleton_fraction,
            size_A=result.size_A,
            spacing_A=result.spacing_A,
            requested_spacing_A=result.requested_spacing_A,
        )
        ag = result.ag_fraction
        ceramic = result.skeleton_fraction
        record = {
            "requested_spacing_A": requested_spacing_A,
            "grid_shape_xyz": list(ag.shape),
            "actual_spacing_A_xyz": result.spacing_A.tolist(),
            "qsgs_source_shape_xyz": list(skeleton_mask.shape),
            "qsgs_voxel_spacing_A_xyz": (size_A / np.asarray(skeleton_mask.shape)).tolist(),
            "total_material_volume_A3": total_volume,
            "ag_volume_fraction_source": source_ag_volume / total_volume,
            "ag_volume_fraction_projected": float(ag.mean()),
            "skeleton_volume_fraction_source": source_skeleton_volume / total_volume,
            "skeleton_volume_fraction_projected": float(ceramic.mean()),
            "ag_volume_error_A3": float(ag.sum() * np.prod(result.spacing_A) - source_ag_volume),
            "skeleton_volume_error_A3": float(ceramic.sum() * np.prod(result.spacing_A) - source_skeleton_volume),
            "ag_pure_cells": int(np.count_nonzero(ag >= 1.0 - 1.0e-12)),
            "skeleton_pure_cells": int(np.count_nonzero(ceramic >= 1.0 - 1.0e-12)),
            "mixed_cells": int(np.count_nonzero((ag > 1.0e-12) & (ceramic > 1.0e-12))),
            "min_ag_fraction": float(ag.min()),
            "max_ag_fraction": float(ag.max()),
            "output_npz": output_path.name,
        }
        records.append(record)

    report = {
        "status": "CONSERVATIVE_QSGS_GEOMETRY_PROJECTION_ONLY",
        "input_npz": args.input.name,
        "source_skeleton_volume_fraction": source_skeleton_volume / total_volume,
        "source_ag_volume_fraction": source_ag_volume / total_volume,
        "projection_method": "exact axis-aligned overlap averaging of source voxels; no binary threshold",
        "records": records,
        "limitations": [
            "This is a geometry map, not a thermal or phase-change simulation.",
            "Continuous cell fractions preserve phase volume but do not prove coarse-grid percolation.",
            "No Ag/scaffold interface heat-transfer coefficient or boundary condition is inferred.",
            "No latent heat, lattice heat capacity, thermal conductivity, or phase kinetics are assigned.",
        ],
    }
    report_path = args.output_dir / "qsgs_projection_report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
