#!/usr/bin/env python3
"""Audit PTM class counts and Ag-only coordination in the corrected Round7 subset."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import ovito
from ovito.io import import_file
from ovito.modifiers import PolyhedralTemplateMatchingModifier, SelectTypeModifier
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
RUNS = HERE.parent
RUN_NEW = RUNS / "ttm_mtp_round7_continue_1p4_to_1p6ps_20261109"
RUN_OLD = RUNS / "ttm_mtp_round7_continue_to_1p5ps_20261109"
TIMES = [("0p6", 0.60), ("0p85", 0.85), ("1p1", 1.10), ("1p35", 1.35), ("1p6", 1.60)]
BOX_XY = 200.0
COORD_CUTOFF_A = 4.0


def read_subset():
    with (HERE / "round7_1p6ps_source_core_ag.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    ids = np.asarray([int(row["atom_id"]) for row in rows], dtype=np.int64)
    candidates = np.asarray([int(row["joint_candidate"]) == 1 for row in rows], dtype=bool)
    return ids, candidates


def checkpoint(tag: str) -> Path:
    filename = f"atoms_{tag}ps.dump"
    path = RUN_NEW / "checkpoints" / filename
    if not path.is_file():
        path = RUN_OLD / "checkpoints" / filename
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def ptm_counts(path: Path, subset_ids: np.ndarray) -> dict[str, int]:
    pipeline = import_file(str(path))
    pipeline.modifiers.append(SelectTypeModifier(types={4}))
    ptm = PolyhedralTemplateMatchingModifier()
    ptm.only_selected = True
    pipeline.modifiers.append(ptm)
    data = pipeline.compute()
    ids = np.asarray(data.particles["Particle Identifier"][:], dtype=np.int64)
    structure = np.asarray(data.particles["Structure Type"][:], dtype=np.int32)
    order = np.argsort(ids)
    ids, structure = ids[order], structure[order]
    ix = np.searchsorted(ids, subset_ids)
    if np.any(ix >= len(ids)) or not np.array_equal(ids[ix], subset_ids):
        raise RuntimeError(f"PTM particle IDs do not match the corrected subset at {path}")
    values = structure[ix]
    return {
        "FCC": int(np.count_nonzero(values == 1)),
        "HCP": int(np.count_nonzero(values == 2)),
        "BCC": int(np.count_nonzero(values == 3)),
        "Other": int(np.count_nonzero(values == 0)),
        "other_codes": {str(int(code)): int(count) for code, count in zip(*np.unique(values, return_counts=True))},
    }


def coordination_counts(path: Path, subset_ids: np.ndarray, candidates: np.ndarray) -> dict[str, dict[str, float | int]]:
    data = np.loadtxt(path, skiprows=9)
    ids = data[:, 0].astype(np.int64)
    types = data[:, 1].astype(np.int32)
    ag_mask = types == 4
    ag_ids = ids[ag_mask]
    ag_xyz = data[ag_mask, 2:5]
    order = np.argsort(ag_ids)
    ag_ids, ag_xyz = ag_ids[order], ag_xyz[order]
    subset_ix = np.searchsorted(ag_ids, subset_ids)
    if np.any(subset_ix >= len(ag_ids)) or not np.array_equal(ag_ids[subset_ix], subset_ids):
        raise RuntimeError(f"Corrected IDs not found among Ag at {path}")

    points = ag_xyz.copy()
    points[:, 0] %= BOX_XY
    points[:, 1] %= BOX_XY
    zlo = float(points[:, 2].min())
    zspan = float(np.ptp(points[:, 2]))
    zperiod = max(zspan + 2 * COORD_CUTOFF_A + 1.0e-6, 4 * COORD_CUTOFF_A)
    points[:, 2] = points[:, 2] - zlo + COORD_CUTOFF_A
    tree = cKDTree(points, boxsize=[BOX_XY, BOX_XY, zperiod])
    counts = np.fromiter(
        (len(neighbors) - 1 for neighbors in tree.query_ball_point(points[subset_ix], COORD_CUTOFF_A)),
        dtype=np.int32,
        count=len(subset_ix),
    )

    result = {}
    for label, mask in (("joint_candidates", candidates), ("other_subset", ~candidates), ("all_subset", np.ones(len(candidates), dtype=bool))):
        subset = counts[mask]
        result[label] = {
            "atoms": int(len(subset)),
            "mean": float(subset.mean()) if len(subset) else 0.0,
            "median": float(np.median(subset)) if len(subset) else 0.0,
            "less_than_8": int(np.count_nonzero(subset < 8)),
            "less_than_10": int(np.count_nonzero(subset < 10)),
            "less_than_12": int(np.count_nonzero(subset < 12)),
        }
    return result


def main():
    subset_ids, candidates = read_subset()
    report = {
        "status": "exploratory_diagnostic",
        "ovito_version": ovito.version_string,
        "subset_atoms": int(len(subset_ids)),
        "joint_candidate_atoms": int(candidates.sum()),
        "coordination_definition": "Ag-only neighbors within 4.0 A; periodic x/y (200 A), open z.",
        "ptm_definition": "OVITO PTM on Ag-selected particles only; runtime default RMSD cutoff retained.",
        "frames": [],
    }
    for tag, time_ps in TIMES:
        path = checkpoint(tag)
        frame = {
            "time_ps": time_ps,
            "ptm_class_counts": ptm_counts(path, subset_ids),
            "ag_coordination": coordination_counts(path, subset_ids, candidates),
        }
        report["frames"].append(frame)
        print(f"diagnosed {time_ps:.2f} ps", flush=True)
    (HERE / "ptm_coordination_diagnostic.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("wrote ptm_coordination_diagnostic.json")


if __name__ == "__main__":
    main()
