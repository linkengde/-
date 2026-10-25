#!/usr/bin/env python3
"""Summarize initial forces by QSGS skeleton voxel-boundary depth."""
from __future__ import annotations
import gzip
import json
from pathlib import Path
import numpy as np
from scipy import ndimage

HERE = Path(__file__).resolve().parent
MODEL = HERE.parents[1]
FORCE_DUMP = HERE / "initial_forces.dump.gz"
MORPH = MODEL / "qsgs_morphology.npz"
OUT = HERE / "force_by_skeleton_surface_depth.json"
SYMBOL = {1: "Ti", 2: "Si", 3: "C", 4: "Ag"}
BINS = [(0.0, 2.0, "0-2"), (2.0, 4.0, "2-4"), (4.0, 8.0, "4-8"), (8.0, float("inf"), "8+")]


def stats(v: np.ndarray) -> dict:
    if not len(v):
        return {"n": 0}
    return {
        "n": int(len(v)),
        "mean_eV_A": float(np.mean(v)),
        "p50_eV_A": float(np.percentile(v, 50)),
        "p95_eV_A": float(np.percentile(v, 95)),
        "p99_eV_A": float(np.percentile(v, 99)),
        "max_eV_A": float(np.max(v)),
        "n_force_gt_2_eV_A": int(np.count_nonzero(v > 2.0)),
        "fraction_force_gt_2": float(np.mean(v > 2.0)),
    }


with gzip.open(FORCE_DUMP, "rt") as f:
    rows = np.loadtxt(f, skiprows=9, dtype=np.float64)
ids = rows[:, 0].astype(np.int64)
types = rows[:, 1].astype(np.int8)
pos = rows[:, 2:5]
forces = np.linalg.norm(rows[:, 5:8], axis=1)
with np.load(MORPH) as npz:
    mask = npz["skeleton"].astype(bool)
    size = npz["size_A"].astype(float)
shape = np.asarray(mask.shape)
# Match build_model.boundary_depth(): tile x/y to respect periodic morphology,
# while keeping z open, then sample the distance on each atom's QSGS voxel.
tiled = np.tile(mask, (3, 3, 1))
depth_grid = ndimage.distance_transform_edt(tiled, sampling=size / shape)
depth_grid = depth_grid[shape[0]:2*shape[0], shape[1]:2*shape[1], :]
ijk = np.floor(pos / size * shape).astype(int)
ijk[:, :2] %= shape[:2]
ijk[:, 2] = np.clip(ijk[:, 2], 0, shape[2] - 1)
depth = depth_grid[ijk[:, 0], ijk[:, 1], ijk[:, 2]]
skeleton = types != 4
mask_at_atoms = mask[ijk[skeleton, 0], ijk[skeleton, 1], ijk[skeleton, 2]]

by_element = {}
by_depth = {}
for t, symbol in SYMBOL.items():
    selected = skeleton & (types == t)
    values = {}
    for lo, hi, label in BINS:
        chosen = selected & (depth >= lo) & (depth < hi)
        values[label] = stats(forces[chosen])
    by_depth[symbol] = values
    by_element[symbol] = stats(forces[selected])

# Record the highest-force backbone atoms and their local morphology-depth bins.
high = np.flatnonzero(skeleton)
high = high[np.argsort(forces[high])[-25:][::-1]]
top_atoms = [
    {
        "id": int(ids[i]), "element": SYMBOL[int(types[i])],
        "position_A": [float(x) for x in pos[i]],
        "voxel_boundary_depth_A": float(depth[i]),
        "force_norm_eV_A": float(forces[i]),
    }
    for i in high
]
result = {
    "status": "STATIC_RUN0_QSGS_MORPHOLOGY_CORRELATION_ONLY",
    "atoms": int(len(rows)),
    "skeleton_atoms": int(np.count_nonzero(skeleton)),
    "voxel_mask_match_fraction_for_skeleton_atoms": float(np.mean(mask_at_atoms)),
    "voxel_grid_shape": [int(x) for x in shape],
    "voxel_size_A": [float(x) for x in size / shape],
    "distance_definition": "Euclidean distance-transform depth inside the cleaned QSGS skeleton voxel mask; x/y tiled periodically, z open.",
    "all_skeleton_force": stats(forces[skeleton]),
    "backbone_by_element": by_element,
    "backbone_by_voxel_boundary_depth_A": by_depth,
    "top_25_backbone_force_atoms": top_atoms,
    "limitations": [
        "Voxel depth is a morphology proxy, not an atomistic bond coordination or interface distance.",
        "The force includes TSCBOP internal and Ag-Morse cross contributions; this analysis does not decompose them.",
        "The snapshot is unrelaxed, static, and conditional on the current potential files.",
    ],
}
OUT.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({
    "mask_match_fraction": result["voxel_mask_match_fraction_for_skeleton_atoms"],
    "all_skeleton": result["all_skeleton_force"],
    "by_depth": by_depth,
    "output": str(OUT),
}, indent=2))
