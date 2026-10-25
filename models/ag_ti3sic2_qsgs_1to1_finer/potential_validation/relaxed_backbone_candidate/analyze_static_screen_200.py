#!/usr/bin/env python3
"""Analyze force/displacement and phase connectivity for the 100-step screen."""
from __future__ import annotations
import gzip
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

HERE = Path(__file__).resolve().parent
MODEL = HERE.parents[1]
BOX = np.asarray(json.loads((HERE / "geometry_report.json").read_text())["simulation_box_A"], dtype=float)
ELEMENT = {1: "Ti", 2: "Si", 3: "C", 4: "Ag"}


def load_dump(path: Path) -> np.ndarray:
    with gzip.open(path, "rt") as stream:
        return np.loadtxt(stream, skiprows=9)


def stats(v: np.ndarray) -> dict:
    if not len(v):
        return {"n": 0}
    return {
        "n": int(len(v)),
        "p50_eV_A": float(np.percentile(v, 50)),
        "p95_eV_A": float(np.percentile(v, 95)),
        "p99_eV_A": float(np.percentile(v, 99)),
        "max_eV_A": float(np.max(v)),
        "n_gt_2_eV_A": int(np.count_nonzero(v > 2.0)),
        "fraction_gt_2": float(np.mean(v > 2.0)),
    }


def phase_topology(points: np.ndarray, ids: np.ndarray, labels: np.ndarray, cutoff: float, box: list[float]) -> dict:
    pairs = cKDTree(points, boxsize=box).query_pairs(cutoff, output_type="ndarray")
    graph = coo_matrix(
        (np.ones(len(pairs) * 2, dtype=np.uint8),
         (np.r_[pairs[:, 0], pairs[:, 1]], np.r_[pairs[:, 1], pairs[:, 0]])),
        shape=(len(points), len(points)),
    ).tocsr()
    ncomp, comp = connected_components(graph, directed=False)
    sizes = np.bincount(comp)
    largest = int(sizes.max())
    small_ids = np.flatnonzero(sizes < largest)
    outside = np.flatnonzero(np.isin(comp, small_ids))
    near_d, near_i = cKDTree(points, boxsize=box).query(points[outside], k=2) if len(outside) else (np.empty((0, 2)), np.empty((0, 2), int))
    detached = [
        {"id": int(ids[i]), "element": str(labels[i]), "nearest_distance_A": float(near_d[j, 1]), "nearest_atom_id": int(ids[near_i[j, 1]])}
        for j, i in enumerate(outside)
    ]
    return {"components": int(ncomp), "largest_component_atoms": largest, "atoms_outside_largest": int(len(points) - largest), "cutoff_A": cutoff, "detached_atoms": detached}


initial = load_dump(HERE / "initial_forces.dump.gz")
final = load_dump(HERE / "mobile_minimized_200.dump.gz")
if not (np.array_equal(initial[:, 0], final[:, 0]) and np.array_equal(initial[:, 1], final[:, 1])):
    raise RuntimeError("Initial and endpoint atom IDs/types differ")
ids = initial[:, 0].astype(np.int64)
types = initial[:, 1].astype(np.int8)
init_pos = initial[:, 2:5]
final_pos = final[:, 2:5]
delta = final_pos - init_pos
delta[:, 0] -= BOX[0] * np.rint(delta[:, 0] / BOX[0])
delta[:, 1] -= BOX[1] * np.rint(delta[:, 1] / BOX[1])
displacement = np.linalg.norm(delta, axis=1)
force = np.linalg.norm(final[:, 5:8], axis=1)
mobile = init_pos[:, 2] >= 20.0
mobile_skeleton = mobile & (types != 4)
mobile_ag = mobile & (types == 4)
by_element = {}
for t, element in ELEMENT.items():
    selected = mobile & (types == t)
    by_element[element] = {"displacement_A": {
        "p50": float(np.percentile(displacement[selected], 50)),
        "p95": float(np.percentile(displacement[selected], 95)),
        "max": float(np.max(displacement[selected])),
    }, "force": stats(force[selected])}

xyz = final_pos.copy()
xyz[:, 0] %= BOX[0]
xyz[:, 1] %= BOX[1]
periodic_box = [float(BOX[0]), float(BOX[1]), 1000.0]
select = {name: types == t for t, name in ELEMENT.items()}
trees = {name: cKDTree(xyz[mask], boxsize=periodic_box) for name, mask in select.items()}
min_pair = {}
for a, b in (("Ag", "Ti"), ("Ag", "Si"), ("Ag", "C"), ("Ti", "C")):
    d, _ = trees[b].query(xyz[select[a]], k=1)
    min_pair[f"{a}-{b}"] = float(np.min(d))

symbols = np.array([ELEMENT[int(t)] for t in types])
skeleton_topology = phase_topology(xyz[types != 4], ids[types != 4], symbols[types != 4], 3.25, periodic_box)
ag_topology = phase_topology(xyz[types == 4], ids[types == 4], symbols[types == 4], 3.15, periodic_box)
idx_sk = np.flatnonzero(mobile_skeleton)[np.argmax(force[mobile_skeleton])]
idx_ag = np.flatnonzero(mobile_ag)[np.argmax(force[mobile_ag])]
result = {
    "status": "STATIC_200_FIRE_ITERATION_SCREEN_NOT_A_COLD_BASELINE",
    "atoms": int(len(ids)),
    "max_mobile_force_eV_A": float(np.max(force[mobile])),
    "max_mobile_skeleton_force_eV_A": float(np.max(force[mobile_skeleton])),
    "max_mobile_Ag_force_eV_A": float(np.max(force[mobile_ag])),
    "max_mobile_skeleton_displacement_A": float(np.max(displacement[mobile_skeleton])),
    "by_mobile_element": by_element,
    "minimum_pair_distances_A": min_pair,
    "phase_connectivity": {"Ti3SiC2": skeleton_topology, "Ag": ag_topology},
    "largest_force_atoms": {
        "skeleton": {"id": int(ids[idx_sk]), "element": ELEMENT[int(types[idx_sk])], "force_eV_A": float(force[idx_sk]), "position_A": [float(x) for x in final_pos[idx_sk]], "displacement_A": float(displacement[idx_sk])},
        "Ag": {"id": int(ids[idx_ag]), "force_eV_A": float(force[idx_ag]), "position_A": [float(x) for x in final_pos[idx_ag]], "displacement_A": float(displacement[idx_ag])},
    },
    "limitations": [
        "The two 100-step FIRE minimization segments hit their iteration caps.",
        "Only the 0-20 A bottom slab was fixed; no dynamical or finite-temperature stability was tested.",
        "Current calibration-v0 cross potential and TSCBOP are not production-validated for the intended hot composite.",
    ],
}
(HERE / "static_screen_metrics_200.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
