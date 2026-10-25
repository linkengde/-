#!/usr/bin/env python3
"""Analyze displacement, released forces, and contact minima for the Ag-fixed screen."""
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
GEO = json.loads((MODEL / "geometry_report.json").read_text())
BOX = np.asarray(GEO["simulation_box_A"], dtype=float)
ELEMENT = {1: "Ti", 2: "Si", 3: "C", 4: "Ag"}


def load_dump(path: Path) -> np.ndarray:
    with gzip.open(path, "rt") as stream:
        return np.loadtxt(stream, skiprows=9)


def force_stats(force: np.ndarray) -> dict:
    if not len(force):
        return {"n": 0}
    return {
        "n": int(len(force)),
        "p50_eV_A": float(np.percentile(force, 50)),
        "p95_eV_A": float(np.percentile(force, 95)),
        "max_eV_A": float(np.max(force)),
        "n_gt_2_eV_A": int(np.count_nonzero(force > 2.0)),
    }


initial = load_dump(HERE / "initial_forces.dump.gz")
final = load_dump(HERE / "fixed_ag_skeleton_minimized.dump.gz")
if not (np.array_equal(initial[:, 0], final[:, 0]) and np.array_equal(initial[:, 1], final[:, 1])):
    raise RuntimeError("Atom IDs/types differ between initial and final snapshots")
ids = initial[:, 0].astype(np.int64)
types = initial[:, 1].astype(np.int8)
initial_pos = initial[:, 2:5]
final_pos = final[:, 2:5]
delta = final_pos - initial_pos
delta[:, 0] -= BOX[0] * np.rint(delta[:, 0] / BOX[0])
delta[:, 1] -= BOX[1] * np.rint(delta[:, 1] / BOX[1])
displacement = np.linalg.norm(delta, axis=1)
final_force = np.linalg.norm(final[:, 5:8], axis=1)
# The original static input fixes atoms whose starting z lies in [0, 20) A.
mobile_skeleton = (types != 4) & (initial_pos[:, 2] >= 20.0)
final_xyz = final_pos.copy()
final_xyz[:, 0] %= BOX[0]
final_xyz[:, 1] %= BOX[1]
periodic_box = [float(BOX[0]), float(BOX[1]), 1000.0]  # z is effectively nonperiodic
sel = {name: types == t for t, name in ELEMENT.items()}
trees = {name: cKDTree(final_xyz[mask], boxsize=periodic_box) for name, mask in sel.items()}
minimum_distances = {}
for a, b in (("Ag", "Ti"), ("Ag", "Si"), ("Ag", "C"), ("Ti", "C")):
    d, _ = trees[b].query(final_xyz[sel[a]], k=1, workers=4)
    minimum_distances[f"{a}-{b}"] = float(np.min(d))

# Recheck atomic backbone connectivity with the repository's 3.25 A criterion.
skel_xyz = final_xyz[types != 4]
skel_pairs = cKDTree(skel_xyz, boxsize=periodic_box).query_pairs(3.25, output_type="ndarray")
skel_graph = coo_matrix(
    (np.ones(len(skel_pairs) * 2, dtype=np.uint8),
     (np.r_[skel_pairs[:, 0], skel_pairs[:, 1]], np.r_[skel_pairs[:, 1], skel_pairs[:, 0]])),
    shape=(len(skel_xyz), len(skel_xyz)),
).tocsr()
n_components, labels = connected_components(skel_graph, directed=False)
component_sizes = np.bincount(labels)
largest_component = int(component_sizes.max())
small_components = np.flatnonzero(component_sizes < largest_component)
small_atom_indices = np.flatnonzero(np.isin(labels, small_components))
nearest_d, nearest_idx = cKDTree(skel_xyz, boxsize=periodic_box).query(
    skel_xyz[small_atom_indices], k=2, workers=4
)
detached_atoms = [
    {
        "id": int(ids[np.flatnonzero(types != 4)[i]]),
        "element": ELEMENT[int(types[np.flatnonzero(types != 4)[i]])],
        "nearest_skeleton_distance_A": float(nearest_d[j, 1]),
        "nearest_skeleton_atom_id": int(ids[np.flatnonzero(types != 4)[nearest_idx[j, 1]]]),
    }
    for j, i in enumerate(small_atom_indices)
]

by_element = {}
for t, name in ELEMENT.items():
    mask = (types == t) & (mobile_skeleton if t != 4 else np.ones(len(types), dtype=bool))
    by_element[name] = {
        "displacement_A": {
            "p50": float(np.percentile(displacement[mask], 50)),
            "p95": float(np.percentile(displacement[mask], 95)),
            "max": float(np.max(displacement[mask])),
        },
        "released_force_eV_A": force_stats(final_force[mask]),
    }
max_skel_idx = np.flatnonzero(mobile_skeleton)[np.argmax(final_force[mobile_skeleton])]
max_ag_idx = np.flatnonzero(types == 4)[np.argmax(final_force[types == 4])]
result = {
    "status": "STATIC_MINIMIZATION_SCREEN_STOPPED_AT_100_ITERATION_CAP_NOT_ACCEPTED_AS_COLD_BASELINE",
    "atoms": int(len(ids)),
    "box_A": [float(x) for x in BOX],
    "mobile_skeleton_atoms": int(np.count_nonzero(mobile_skeleton)),
    "maximum_skeleton_displacement_A": float(np.max(displacement[types != 4])),
    "maximum_skeleton_force_after_minimization_eV_A": float(np.max(final_force[mobile_skeleton])),
    "maximum_Ag_force_after_releasing_Ag_constraint_eV_A": float(np.max(final_force[types == 4])),
    "by_element": by_element,
    "minimum_pair_distances_after_minimization_A": minimum_distances,
    "skeleton_connectivity_at_3p25A": {
        "components": int(n_components),
        "largest_component_atoms": largest_component,
        "atoms_outside_largest_component": int(len(skel_xyz) - largest_component),
        "detached_atoms": detached_atoms,
    },
    "maximum_force_atoms": {
        "skeleton": {
            "id": int(ids[max_skel_idx]), "element": ELEMENT[int(types[max_skel_idx])],
            "force_eV_A": float(final_force[max_skel_idx]),
            "position_A": [float(x) for x in final_pos[max_skel_idx]],
            "displacement_A": float(displacement[max_skel_idx]),
        },
        "Ag": {
            "id": int(ids[max_ag_idx]), "force_eV_A": float(final_force[max_ag_idx]),
            "position_A": [float(x) for x in final_pos[max_ag_idx]],
        },
    },
    "limitations": [
        "Static 100-iteration FIRE screen only; minimization stopped at its iteration cap.",
        "Ag was fixed during skeleton relaxation; Ag forces are reported only after removing that constraint.",
        "Bottom atoms remained fixed; no MD, TTM, heat source, or topology revalidation was performed.",
        "A low force on most atoms does not establish structural or thermal stability.",
    ],
}
out = HERE / "fixed_ag_skeleton_minimize_metrics.json"
out.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
