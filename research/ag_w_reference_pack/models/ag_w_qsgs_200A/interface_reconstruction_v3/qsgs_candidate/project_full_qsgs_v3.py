#!/usr/bin/env python3
"""Geometry-only all-atom Ag/W QSGS projection using v3 local-slab distance bounds."""
from __future__ import annotations

import gzip
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[4]
MODEL = ROOT / "models/ag_w_qsgs_200A"
SOURCE = MODEL / "ag_w_qsgs.data.gz"
OUT = Path(__file__).resolve().parent
BOX_XY = 200.0
MATERIAL_TOP_A = 200.0
TREE_BOX = [200.0, 200.0, 1200.0]
Z_SHIFT = 500.0

# Geometric screening floors are taken from the extrema of the v3 relaxed local slabs.
# They are preconditioner targets only; unconstrained EFS remains the final force check.
CROSS_MIN_A = 2.54
SAME_MIN_A = {1: 2.60, 2: 2.48}  # type 1 W, type 2 Ag
RELAXATION_FACTOR = 0.8
MAX_STEP_A = 0.15
MAX_ITERATIONS = 500
REPORT_EVERY = 10


def read_atoms(path: Path) -> np.ndarray:
    with gzip.open(path, "rt") as stream:
        lines = stream.readlines()
    section = next(i for i, line in enumerate(lines) if line.startswith("Atoms")) + 1
    rows = []
    for line in lines[section:]:
        fields = line.split()
        if len(fields) < 5:
            if rows:
                break
            continue
        rows.append((int(fields[0]), int(fields[1]), *map(float, fields[2:5])))
    return np.asarray(rows, dtype=float)


def tree_coords(pos: np.ndarray) -> np.ndarray:
    result = pos.copy()
    result[:, 0] %= BOX_XY
    result[:, 1] %= BOX_XY
    result[:, 2] += Z_SHIFT
    return result


def distance_and_pairs(pos: np.ndarray, types: np.ndarray):
    p = tree_coords(pos)
    w = np.flatnonzero(types == 1)
    ag = np.flatnonzero(types == 2)
    wt = cKDTree(p[w], boxsize=TREE_BOX)
    at = cKDTree(p[ag], boxsize=TREE_BOX)
    cross = wt.sparse_distance_matrix(at, CROSS_MIN_A, output_type="coo_matrix")
    keep = cross.data < CROSS_MIN_A - 1.0e-10
    cross_w = w[cross.row[keep]]
    cross_ag = ag[cross.col[keep]]
    same_w_local = wt.query_pairs(SAME_MIN_A[1], output_type="ndarray")
    same_ag_local = at.query_pairs(SAME_MIN_A[2], output_type="ndarray")
    same_w = np.column_stack((w[same_w_local[:, 0]], w[same_w_local[:, 1]])) if len(same_w_local) else np.empty((0, 2), int)
    same_ag = np.column_stack((ag[same_ag_local[:, 0]], ag[same_ag_local[:, 1]])) if len(same_ag_local) else np.empty((0, 2), int)
    min_cross = float(min(wt.query(p[ag], workers=-1)[0].min(), at.query(p[w], workers=-1)[0].min()))
    min_w = float(wt.query(p[w], k=2, workers=-1)[0][:, 1].min())
    min_ag = float(at.query(p[ag], k=2, workers=-1)[0][:, 1].min())
    return w, ag, (cross_w, cross_ag), same_w, same_ag, min_cross, min_w, min_ag


def add_corrections(
    pos: np.ndarray,
    left: np.ndarray,
    right: np.ndarray,
    target: float,
    corrections: np.ndarray,
    weights: np.ndarray,
) -> int:
    if len(left) == 0:
        return 0
    delta = pos[right] - pos[left]
    delta[:, 0] = (delta[:, 0] + BOX_XY / 2.0) % BOX_XY - BOX_XY / 2.0
    delta[:, 1] = (delta[:, 1] + BOX_XY / 2.0) % BOX_XY - BOX_XY / 2.0
    distance = np.linalg.norm(delta, axis=1)
    violation = distance < target - 1.0e-10
    if not np.any(violation):
        return 0
    left = left[violation]
    right = right[violation]
    delta = delta[violation]
    distance = distance[violation]
    unit = np.zeros_like(delta)
    nonzero = distance > 1.0e-12
    unit[nonzero] = delta[nonzero] / distance[nonzero, None]
    if np.any(~nonzero):
        axis = (left[~nonzero] + right[~nonzero]) % 3
        unit[np.flatnonzero(~nonzero), axis] = 1.0
    correction = 0.5 * (target - distance)[:, None] * unit
    np.add.at(corrections, left, -correction)
    np.add.at(corrections, right, correction)
    np.add.at(weights, left, 1.0)
    np.add.at(weights, right, 1.0)
    return int(len(left))


def topology_summary(pos: np.ndarray, types: np.ndarray) -> dict:
    sys.path.insert(0, str(ROOT / "scripts"))
    import build_ag_w_qsgs as builder

    w = np.flatnonzero(types == 1)
    ag = np.flatnonzero(types == 2)
    return {
        "W": builder.topology(pos[w], "W"),
        "Ag": builder.topology(pos[ag], "Ag"),
    }


def write_lammps_data(path: Path, atoms: np.ndarray, types: np.ndarray, pos: np.ndarray) -> None:
    with gzip.open(path, "wt", compresslevel=6) as stream:
        stream.write("Ag-W QSGS v3 geometry-only candidate from collision-aware projection\n\n")
        stream.write(f"{len(atoms)} atoms\n\n2 atom types\n\n")
        stream.write("0.0 200.0 xlo xhi\n0.0 200.0 ylo yhi\n0.0 300.0 zlo zhi\n\n")
        stream.write("Masses\n\n1 183.83999995270821 # W\n2 107.86819997225152 # Ag\n\n")
        stream.write("Atoms # atomic\n\n")
        for atom_id, atom_type, xyz in zip(atoms[:, 0].astype(int), types, pos):
            stream.write(f"{atom_id} {int(atom_type)} {xyz[0]:.15g} {xyz[1]:.15g} {xyz[2]:.15g}\n")


def main() -> None:
    atoms = read_atoms(SOURCE)
    original = atoms[:, 2:5].copy()
    pos = original.copy()
    types = atoms[:, 1].astype(np.int8)
    history = []
    for iteration in range(MAX_ITERATIONS + 1):
        w, ag, cross, same_w, same_ag, min_cross, min_w, min_ag = distance_and_pairs(pos, types)
        nviol = len(cross[0]) + len(same_w) + len(same_ag)
        entry = {
            "iteration": iteration,
            "Ag_W_pairs_below_2_54_A": int(len(cross[0])),
            "W_W_pairs_below_2_60_A": int(len(same_w)),
            "Ag_Ag_pairs_below_2_48_A": int(len(same_ag)),
            "minimum_distances_A": {"Ag_W": min_cross, "W_W": min_w, "Ag_Ag": min_ag},
        }
        if iteration == 0 or iteration % REPORT_EVERY == 0 or nviol == 0 or iteration == MAX_ITERATIONS:
            history.append(entry)
            print(json.dumps(entry), flush=True)
        if nviol == 0 or iteration == MAX_ITERATIONS:
            break

        corrections = np.zeros_like(pos)
        weights = np.zeros(len(pos), dtype=float)
        n_cross = add_corrections(pos, cross[0], cross[1], CROSS_MIN_A, corrections, weights)
        n_same_w = add_corrections(pos, same_w[:, 0], same_w[:, 1], SAME_MIN_A[1], corrections, weights)
        n_same_ag = add_corrections(pos, same_ag[:, 0], same_ag[:, 1], SAME_MIN_A[2], corrections, weights)
        touched = weights > 0
        delta = np.zeros_like(pos)
        delta[touched] = RELAXATION_FACTOR * corrections[touched] / weights[touched, None]
        norms = np.linalg.norm(delta, axis=1)
        delta *= np.minimum(1.0, MAX_STEP_A / np.maximum(norms, 1.0e-30))[:, None]
        pos += delta
        pos[:, :2] %= BOX_XY
        pos[:, 2] = np.clip(pos[:, 2], 0.0, MATERIAL_TOP_A)
        if iteration % REPORT_EVERY == 0:
            history[-1]["projected_pairs"] = {"Ag_W": n_cross, "W_W": n_same_w, "Ag_Ag": n_same_ag}
            history[-1]["maximum_step_A"] = float(norms.max())

    w, ag, cross, same_w, same_ag, min_cross, min_w, min_ag = distance_and_pairs(pos, types)
    conflict_free = len(cross[0]) == len(same_w) == len(same_ag) == 0
    displacement = pos - original
    displacement[:, :2] = (displacement[:, :2] + BOX_XY / 2) % BOX_XY - BOX_XY / 2
    displacement_norm = np.linalg.norm(displacement, axis=1)
    topology = topology_summary(pos, types)
    topology_pass = all(
        phase["xy_winding_rank"] == 2
        and phase["largest_component_spans_z"]
        and phase["largest_component_share"] >= 0.98
        for phase in topology.values()
    )
    ratio = int(np.count_nonzero(types == 2)) * 107.86819997225152 / (
        int(np.count_nonzero(types == 1)) * 183.83999995270821
    )
    report = {
        "title": "Ag-W v3 exploratory full-size QSGS geometric candidate",
        "status": "GEOMETRY_CANDIDATE_REQUIRES_UNCONSTRAINED_EFS_STATIC_CHECKS" if conflict_free and topology_pass else "REJECTED_GEOMETRY_PRECHECK",
        "method": "simultaneous all-atom distance projection using minima observed in the v3 relaxed local slabs",
        "source_file": str(SOURCE.relative_to(ROOT)),
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "qsgs_voxel_mask": "unchanged; all atoms retained and no type changes",
        "box_A": [200.0, 200.0, 300.0],
        "material_z_A": [0.0, 200.0],
        "boundary": "p p f",
        "screening_distances_A": {"Ag_W": CROSS_MIN_A, "W_W": SAME_MIN_A[1], "Ag_Ag": SAME_MIN_A[2]},
        "screening_rationale": "geometric preconditioner floors only; the final gate is unconstrained Dai EFS static relaxation",
        "projection_parameters": {
            "relaxation_factor": RELAXATION_FACTOR,
            "maximum_atom_step_A": MAX_STEP_A,
            "maximum_iterations": MAX_ITERATIONS,
        },
        "history": history,
        "remaining_violations": {"Ag_W": int(len(cross[0])), "W_W": int(len(same_w)), "Ag_Ag": int(len(same_ag))},
        "atom_counts": {"W": int(len(w)), "Ag": int(len(ag)), "total": int(len(atoms))},
        "mass_ratio_Ag_over_W": ratio,
        "mass_ratio_relative_error": abs(ratio - 1.0),
        "minimum_distances_A": {"Ag_W": min_cross, "W_W": min_w, "Ag_Ag": min_ag},
        "maximum_displacement_A": float(displacement_norm.max()),
        "displacement_percentiles_A": {str(q): float(np.percentile(displacement_norm, q)) for q in (50, 90, 95, 99)},
        "atoms_displaced_over_0_5_A": int(np.count_nonzero(displacement_norm > 0.5)),
        "atoms_displaced_over_1_A": int(np.count_nonzero(displacement_norm > 1.0)),
        "topology": topology,
        "both_phases_bicontinuous": topology_pass,
        "static_dynamics": "NOT_RUN",
        "interpretation": "candidate geometry only; not an MD initial state until the EFS static report passes",
    }
    data_path = OUT / "ag_w_qsgs_v3_exploratory_geometry.data.gz"
    write_lammps_data(data_path, atoms, types, pos)
    report["candidate_file"] = data_path.name
    report["candidate_sha256"] = hashlib.sha256(data_path.read_bytes()).hexdigest()
    (OUT / "geometry_projection_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "history"}, indent=2))


if __name__ == "__main__":
    main()
