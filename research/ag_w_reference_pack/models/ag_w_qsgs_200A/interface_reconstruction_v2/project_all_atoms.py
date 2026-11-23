#!/usr/bin/env python3
"""Project all original lattice atoms out of short contacts without deleting sites.

This is a geometry-only reconstruction preconditioner. It repeatedly projects
Ag-W pairs below the EFS-informed screening distance apart, and prevents any
same-phase pair from becoming shorter than 2.60 Å. All original atoms and
species are retained; the QSGS mask and bicontinuous topology are unchanged.
The result is only a candidate for unconstrained Dai EFS static relaxation.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


PROJECT = Path(__file__).resolve().parents[3]
MODEL = PROJECT / "models" / "ag_w_qsgs_200A"
SOURCE = MODEL / "ag_w_qsgs.data.gz"
OUT = Path(__file__).resolve().parent
BOX_XY = 200.0
Z_SHIFT = 500.0
BOX = [200.0, 200.0, 1200.0]
CROSS_MIN = 2.8
SAME_MIN = 2.6
RELAX = 0.8
MAX_ATOM_STEP = 0.12
MAX_ITER = 160


def read_atoms(path):
    with gzip.open(path, "rt") as stream:
        lines = stream.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("Atoms")) + 1
    rows = []
    for line in lines[start:]:
        fields = line.split()
        if len(fields) < 5:
            if rows:
                break
            continue
        rows.append((int(fields[0]), int(fields[1]), *map(float, fields[2:5])))
    return np.asarray(rows, dtype=float)


def tree_positions(pos):
    p = pos.copy()
    p[:, 0] %= BOX_XY
    p[:, 1] %= BOX_XY
    p[:, 2] += Z_SHIFT
    return p


def add_pair_corrections(pos, left, right, target, corrections, weights):
    if len(left) == 0:
        return 0
    delta = pos[right] - pos[left]
    delta[:, 0] = (delta[:, 0] + BOX_XY / 2) % BOX_XY - BOX_XY / 2
    delta[:, 1] = (delta[:, 1] + BOX_XY / 2) % BOX_XY - BOX_XY / 2
    distance = np.linalg.norm(delta, axis=1)
    violation = distance < target
    if not np.any(violation):
        return 0
    left = left[violation]
    right = right[violation]
    delta = delta[violation]
    distance = distance[violation]
    # A deterministic direction handles exact coincidences without random jitter.
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


def measure(pos, types):
    p = tree_positions(pos)
    w = np.flatnonzero(types == 1)
    ag = np.flatnonzero(types == 2)
    wt = cKDTree(p[w], boxsize=BOX)
    at = cKDTree(p[ag], boxsize=BOX)
    cross = wt.sparse_distance_matrix(at, CROSS_MIN, output_type="coo_matrix")
    same_w = wt.query_pairs(SAME_MIN, output_type="ndarray")
    same_ag = at.query_pairs(SAME_MIN, output_type="ndarray")
    min_cross = float(min(wt.query(p[ag], workers=-1)[0].min(), at.query(p[w], workers=-1)[0].min()))
    min_w = float(wt.query(p[w], k=2, workers=-1)[0][:, 1].min())
    min_ag = float(at.query(p[ag], k=2, workers=-1)[0][:, 1].min())
    return w, ag, cross, same_w, same_ag, min_cross, min_w, min_ag


def write_data(path, atoms):
    with gzip.open(path, "wt", compresslevel=3) as stream:
        stream.write("Ag-W QSGS all-atom projected geometry candidate; geometry only\n\n")
        stream.write(f"{len(atoms)} atoms\n\n2 atom types\n\n")
        stream.write("0.0 200.0 xlo xhi\n0.0 200.0 ylo yhi\n0.0 300.0 zlo zhi\n\n")
        stream.write("Masses\n\n1 183.83999995270821 # W\n2 107.86819997225152 # Ag\n\n")
        stream.write("Atoms # atomic\n\n")
        for atom_id, typ, x, y, z in atoms:
            stream.write(f"{int(atom_id)} {int(typ)} {x:.15g} {y:.15g} {z:.15g}\n")


def main():
    atoms = read_atoms(SOURCE)
    original = atoms[:, 2:5].copy()
    pos = original.copy()
    types = atoms[:, 1].astype(np.int8)
    history = []

    for iteration in range(MAX_ITER + 1):
        w, ag, cross, same_w, same_ag, min_cross, min_w, min_ag = measure(pos, types)
        nviol = cross.nnz + len(same_w) + len(same_ag)
        history.append({
            "iteration": iteration,
            "cross_pairs_below_2_8_A": int(cross.nnz),
            "same_W_pairs_below_2_6_A": int(len(same_w)),
            "same_Ag_pairs_below_2_6_A": int(len(same_ag)),
            "minimum_cross_A": min_cross,
            "minimum_W_W_A": min_w,
            "minimum_Ag_Ag_A": min_ag,
        })
        if nviol == 0 or iteration == MAX_ITER:
            break

        corrections = np.zeros_like(pos)
        weights = np.zeros(len(pos), dtype=float)
        cross_left, cross_right = w[cross.row], ag[cross.col]
        n_cross = add_pair_corrections(pos, cross_left, cross_right, CROSS_MIN, corrections, weights)
        n_same_w = add_pair_corrections(pos, w[same_w[:, 0]], w[same_w[:, 1]], SAME_MIN, corrections, weights)
        n_same_ag = add_pair_corrections(pos, ag[same_ag[:, 0]], ag[same_ag[:, 1]], SAME_MIN, corrections, weights)
        touched = weights > 0
        delta = np.zeros_like(pos)
        delta[touched] = RELAX * corrections[touched] / weights[touched, None]
        norm = np.linalg.norm(delta, axis=1)
        scale = np.minimum(1.0, MAX_ATOM_STEP / np.maximum(norm, 1.0e-30))
        delta *= scale[:, None]
        pos += delta
        pos[:, 0] %= BOX_XY
        pos[:, 1] %= BOX_XY
        # Keep the starting material slab inside its specified free-surface region.
        pos[:, 2] = np.clip(pos[:, 2], 0.0, 200.0)
        history[-1]["projected_cross_pairs"] = n_cross
        history[-1]["projected_same_W_pairs"] = n_same_w
        history[-1]["projected_same_Ag_pairs"] = n_same_ag
        history[-1]["maximum_atom_step_A"] = float(norm.max())

    w, ag, cross, same_w, same_ag, min_cross, min_w, min_ag = measure(pos, types)
    conflict_free = cross.nnz == 0 and len(same_w) == 0 and len(same_ag) == 0
    displacement = pos - original
    displacement[:, :2] = (displacement[:, :2] + BOX_XY / 2) % BOX_XY - BOX_XY / 2
    displacement_norm = np.linalg.norm(displacement, axis=1)
    ratio = int(np.count_nonzero(types == 2)) * 107.86819997225152 / (
        int(np.count_nonzero(types == 1)) * 183.83999995270821
    )
    report = {
        "status": "GEOMETRY_CANDIDATE_ONLY_REQUIRES_DAI_EFS_STATIC_RELAXATION" if conflict_free else "REJECTED_PRECONDITION_NOT_CONFLICT_FREE",
        "method": "iterative pair-distance projection; no atoms deleted or reassigned",
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "qsgs_mask": "unchanged; all source atoms and species retained",
        "screening_distances_A": {"Ag_W": CROSS_MIN, "same_phase_floor": SAME_MIN},
        "projection_parameters": {"relaxation_factor": RELAX, "maximum_atom_step_A": MAX_ATOM_STEP, "maximum_iterations": MAX_ITER},
        "projection_history": history,
        "atom_counts": {"W": int(np.count_nonzero(types == 1)), "Ag": int(np.count_nonzero(types == 2))},
        "total_atoms": int(len(atoms)),
        "mass_ratio_Ag_over_W": ratio,
        "mass_ratio_relative_error": abs(ratio - 1.0),
        "minimum_distances_A": {"Ag_W": min_cross, "W_W": min_w, "Ag_Ag": min_ag},
        "maximum_displacement_A": float(displacement_norm.max()),
        "displacement_percentiles_A": {str(q): float(np.percentile(displacement_norm, q)) for q in (50, 90, 95, 99)},
        "atoms_displaced_over_0_5_A": int(np.count_nonzero(displacement_norm > 0.5)),
        "atoms_displaced_over_1_A": int(np.count_nonzero(displacement_norm > 1.0)),
        "static_dynamics": "NOT_RUN",
        "decision": "PENDING_DAI_EFS_STATIC_RELAXATION" if conflict_free else "REJECTED",
    }
    data_path = OUT / "ag_w_qsgs_interface_v2.projected.data.gz"
    write_data(data_path, np.column_stack((atoms[:, 0], types, pos)))
    report["candidate_file"] = data_path.name
    report["candidate_sha256"] = hashlib.sha256(data_path.read_bytes()).hexdigest()
    (OUT / "projection_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
