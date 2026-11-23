#!/usr/bin/env python3
"""Sequentially project every atom-pair overlap while retaining all atoms."""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from project_all_atoms import BOX, BOX_XY, CROSS_MIN, OUT, SOURCE, Z_SHIFT, measure, read_atoms


SAME_MIN = 2.6
MAX_SWEEPS = 80
RELAX = 1.0
SEED = 20261123


def pair_lists(pos, types):
    p = pos.copy()
    p[:, 0] %= BOX_XY
    p[:, 1] %= BOX_XY
    p[:, 2] += Z_SHIFT
    w = np.flatnonzero(types == 1)
    ag = np.flatnonzero(types == 2)
    wt = cKDTree(p[w], boxsize=BOX)
    at = cKDTree(p[ag], boxsize=BOX)
    cross = wt.sparse_distance_matrix(at, CROSS_MIN, output_type="coo_matrix")
    cross_keep = cross.data < CROSS_MIN - 1.0e-8
    cross_row, cross_col = cross.row[cross_keep], cross.col[cross_keep]
    same_w = wt.query_pairs(SAME_MIN, output_type="ndarray")
    same_ag = at.query_pairs(SAME_MIN, output_type="ndarray")
    if len(same_w):
        dw = p[w[same_w[:, 1]]] - p[w[same_w[:, 0]]]
        dw[:, :2] = (dw[:, :2] + BOX_XY / 2) % BOX_XY - BOX_XY / 2
        same_w = same_w[np.linalg.norm(dw, axis=1) < SAME_MIN - 1.0e-8]
    if len(same_ag):
        da = p[ag[same_ag[:, 1]]] - p[ag[same_ag[:, 0]]]
        da[:, :2] = (da[:, :2] + BOX_XY / 2) % BOX_XY - BOX_XY / 2
        same_ag = same_ag[np.linalg.norm(da, axis=1) < SAME_MIN - 1.0e-8]
    edges_left = np.concatenate((w[cross_row], w[same_w[:, 0]], ag[same_ag[:, 0]]))
    edges_right = np.concatenate((ag[cross_col], w[same_w[:, 1]], ag[same_ag[:, 1]]))
    targets = np.concatenate((np.full(len(cross_row), CROSS_MIN), np.full(len(same_w) + len(same_ag), SAME_MIN)))
    return edges_left, edges_right, targets, int(len(cross_row)), int(len(same_w)), int(len(same_ag))


def one_sweep(pos, left, right, targets, rng, reverse):
    order = rng.permutation(len(left))
    if reverse:
        order = order[::-1]
    max_gap = 0.0
    for edge in order:
        i, j = int(left[edge]), int(right[edge])
        delta = pos[j] - pos[i]
        delta[0] = (delta[0] + BOX_XY / 2) % BOX_XY - BOX_XY / 2
        delta[1] = (delta[1] + BOX_XY / 2) % BOX_XY - BOX_XY / 2
        distance = float(np.linalg.norm(delta))
        gap = float(targets[edge] - distance)
        if gap <= 1.0e-10:
            continue
        max_gap = max(max_gap, gap)
        if distance > 1.0e-12:
            direction = delta / distance
        else:
            direction = np.zeros(3)
            direction[(i + j) % 3] = 1.0
        step = 0.5 * RELAX * gap * direction
        pos[i] -= step
        pos[j] += step
        pos[i, :2] %= BOX_XY
        pos[j, :2] %= BOX_XY
        pos[i, 2] = np.clip(pos[i, 2], 0.0, 200.0)
        pos[j, 2] = np.clip(pos[j, 2], 0.0, 200.0)
    return max_gap


def write_data(path, atoms):
    with gzip.open(path, "wt", compresslevel=3) as stream:
        stream.write("Ag-W QSGS sequential all-atom projection candidate; geometry only\n\n")
        stream.write(f"{len(atoms)} atoms\n\n2 atom types\n\n")
        stream.write("0.0 200.0 xlo xhi\n0.0 200.0 ylo yhi\n0.0 300.0 zlo zhi\n\n")
        stream.write("Masses\n\n1 183.83999995270821 # W\n2 107.86819997225152 # Ag\n\nAtoms # atomic\n\n")
        for atom_id, typ, x, y, z in atoms:
            stream.write(f"{int(atom_id)} {int(typ)} {x:.15g} {y:.15g} {z:.15g}\n")


def main():
    source = read_atoms(SOURCE)
    original = source[:, 2:5].copy()
    pos = original.copy()
    types = source[:, 1].astype(np.int8)
    rng = np.random.default_rng(SEED)
    history = []
    for sweep in range(MAX_SWEEPS + 1):
        left, right, targets, nc, nw, na = pair_lists(pos, types)
        if len(left):
            delta = pos[right] - pos[left]
            delta[:, 0] = (delta[:, 0] + BOX_XY / 2) % BOX_XY - BOX_XY / 2
            delta[:, 1] = (delta[:, 1] + BOX_XY / 2) % BOX_XY - BOX_XY / 2
            min_gap = float(np.max(targets - np.linalg.norm(delta, axis=1)))
        else:
            min_gap = 0.0
        entry = {"sweep": sweep, "cross_pairs_below_2_8_A": nc, "same_W_pairs_below_2_6_A": nw,
                 "same_Ag_pairs_below_2_6_A": na, "maximum_current_penetration_A": min_gap}
        history.append(entry)
        print(json.dumps(entry), flush=True)
        if len(left) == 0 or sweep == MAX_SWEEPS:
            break
        entry["maximum_pair_penetration_before_sweep_A"] = one_sweep(pos, left, right, targets, rng, bool(sweep % 2))

    left, right, targets, nc, nw, na = pair_lists(pos, types)
    displacement = pos - original
    displacement[:, :2] = (displacement[:, :2] + BOX_XY / 2) % BOX_XY - BOX_XY / 2
    dn = np.linalg.norm(displacement, axis=1)
    wi, agi, _, _, _, min_cross, min_w, min_ag = measure(pos, types)
    import sys
    sys.path.insert(0, str(SOURCE.parents[2] / "scripts"))
    import build_ag_w_qsgs as builder
    topology = {
        "W": builder.topology(pos[wi], "W"),
        "Ag": builder.topology(pos[agi], "Ag"),
    }
    topology_pass = all(
        item["xy_winding_rank"] == 2 and item["largest_component_spans_z"]
        and item["largest_component_share"] >= 0.98 for item in topology.values()
    )
    geometry_pass = len(left) == 0 and topology_pass
    report = {
        "status": "GEOMETRY_CANDIDATE_ONLY_REQUIRES_DAI_EFS_STATIC_RELAXATION" if geometry_pass else "REJECTED_PROJECTION_OR_TOPOLOGY",
        "method": "sequential pair-distance projection with alternating pair order; all atoms retained",
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "qsgs_mask": "unchanged; original species counts retained",
        "screening_distances_A": {"Ag_W": CROSS_MIN, "same_phase_floor": SAME_MIN},
        "history": history,
        "remaining_violations": {"Ag_W_pairs": nc, "W_W_pairs": nw, "Ag_Ag_pairs": na},
        "minimum_distances_A": {"Ag_W": min_cross, "W_W": min_w, "Ag_Ag": min_ag},
        "atom_counts": {"W": int(np.count_nonzero(types == 1)), "Ag": int(np.count_nonzero(types == 2))},
        "mass_ratio_Ag_over_W": int(np.count_nonzero(types == 2)) * 107.86819997225152 / (int(np.count_nonzero(types == 1)) * 183.83999995270821),
        "minimum_atom_displacement_A": {"median": float(np.median(dn)), "p95": float(np.percentile(dn, 95)), "p99": float(np.percentile(dn, 99)), "max": float(dn.max())},
        "static_dynamics": "NOT_RUN",
        "topology": topology,
    }
    out = OUT / "ag_w_qsgs_interface_v2.sequential.data.gz"
    write_data(out, np.column_stack((source[:, 0], types, pos)))
    report["candidate_file"] = out.name
    report["candidate_sha256"] = hashlib.sha256(out.read_bytes()).hexdigest()
    (OUT / "sequential_projection_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "history"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
