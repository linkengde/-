#!/usr/bin/env python3
"""Build a collision-screened Ag-W candidate without moving lattice sites.

The original QSGS mask and phase-specific BCC/FCC sites are retained. A
bipartite conflict graph marks W-Ag site pairs closer than a screening radius.
A deterministic mass-balanced greedy vertex cover removes a small boundary
subset so no marked pair remains. This is a geometry construction step, not a
physical relaxation; the resulting candidate must pass Dai EFS static checks.
"""
from __future__ import annotations

import gzip
import heapq
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
R_SCREEN = 2.8  # Angstrom; geometric pre-screen only, not a force acceptance gate
BOX_XY = 200.0
Z_SHIFT = 500.0
BOX_TREE = [BOX_XY, BOX_XY, 1200.0]
MASS = {1: 183.83999995270821, 2: 107.86819997225152}


def read_atomic_data(path: Path) -> np.ndarray:
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


def top_vertex(heap, degree, active):
    while heap:
        neg_degree, vertex = heapq.heappop(heap)
        if active[vertex] and -neg_degree == degree[vertex] and degree[vertex] > 0:
            return vertex
    return None


def mass_balanced_vertex_cover(n_w, n_ag, wi, ai):
    """Greedy bipartite cover; balance fractional removals in both phases."""
    degree_w = np.bincount(wi, minlength=n_w).astype(np.int32)
    degree_ag = np.bincount(ai, minlength=n_ag).astype(np.int32)
    active_w = np.ones(n_w, dtype=bool)
    active_ag = np.ones(n_ag, dtype=bool)
    heap_w = [(-int(d), int(i)) for i, d in enumerate(degree_w) if d]
    heap_ag = [(-int(d), int(i)) for i, d in enumerate(degree_ag) if d]
    heapq.heapify(heap_w)
    heapq.heapify(heap_ag)

    w_order = np.argsort(wi, kind="stable")
    wi_w, ai_w = wi[w_order], ai[w_order]
    start_w = np.searchsorted(wi_w, np.arange(n_w + 1))
    ag_order = np.argsort(ai, kind="stable")
    ai_ag, wi_ag = ai[ag_order], wi[ag_order]
    start_ag = np.searchsorted(ai_ag, np.arange(n_ag + 1))

    remaining_edges = len(wi)
    removed_w = removed_ag = 0
    while remaining_edges:
        # Preserve the input composition by removing approximately the same
        # fraction of each phase, rather than equal numbers of unlike atoms.
        w_fraction = removed_w / n_w
        ag_fraction = removed_ag / n_ag
        prefer_w = w_fraction <= ag_fraction
        if prefer_w:
            vertex = top_vertex(heap_w, degree_w, active_w)
            if vertex is not None:
                active_w[vertex] = False
                n_edges = int(degree_w[vertex])
                degree_w[vertex] = 0
                remaining_edges -= n_edges
                removed_w += 1
                for other in ai_w[start_w[vertex] : start_w[vertex + 1]]:
                    if active_ag[other] and degree_ag[other] > 0:
                        degree_ag[other] -= 1
                        if degree_ag[other]:
                            heapq.heappush(heap_ag, (-int(degree_ag[other]), int(other)))
                continue
        vertex = top_vertex(heap_ag, degree_ag, active_ag)
        if vertex is not None:
            active_ag[vertex] = False
            n_edges = int(degree_ag[vertex])
            degree_ag[vertex] = 0
            remaining_edges -= n_edges
            removed_ag += 1
            for other in wi_ag[start_ag[vertex] : start_ag[vertex + 1]]:
                if active_w[other] and degree_w[other] > 0:
                    degree_w[other] -= 1
                    if degree_w[other]:
                        heapq.heappush(heap_w, (-int(degree_w[other]), int(other)))
            continue
        vertex = top_vertex(heap_w, degree_w, active_w)
        if vertex is None:
            raise RuntimeError(f"Conflict graph stalled with {remaining_edges} uncovered edges")
        active_w[vertex] = False
        n_edges = int(degree_w[vertex])
        degree_w[vertex] = 0
        remaining_edges -= n_edges
        removed_w += 1
        for other in ai_w[start_w[vertex] : start_w[vertex + 1]]:
            if active_ag[other] and degree_ag[other] > 0:
                degree_ag[other] -= 1
                if degree_ag[other]:
                    heapq.heappush(heap_ag, (-int(degree_ag[other]), int(other)))
    return active_w, active_ag, removed_w, removed_ag


def nearest_distance(source, target):
    src = source.copy()
    tgt = target.copy()
    src[:, 2] += Z_SHIFT
    tgt[:, 2] += Z_SHIFT
    tree = cKDTree(tgt, boxsize=BOX_TREE)
    return tree.query(src, k=1, workers=-1)[0]


def same_phase_nearest(points):
    tree_points = points.copy()
    tree_points[:, 2] += Z_SHIFT
    return cKDTree(tree_points, boxsize=BOX_TREE).query(
        tree_points, k=2, workers=-1
    )[0][:, 1]


def pair_force_magnitude(r):
    # Dai EFS Ag-W radial pair term: V=(r-c)^2 sum_k c_k r^k.
    c = 4.3
    coeff = (25.618, -26.3407, 9.1414, -1.0693, 0.0)
    poly = sum(value * r**k for k, value in enumerate(coeff))
    dpoly = sum(k * value * r ** (k - 1) for k, value in enumerate(coeff) if k)
    dV = 2.0 * (r - c) * poly + (r - c) ** 2 * dpoly
    return abs(dV)


def write_lammps(path, atoms):
    with gzip.open(path, "wt", compresslevel=3) as stream:
        stream.write("Ag-W QSGS collision-screened geometry candidate; geometry only\n\n")
        stream.write(f"{len(atoms)} atoms\n\n2 atom types\n\n")
        stream.write("0.0 200.0 xlo xhi\n0.0 200.0 ylo yhi\n0.0 300.0 zlo zhi\n\n")
        stream.write("Masses\n\n")
        stream.write(f"1 {MASS[1]:.14g} # W\n2 {MASS[2]:.14g} # Ag\n\n")
        stream.write("Atoms # atomic\n\n")
        for atom_id, atom_type, x, y, z in atoms:
            stream.write(f"{int(atom_id)} {int(atom_type)} {x:.15g} {y:.15g} {z:.15g}\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    source = read_atomic_data(SOURCE)
    w_rows = np.flatnonzero(source[:, 1] == 1)
    ag_rows = np.flatnonzero(source[:, 1] == 2)
    w_pos = source[w_rows, 2:5]
    ag_pos = source[ag_rows, 2:5]

    w_tree_pos = w_pos.copy()
    ag_tree_pos = ag_pos.copy()
    w_tree_pos[:, 2] += Z_SHIFT
    ag_tree_pos[:, 2] += Z_SHIFT
    wt = cKDTree(w_tree_pos, boxsize=BOX_TREE)
    at = cKDTree(ag_tree_pos, boxsize=BOX_TREE)
    conflicts = wt.sparse_distance_matrix(
        at, R_SCREEN, output_type="coo_matrix"
    )
    keep_w, keep_ag, removed_w, removed_ag = mass_balanced_vertex_cover(
        len(w_pos), len(ag_pos), conflicts.row, conflicts.col
    )
    if np.any(keep_w[conflicts.row] & keep_ag[conflicts.col]):
        raise RuntimeError("The conflict cover left at least one screened Ag-W pair")

    retained_rows = np.concatenate((w_rows[keep_w], ag_rows[keep_ag]))
    candidate = source[retained_rows]
    candidate_w = candidate[candidate[:, 1] == 1, 2:5]
    candidate_ag = candidate[candidate[:, 1] == 2, 2:5]
    cross_nn = nearest_distance(candidate_w, candidate_ag)
    w_nn = same_phase_nearest(candidate_w)
    ag_nn = same_phase_nearest(candidate_ag)
    report = {
        "status": "GEOMETRY_CANDIDATE_ONLY_REQUIRES_DAI_EFS_STATIC_RELAXATION",
        "method": "mass-balanced greedy vertex cover on W-Ag lattice-site conflict graph",
        "source_file": "../ag_w_qsgs.data.gz",
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "qsgs_mask": "retained unchanged from ../qsgs_morphology.npz",
        "box_A": [200.0, 200.0, 300.0],
        "material_z_A": [0.0, 200.0],
        "boundary": "p p f",
        "screening_cutoff_A": R_SCREEN,
        "screening_rationale": {
            "note": "Geometric pruning radius only; final acceptance is based on unconstrained Dai EFS static forces and structure checks.",
            "dai_efs_AgW_radial_pair_force_at_cutoff_eV_per_A": pair_force_magnitude(R_SCREEN),
            "radial_pair_force_at_2_3_A_eV_per_A": pair_force_magnitude(2.3),
            "radial_pair_force_definition": "absolute derivative of V_AgW(r); embedding-density terms are not included",
        },
        "conflict_graph_edges_removed": int(conflicts.nnz),
        "removed_atom_counts": {"W": removed_w, "Ag": removed_ag},
        "atom_counts": {"W": int(keep_w.sum()), "Ag": int(keep_ag.sum())},
        "total_atoms": int(len(candidate)),
        "atom_fractions": {
            "W": float(keep_w.sum() / len(candidate)),
            "Ag": float(keep_ag.sum() / len(candidate)),
        },
        "mass_ratio_Ag_over_W": float(keep_ag.sum() * MASS[2] / (keep_w.sum() * MASS[1])),
        "mass_ratio_relative_error": abs(float(keep_ag.sum() * MASS[2] / (keep_w.sum() * MASS[1])) - 1.0),
        "minimum_cross_distance_A": float(cross_nn.min()),
        "same_phase_minimum_neighbor_distance_A": {"W": float(w_nn.min()), "Ag": float(ag_nn.min())},
        "topology": {},
        "static_dynamics": "NOT_RUN",
        "decision": "PENDING_DAI_EFS_STATIC_RELAXATION",
    }

    sys.path.insert(0, str(PROJECT / "scripts"))
    import build_ag_w_qsgs as builder

    report["topology"] = {
        "W": builder.topology(candidate_w, "W"),
        "Ag": builder.topology(candidate_ag, "Ag"),
    }
    topology_pass = all(
        item["xy_winding_rank"] == 2
        and item["largest_component_spans_z"]
        and item["largest_component_share"] >= 0.98
        for item in report["topology"].values()
    )
    report["geometric_prechecks"] = {
        "conflict_graph_covered": True,
        "mass_ratio_error_below_0_1_percent": report["mass_ratio_relative_error"] < 0.001,
        "both_phases_bicontinuous": topology_pass,
    }
    report["decision"] = (
        "GEOMETRY_CANDIDATE_READY_FOR_STATIC_EFS_CHECKS"
        if topology_pass and report["mass_ratio_relative_error"] < 0.001
        else "REJECTED_GEOMETRY_PRECHECK"
    )

    data_path = OUT / "ag_w_qsgs_interface_v2.prepruned.data.gz"
    write_lammps(data_path, candidate)
    report["candidate_file"] = data_path.name
    report["candidate_sha256"] = hashlib.sha256(data_path.read_bytes()).hexdigest()
    (OUT / "geometry_precheck_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
