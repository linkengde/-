#!/usr/bin/env python3
"""Try a geometry-only crystal mapping with local distance constraints.

This is a hard-constraint clearance step, not a force-field minimization.
It writes an atomic candidate only if cross-phase and bulk-lattice minimum
distances can be satisfied without deleting atoms.
"""
from __future__ import annotations

import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

import build_phasefield_candidate as pf

ROOT = Path(__file__).resolve().parent
REGISTRY = (0.4, 0.6, -0.2)  # Ag(111) surface primitive u/v and d111 fractions
CROSS_MIN_A = 2.02
SAME_PHASE_COMPRESSION = 0.95
MAX_ATOM_SHIFT_A = 2.0
MAX_ITERATIONS = 220
CONVERGENCE_A = 1.0e-5


def wrap_xy(xyz, lx, ly, xy):
    xyz = np.asarray(xyz, float).copy()
    v = (xyz[:, 1] / ly) % 1.0
    u = ((xyz[:, 0] - xy * (xyz[:, 1] / ly)) / lx) % 1.0
    xyz[:, 0] = u * lx + v * xy
    xyz[:, 1] = v * ly
    return xyz


def cell_translations(lx, ly, xy):
    a = np.array([lx, 0.0, 0.0])
    b = np.array([xy, ly, 0.0])
    return np.array([i * a + j * b for i in (-1, 0, 1) for j in (-1, 0, 1)])


def phase_at(xyz, mask, lx, ly, xy, height):
    nx, ny, nz = mask.shape
    v = xyz[:, 1] / ly
    u = (xyz[:, 0] - xy * v) / lx
    i = np.floor((u % 1.0) * nx).astype(int) % nx
    j = np.floor((v % 1.0) * ny).astype(int) % ny
    valid = (xyz[:, 2] >= 0.0) & (xyz[:, 2] < height)
    k = np.clip(np.floor(xyz[:, 2] / height * nz).astype(int), 0, nz - 1)
    return valid, mask[i, j, k]


def pair_list(coords_a, coords_b, types_a, types_b, shifts, radius, same=False):
    images = np.concatenate([coords_b + shift for shift in shifts], axis=0)
    sparse = cKDTree(coords_a).sparse_distance_matrix(
        cKDTree(images), radius + 1e-8, output_type="ndarray")
    best = {}
    for row in sparse:
        i = int(row["i"])
        image = int(row["j"])
        j = image % len(coords_b)
        if same and i == j:
            continue
        key = tuple(sorted((i, j))) if same else (i, j)
        dist = float(row["v"])
        if dist < best.get(key, (0, 0, 0, math.inf))[3]:
            best[key] = (i, j, image // len(coords_b), dist)
    return list(best.values())


def reference_minima(coords, types, shifts, radius):
    pairs = pair_list(coords, coords, types, types, shifts, radius, same=True)
    minima = {}
    for i, j, _s, dist in pairs:
        key = tuple(sorted((int(types[i]), int(types[j]))))
        minima[key] = min(minima.get(key, math.inf), dist)
    return minima


def push_pair(coords_a, disp_a, coords_b, disp_b, i, j, shift_index,
              shifts, target, cap):
    vector = coords_b[j] + shifts[shift_index] - coords_a[i]
    distance = float(np.linalg.norm(vector))
    if distance >= target - 1e-10 or distance < 1e-12:
        return
    direction = vector / distance
    gap = target - distance
    available_a = max(0.0, cap - float(np.linalg.norm(disp_a[i])))
    available_b = max(0.0, cap - float(np.linalg.norm(disp_b[j])))
    move_a = min(0.5 * gap, available_a)
    move_b = min(0.5 * gap, available_b)
    remaining = gap - move_a - move_b
    if remaining > 0.0:
        extra = min(remaining, available_a - move_a)
        move_a += extra
        remaining -= extra
    if remaining > 0.0:
        extra = min(remaining, available_b - move_b)
        move_b += extra
    coords_a[i] -= move_a * direction
    disp_a[i] -= move_a * direction
    coords_b[j] += move_b * direction
    disp_b[j] += move_b * direction


def repair(coords_t, types_t, coords_ag, shifts, ref_t, ref_ag):
    disp_t = np.zeros_like(coords_t)
    disp_ag = np.zeros_like(coords_ag)
    t_targets = {key: value * SAME_PHASE_COMPRESSION for key, value in ref_t.items()}
    ag_targets = {key: value * SAME_PHASE_COMPRESSION for key, value in ref_ag.items()}
    history = []

    def max_violations():
        cross = pair_list(coords_t, coords_ag, types_t,
                          np.full(len(coords_ag), 4, int), shifts, CROSS_MIN_A)
        tsc = pair_list(coords_t, coords_t, types_t, types_t,
                        shifts, max(t_targets.values()), same=True)
        ag = pair_list(coords_ag, coords_ag, np.full(len(coords_ag), 4, int),
                       np.full(len(coords_ag), 4, int), shifts, max(ag_targets.values()), same=True)
        cross_v = max((CROSS_MIN_A - item[3] for item in cross), default=0.0)
        t_v = max((t_targets.get(tuple(sorted((int(types_t[i]), int(types_t[j])))), 0.0) - d
                   for i, j, _s, d in tsc), default=0.0)
        ag_v = max((ag_targets[(4, 4)] - d for _i, _j, _s, d in ag), default=0.0)
        return max(0.0, cross_v), max(0.0, t_v), max(0.0, ag_v)

    for iteration in range(1, MAX_ITERATIONS + 1):
        cross = pair_list(coords_t, coords_ag, types_t,
                          np.full(len(coords_ag), 4, int), shifts, CROSS_MIN_A)
        cross.sort(key=lambda row: row[3])
        for i, j, si, _d in cross:
            push_pair(coords_t, disp_t, coords_ag, disp_ag, i, j, si,
                      shifts, CROSS_MIN_A, MAX_ATOM_SHIFT_A)

        tsc_pairs = pair_list(coords_t, coords_t, types_t, types_t,
                              shifts, max(t_targets.values()), same=True)
        tsc_pairs.sort(key=lambda row: row[3])
        for i, j, si, _d in tsc_pairs:
            target = t_targets.get(tuple(sorted((int(types_t[i]), int(types_t[j])))))
            if target is not None:
                push_pair(coords_t, disp_t, coords_t, disp_t, i, j, si,
                          shifts, target, MAX_ATOM_SHIFT_A)

        ag_types = np.full(len(coords_ag), 4, int)
        ag_pairs = pair_list(coords_ag, coords_ag, ag_types, ag_types,
                             shifts, max(ag_targets.values()), same=True)
        ag_pairs.sort(key=lambda row: row[3])
        for i, j, si, _d in ag_pairs:
            push_pair(coords_ag, disp_ag, coords_ag, disp_ag, i, j, si,
                      shifts, ag_targets[(4, 4)], MAX_ATOM_SHIFT_A)

        if iteration == 1 or iteration % 10 == 0:
            violations = max_violations()
            history.append({"iteration": iteration,
                             "cross_violation_A": violations[0],
                             "TSC_lattice_violation_A": violations[1],
                             "Ag_lattice_violation_A": violations[2],
                             "max_TSC_shift_A": float(np.linalg.norm(disp_t, axis=1).max()),
                             "max_Ag_shift_A": float(np.linalg.norm(disp_ag, axis=1).max())})
            print(json.dumps(history[-1]))
            if max(violations) < CONVERGENCE_A:
                break
    return coords_t, disp_t, coords_ag, disp_ag, history, max_violations()


def all_pair_stats(coords_a, types_a, coords_b, types_b, shifts, radius, same=False):
    pairs = pair_list(coords_a, coords_b, types_a, types_b, shifts, radius, same=same)
    by_type = defaultdict(list)
    for i, j, _s, dist in pairs:
        pair_type = tuple(sorted((int(types_a[i]), int(types_b[j]))))
        by_type[pair_type].append(dist)
    return {f"{a}-{b}": {
        "pairs_within_radius": len(ds), "min_A": float(min(ds)),
        "below_1p5_A": int(sum(d < 1.5 for d in ds)),
        "below_2p0_A": int(sum(d < 2.0 for d in ds)),
        "below_2p5_A": int(sum(d < 2.5 for d in ds))}
        for (a, b), ds in sorted(by_type.items())}


def component_summary(coords, cutoff, shifts):
    n = len(coords)
    pairs = pair_list(coords, coords, np.zeros(n, int), np.zeros(n, int),
                      shifts, cutoff, same=True)
    parent = np.arange(n)

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, j, _s, _d in pairs:
        ri, rj = root(i), root(j)
        if ri != rj:
            parent[ri] = rj
    labels = np.array([root(i) for i in range(n)])
    _, counts = np.unique(labels, return_counts=True)
    largest = int(counts.max())
    label = np.unique(labels, return_counts=True)[0][int(np.argmax(counts))]
    z = coords[labels == label, 2]
    return {"cutoff_A": cutoff, "components": int(len(counts)),
            "largest_fraction": largest / n, "largest_atoms": largest,
            "largest_z_min_A": float(z.min()), "largest_z_max_A": float(z.max()),
            "largest_z_span_A": float(z.max() - z.min())}


def write_lammps_data(path, atoms, lx, ly, xy, box_z):
    masses = {1: (47.867, "Ti"), 2: (28.0855, "Si"),
              3: (12.0107, "C"), 4: (107.8682, "Ag")}
    with path.open("w") as f:
        f.write("Phase-field bicontinuous Ag-Ti3SiC2; locally clearance-screened geometry, not relaxed\n\n")
        f.write(f"{len(atoms)} atoms\n4 atom types\n\n")
        f.write(f"0.0 {lx:.12f} xlo xhi\n0.0 {ly:.12f} ylo yhi\n0.0 {box_z:.12f} zlo zhi\n")
        f.write(f"{xy:.12f} 0.0 0.0 xy xz yz\n\nMasses\n\n")
        for atom_type, (mass, name) in masses.items():
            f.write(f"{atom_type} {mass:.7f} # {name}\n")
        f.write("\nAtoms # atomic\n\n")
        for atom_id, atom in enumerate(atoms, 1):
            f.write(f"{atom_id} {atom['type']} {atom['x']:.10f} {atom['y']:.10f} {atom['z']:.10f}\n")


def main():
    mask = np.load(ROOT / "phasefield_mask_tsc.npy").astype(bool)
    bounds, tilt, source = pf.parse_data(pf.SOURCE)
    lx = bounds["x"][1] - bounds["x"][0]
    ly = bounds["y"][1] - bounds["y"][0]
    xy = tilt[0]
    shifts = cell_translations(lx, ly, xy)
    area = lx * ly
    c_tsc, d111, ag_z0, tsc_motif, ag_motifs = pf.crystal_lattices(source, lx, ly, xy)

    tsc_full, tsc_full_types = [], []
    for repeat in range(math.ceil(pf.HEIGHT / c_tsc)):
        for row in tsc_motif:
            z = float(row[4]) + repeat * c_tsc
            if z < pf.HEIGHT:
                tsc_full.append((float(row[2]), float(row[3]), z))
                tsc_full_types.append(int(row[1]))
    tsc_full = wrap_xy(np.asarray(tsc_full), lx, ly, xy)
    tsc_full_types = np.asarray(tsc_full_types, int)

    fu, fv, fw = REGISTRY
    ag_translation = (fu * np.array([lx / 13.0, 0.0, 0.0])
                      + fv * np.array([-lx / 26.0, ly / 26.0, 0.0])
                      + np.array([0.0, 0.0, fw * d111]))
    ag_full, ag_full_types = [], []
    for layer in range(math.floor(-ag_z0 / d111) - 2,
                       math.ceil((pf.HEIGHT - ag_z0) / d111) + 2):
        z = ag_z0 + layer * d111 + ag_translation[2]
        if 0.0 <= z < pf.HEIGHT:
            for row in ag_motifs[layer % 3]:
                ag_full.append((float(row[2]) + ag_translation[0],
                                float(row[3]) + ag_translation[1], z))
                ag_full_types.append(4)
    ag_full = wrap_xy(np.asarray(ag_full), lx, ly, xy)
    ag_full_types = np.asarray(ag_full_types, int)

    valid_t, in_t = phase_at(tsc_full, mask, lx, ly, xy, pf.HEIGHT)
    coords_t = tsc_full[valid_t & in_t].copy()
    types_t = tsc_full_types[valid_t & in_t].copy()
    valid_a, in_a = phase_at(ag_full, mask, lx, ly, xy, pf.HEIGHT)
    coords_ag = ag_full[valid_a & ~in_a].copy()
    types_ag = ag_full_types[valid_a & ~in_a].copy()
    coords_t0, coords_ag0 = coords_t.copy(), coords_ag.copy()
    atom_types = np.concatenate((types_t, types_ag))
    if not len(atom_types) or not np.any(atom_types == 4) or not np.any(np.isin(atom_types, (1, 2, 3))):
        raise RuntimeError("Atomic candidate must contain both Ag and Ti3SiC2 phases")

    ref_t = reference_minima(tsc_full, tsc_full_types, shifts, 5.0)
    ref_ag = reference_minima(ag_full, ag_full_types, shifts, 3.5)
    coords_t, disp_t, coords_ag, disp_ag, history, residuals = repair(
        coords_t, types_t, coords_ag, shifts, ref_t, ref_ag)
    if max(residuals) >= CONVERGENCE_A:
        result = {"status": "no candidate written; geometric constraints did not converge",
                  "residual_violations_A": {"cross": residuals[0], "TSC_lattice": residuals[1], "Ag_lattice": residuals[2]},
                  "history": history}
        (ROOT / "atomic_mapping_attempt.json").write_text(json.dumps(result, indent=2) + "\n")
        raise RuntimeError(f"Clearance constraints not converged: {result['residual_violations_A']}")

    # Apply only a common z-origin shift, preserving all relative positions.
    zmin = min(float(coords_t[:, 2].min()), float(coords_ag[:, 2].min()))
    zshift = max(0.0, -zmin)
    coords_t[:, 2] += zshift
    coords_ag[:, 2] += zshift
    coords_t_out = wrap_xy(coords_t, lx, ly, xy)
    coords_ag_out = wrap_xy(coords_ag, lx, ly, xy)
    zmax = max(float(coords_t_out[:, 2].max()), float(coords_ag_out[:, 2].max()))
    box_z = pf.HEIGHT + pf.VACUUM
    if zmax >= box_z:
        raise RuntimeError(f"Moved atom exceeds box top: zmax={zmax}")

    atom_rows = []
    for xyz, typ in zip(coords_t_out, types_t):
        atom_rows.append({"type": int(typ), "x": float(xyz[0]), "y": float(xyz[1]), "z": float(xyz[2])})
    for xyz, typ in zip(coords_ag_out, types_ag):
        atom_rows.append({"type": int(typ), "x": float(xyz[0]), "y": float(xyz[1]), "z": float(xyz[2])})
    data_path = ROOT / "AgTi3SiC2_phasefield_atomic_candidate.data"
    write_lammps_data(data_path, atom_rows, lx, ly, xy, box_z)

    counts = Counter(atom_types)
    masses = {1: 47.867, 2: 28.0855, 3: 12.0107, 4: 107.8682}
    total_mass = sum(masses[int(t)] for t in atom_types)
    ag_mass_fraction = counts[4] * masses[4] / total_mass
    rho_tsc = 1728.0 / (area * c_tsc)
    rho_ag = 169.0 / (area * d111)
    vol_tsc = sum(counts[t] for t in (1, 2, 3)) / rho_tsc
    vol_ag = counts[4] / rho_ag
    dis_t = np.linalg.norm(disp_t, axis=1)
    dis_ag = np.linalg.norm(disp_ag, axis=1)
    pair_t = all_pair_stats(coords_t_out, types_t, coords_t_out, types_t, shifts, 5.0, same=True)
    pair_ag = all_pair_stats(coords_ag_out, types_ag, coords_ag_out, types_ag, shifts, 4.0, same=True)
    pair_cross = all_pair_stats(coords_t_out, types_t, coords_ag_out, types_ag, shifts, 3.2, same=False)
    comp_t = component_summary(coords_t_out, 3.35, shifts)
    comp_ag = component_summary(coords_ag_out, 3.20, shifts)

    atom_id = 0
    with (ROOT / "atomic_displacements.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["atom_id", "type", "ideal_x_A", "ideal_y_A", "ideal_z_A",
                    "final_x_A", "final_y_A", "final_z_A", "shift_A"])
        for original, final, typ, displacement in zip(coords_t0, coords_t_out, types_t, dis_t):
            atom_id += 1
            w.writerow([atom_id, int(typ), *[f"{x:.8f}" for x in original],
                        *[f"{x:.8f}" for x in final], f"{displacement:.8f}"])
        for original, final, typ, displacement in zip(coords_ag0, coords_ag_out, types_ag, dis_ag):
            atom_id += 1
            w.writerow([atom_id, int(typ), *[f"{x:.8f}" for x in original],
                        *[f"{x:.8f}" for x in final], f"{displacement:.8f}"])

    audit = {
        "status": "geometry-only, locally clearance-screened atomic candidate; no pair coefficients, force evaluation, minimization, or dynamics",
        "method": "phase-field mask crystal-site assignment, then local simultaneous hard-distance projection; no atoms deleted",
        "atom_count_policy": "emergent from phase-mask occupancy, lattice sites, and surface termination; no target-count filter",
        "source_sha256": pf.sha256(pf.SOURCE),
        "dimensions_A": {"a_vector": [lx, 0.0], "b_vector": [xy, ly],
            "solid_reference_height": pf.HEIGHT, "z_origin_shift_A": zshift,
            "atom_z_min_A": min(float(coords_t_out[:,2].min()), float(coords_ag_out[:,2].min())),
            "atom_z_max_A": zmax, "box_z_A": box_z, "top_empty_space_A": box_z-zmax,
            "in_plane_angle_deg": 120.0, "boundary": "p p f"},
        "phase_mask_topology_6_neighbor": {
            "TSC": pf.topology(mask), "Ag": pf.topology(~mask)},
        "registry": {"Ag_surface_u_fraction": fu, "Ag_surface_v_fraction": fv,
            "Ag_layer_phase_fraction_d111": fw, "translation_A": ag_translation.tolist()},
        "composition": {"counts_by_type": {str(k): int(v) for k, v in sorted(counts.items())},
            "TSC_atoms": int(sum(counts[t] for t in (1,2,3)),), "Ag_atoms": int(counts[4]),
            "total_atoms": int(len(atom_types)), "Ag_mass_fraction": ag_mass_fraction,
            "estimated_Ag_volume_fraction_from_bulk_site_density": vol_ag/(vol_tsc+vol_ag),
            "estimated_TSC_volume_fraction_from_bulk_site_density": vol_tsc/(vol_tsc+vol_ag)},
        "clearance_projection": {"cross_phase_minimum_A": CROSS_MIN_A,
            "bulk_same_phase_minimum_fraction": SAME_PHASE_COMPRESSION,
            "maximum_per_atom_shift_A": MAX_ATOM_SHIFT_A,
            "iterations": history[-1]["iteration"] if history else 0,
            "max_residual_A": max(residuals), "history": history,
            "TSC_displacement_A": {"p50": float(np.percentile(dis_t,50)), "p90": float(np.percentile(dis_t,90)),
                "p95": float(np.percentile(dis_t,95)), "p99": float(np.percentile(dis_t,99)), "max": float(dis_t.max()),
                "count_over_0p5_A": int((dis_t>0.5).sum()), "count_over_1p0_A": int((dis_t>1.0).sum())},
            "Ag_displacement_A": {"p50": float(np.percentile(dis_ag,50)), "p90": float(np.percentile(dis_ag,90)),
                "p95": float(np.percentile(dis_ag,95)), "p99": float(np.percentile(dis_ag,99)), "max": float(dis_ag.max()),
                "count_over_0p5_A": int((dis_ag>0.5).sum()), "count_over_1p0_A": int((dis_ag>1.0).sum())}},
        "interatomic_pair_audit": {"TSC_type_pairs_within_5A": pair_t,
            "Ag_Ag_pairs_within_4A": pair_ag, "Ag_TSC_pairs_within_3p2A": pair_cross},
        "phase_atom_graph_connectivity": {"TSC_cutoff_3p35A": comp_t, "Ag_cutoff_3p20A": comp_ag},
        "caveat": "Clearance thresholds are geometric screening values, not validated equilibrium distances. This candidate is not statically or dynamically validated; inspect per-atom shifts and use only for morphology review until cross-potential validation."
    }
    (ROOT / "atomic_mapping_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps({"candidate": str(data_path), "atoms": len(atom_types),
        "atom_components_TSC": comp_t, "atom_components_Ag": comp_ag,
        "cross_min": min(v["min_A"] for v in pair_cross.values()),
        "TSC_shift_max_A": float(dis_t.max()), "Ag_shift_max_A": float(dis_ag.max()),
        "Ag_volume_fraction_est": vol_ag/(vol_tsc+vol_ag), "Ag_mass_fraction": ag_mass_fraction}, indent=2))


if __name__ == "__main__":
    main()
