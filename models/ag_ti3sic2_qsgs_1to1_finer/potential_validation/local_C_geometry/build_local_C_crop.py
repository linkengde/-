#!/usr/bin/env python3
"""Extract a source-centered local geometry from the existing QSGS atom model.

This creates geometry only. It intentionally does not relax atoms or select a
cross-phase potential; those steps are blocked until a defensible interaction
model is identified and passes the local cold-state gate.
"""

from __future__ import annotations

import csv
import hashlib
import gzip
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import connected_components

MODEL = Path(__file__).resolve().parents[2]
SOURCE = MODEL / "potential_validation/relaxed_backbone_candidate/ag_ti3sic2_qsgs.data.gz"
GEOMETRY = MODEL / "potential_validation/relaxed_backbone_candidate/geometry_report.json"
OUT = Path(__file__).resolve().parent
HALF_WIDTH_A = 80.0
SIDE_BUFFER_A = 20.0
BOTTOM_FIXED_A = 20.0
HOT_RADIUS_A = 30.0
SKIN_DEPTH_A = 20.0
HOT_DEPTH_A = 60.0
MASSES = {1: 47.867, 2: 28.085, 3: 12.011, 4: 107.8682}
SYMBOLS = {1: "Ti", 2: "Si", 3: "C", 4: "Ag"}


def parse_lammps_data(path: Path):
    with gzip.open(path, "rt") as f:
        lines = f.readlines()
    n_atoms = None
    bounds = {}
    for line in lines:
        fields = line.split()
        if len(fields) >= 2 and fields[1] == "atoms":
            n_atoms = int(fields[0])
        if len(fields) >= 4 and fields[-2:] == ["xlo", "xhi"]:
            bounds[0] = (float(fields[0]), float(fields[1]))
        if len(fields) >= 4 and fields[-2:] == ["ylo", "yhi"]:
            bounds[1] = (float(fields[0]), float(fields[1]))
        if len(fields) >= 4 and fields[-2:] == ["zlo", "zhi"]:
            bounds[2] = (float(fields[0]), float(fields[1]))
    atom_line = next(i for i, line in enumerate(lines) if line.startswith("Atoms"))
    rows = []
    for line in lines[atom_line + 1 :]:
        if not line.strip():
            if rows:
                break
            continue
        fields = line.split()
        if len(fields) < 5:
            break
        rows.append((int(fields[0]), int(fields[1]), *map(float, fields[2:5])))
    data = np.asarray(rows, dtype=float)
    if n_atoms != len(data) or len(bounds) != 3:
        raise ValueError(f"LAMMPS data parse mismatch: expected {n_atoms}, read {len(data)}")
    return data, bounds


def component_summary(xyz: np.ndarray, cutoff: float):
    if len(xyz) == 0:
        return {"components": 0, "largest_atoms": 0, "atoms_outside_largest": 0}
    pairs = cKDTree(xyz).query_pairs(cutoff, output_type="ndarray")
    if len(pairs) == 0:
        return {"components": len(xyz), "largest_atoms": 1, "atoms_outside_largest": len(xyz) - 1}
    vals = np.ones(len(pairs), dtype=np.uint8)
    graph = sparse.coo_matrix(
        (np.r_[vals, vals], (np.r_[pairs[:, 0], pairs[:, 1]], np.r_[pairs[:, 1], pairs[:, 0]])),
        shape=(len(xyz), len(xyz)),
    ).tocsr()
    ncomp, labels = connected_components(graph, directed=False)
    counts = np.bincount(labels)
    return {
        "components": int(ncomp),
        "largest_atoms": int(counts.max()),
        "atoms_outside_largest": int(len(xyz) - counts.max()),
        "cutoff_A": float(cutoff),
    }


def distance_stats(values: np.ndarray):
    return {
        "minimum_A": float(np.min(values)),
        "p01_A": float(np.percentile(values, 1)),
        "p05_A": float(np.percentile(values, 5)),
        "median_A": float(np.median(values)),
        "p95_A": float(np.percentile(values, 95)),
        "count_below_2p0A": int(np.count_nonzero(values < 2.0)),
        "count_below_2p2A": int(np.count_nonzero(values < 2.2)),
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    geom = json.loads(GEOMETRY.read_text())
    lx, ly, material_lz = map(float, geom["material_dimensions_A"])
    box_lz = float(geom["simulation_box_A"][2])
    raw, bounds = parse_lammps_data(SOURCE)
    ids = raw[:, 0].astype(np.int64)
    types = raw[:, 1].astype(np.int8)
    xyz = raw[:, 2:5]

    cx, cy = lx / 2.0, ly / 2.0
    dx = (xyz[:, 0] - cx + lx / 2.0) % lx - lx / 2.0
    dy = (xyz[:, 1] - cy + ly / 2.0) % ly - ly / 2.0
    keep = (np.abs(dx) <= HALF_WIDTH_A) & (np.abs(dy) <= HALF_WIDTH_A) & (xyz[:, 2] <= material_lz + 1e-8)
    local_ids = np.flatnonzero(keep)
    local_types = types[keep]
    local_xyz = np.column_stack((dx[keep] + HALF_WIDTH_A, dy[keep] + HALF_WIDTH_A, xyz[keep, 2]))
    local_orig_ids = ids[keep]

    # Make a local atomic-style LAMMPS data file with a fixed-width side buffer.
    out_data = OUT / "C_local_crop.data.gz"
    with gzip.open(out_data, "wt") as f:
        f.write("Source-centered QSGS local crop; geometry only, not cold-state validated\n\n")
        f.write(f"{len(local_ids)} atoms\n4 atom types\n\n")
        f.write(f"0.0 {2*HALF_WIDTH_A:.8f} xlo xhi\n")
        f.write(f"0.0 {2*HALF_WIDTH_A:.8f} ylo yhi\n")
        f.write(f"0.0 {box_lz:.8f} zlo zhi\n\nMasses\n\n")
        for atom_type in range(1, 5):
            f.write(f"{atom_type} {MASSES[atom_type]:.8f} # {SYMBOLS[atom_type]}\n")
        f.write("\nAtoms # atomic\n\n")
        for new_id, (atom_type, pos) in enumerate(zip(local_types, local_xyz), 1):
            f.write(f"{new_id} {int(atom_type)} {pos[0]:.10f} {pos[1]:.10f} {pos[2]:.10f}\n")

    id_map = OUT / "C_local_crop_id_map.csv.gz"
    with gzip.open(id_map, "wt", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(("local_id", "source_id", "type", "element", "x_A", "y_A", "z_A"))
        for new_id, (source_id, atom_type, pos) in enumerate(zip(local_orig_ids, local_types, local_xyz), 1):
            writer.writerow((new_id, int(source_id), int(atom_type), SYMBOLS[int(atom_type)], *[f"{v:.10f}" for v in pos]))

    counts = {SYMBOLS[t]: int(np.count_nonzero(local_types == t)) for t in range(1, 5)}
    ag = local_xyz[local_types == 4]
    ceramic = local_xyz[np.isin(local_types, (1, 2, 3))]
    c = local_xyz[local_types == 3]
    ceramic_mass = sum(counts[s] * MASSES[t] for t, s in ((1, "Ti"), (2, "Si"), (3, "C")))
    ag_mass = counts["Ag"] * MASSES[4]
    nearest_ceramic, _ = cKDTree(ceramic).query(ag, k=1)
    nearest_carbon, _ = cKDTree(c).query(ag, k=1)
    nearest_by_species = {}
    for atom_type in (1, 2, 3):
        target = local_xyz[local_types == atom_type]
        nearest, _ = cKDTree(target).query(ag, k=1)
        nearest_by_species[SYMBOLS[atom_type]] = distance_stats(nearest)
    relx, rely = local_xyz[:, 0] - HALF_WIDTH_A, local_xyz[:, 1] - HALF_WIDTH_A
    r = np.sqrt(relx**2 + rely**2)
    depth = material_lz - local_xyz[:, 2]
    agmask = local_types == 4
    hot_core = agmask & (r <= HOT_RADIUS_A) & (depth >= -1e-8) & (depth <= HOT_DEPTH_A)
    hot_core_xyz = local_xyz[hot_core]
    hot_core_nearest_by_species = {}
    for atom_type in (1, 2, 3):
        target = local_xyz[local_types == atom_type]
        nearest, _ = cKDTree(target).query(hot_core_xyz, k=1)
        hot_core_nearest_by_species[SYMBOLS[atom_type]] = distance_stats(nearest)
    core_extent = agmask & (r <= 60.0) & (depth >= -1e-8) & (depth <= 100.0)
    side_fixed = (local_xyz[:, 0] < SIDE_BUFFER_A) | (local_xyz[:, 0] > 2*HALF_WIDTH_A-SIDE_BUFFER_A) | (local_xyz[:, 1] < SIDE_BUFFER_A) | (local_xyz[:, 1] > 2*HALF_WIDTH_A-SIDE_BUFFER_A)
    bottom_fixed = local_xyz[:, 2] <= BOTTOM_FIXED_A

    report = {
        "status": "C_LOCAL_GEOMETRY_CROP_NOT_COLD_VALIDATED",
        "source_data": str(SOURCE.relative_to(MODEL)),
        "source_data_sha256": sha256(SOURCE),
        "source_geometry_report_sha256": sha256(GEOMETRY),
        "source_geometry_status": geom.get("status"),
        "source_model_atoms": int(len(raw)),
        "qsgs_seed": geom.get("qsgs", {}).get("seed"),
        "crop_center_source_xy_A": [cx, cy],
        "crop_xy_dimensions_A": [2*HALF_WIDTH_A, 2*HALF_WIDTH_A],
        "crop_material_z_A": material_lz,
        "crop_simulation_box_z_A": box_lz,
        "lateral_boundary": "nonperiodic cut with 20 A fixed side buffer planned",
        "bottom_boundary": f"bottom {BOTTOM_FIXED_A:g} A planned fixed",
        "free_surface": "retained at original material z top with source vacuum",
        "atom_count": int(len(local_ids)),
        "element_counts": counts,
        "local_Ag_to_ceramic_mass_ratio": float(ag_mass / ceramic_mass),
        "hot_source_core_Ag_atoms_r_le_30_depth_le_60": int(hot_core.sum()),
        "source_buffer_Ag_atoms_r_le_60_depth_le_100": int(core_extent.sum()),
        "planned_fixed_side_buffer_atoms": int(side_fixed.sum()),
        "planned_fixed_bottom_atoms": int(bottom_fixed.sum()),
        "phase_geometry_connectivity": {
            "Ag_cutoff_3p98_A": component_summary(ag, 3.98),
            "Ti3SiC2_cutoff_3p25_A": component_summary(ceramic, 3.25),
        },
        "nearest_Ag_to_any_skeleton_distance": distance_stats(nearest_ceramic),
        "nearest_Ag_to_skeleton_species": nearest_by_species,
        "nearest_Ag_to_C_distance": distance_stats(nearest_carbon),
        "hot_core_Ag_nearest_to_skeleton_species": hot_core_nearest_by_species,
        "interaction_model_status": "BLOCKED: no repository Ag-Ti/Si/C model passes current force/separation-force development gates for MD/TTM use",
        "limitations": [
            "The crop preserves the source QSGS coordinates; it is not a reconstructed or relaxed local interface.",
            "The side and bottom fixed layers are planned boundary groups, not yet dynamically validated.",
            "No potential evaluation, minimization, thermostat, TTM, heating, or phase claim is included.",
        ],
        "artifacts": {
            "data": out_data.name,
            "id_map": id_map.name,
        },
    }
    (OUT / "C_local_geometry_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
