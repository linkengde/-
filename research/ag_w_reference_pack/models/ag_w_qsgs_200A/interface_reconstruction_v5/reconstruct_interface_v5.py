#!/usr/bin/env python3
"""Reconstruct QSGS Ag/W interfaces using a signed-distance normal field.

The original QSGS mask is regularized without changing phase volume, then
repopulated using separate Ag-FCC and W-BCC lattice sites. The discrete site
counts can change slightly. Each phase's interface-band atoms are displaced
along a smoothed local normal, away from the other phase;
the displacement decays continuously into the same-phase lattice. This script
creates a geometry precondition only. Its output is not an MD start until the
staged spring-constrained relaxation and final spring-free EFS audit pass.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.ndimage import distance_transform_edt, gaussian_filter
from scipy.spatial import cKDTree

PACK = Path(__file__).resolve().parents[2]
MODEL = PACK / "models/ag_w_qsgs_200A"
SOURCE_DATA = MODEL / "ag_w_qsgs.data.gz"
SOURCE_MASK = MODEL / "qsgs_morphology.npz"
OUT_ROOT = Path(__file__).resolve().parent
BOX_XY_A = 200.0
MATERIAL_Z_A = 200.0
BOX_Z_A = 300.0
VOXELS_PER_SIDE = 80
VOXEL_A = 2.5
MASS_W = 183.84
MASS_AG = 107.8682
NEIGHBOR_CUTOFF_A = {1: 2.95, 2: 3.15}
DISTANCE_FLOORS_A = {"Ag_W": 2.54, "W_W": 2.60, "Ag_Ag": 2.48}


def read_lammps_atomic(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        lines = stream.readlines()
    box = []
    for axis in "xyz":
        row = next(line for line in lines if line.strip().endswith(f"{axis}lo {axis}hi"))
        lo, hi = map(float, row.split()[:2])
        box.append((lo, hi))
    section = next(i for i, line in enumerate(lines) if line.strip().startswith("Atoms")) + 1
    rows = []
    for line in lines[section:]:
        fields = line.split()
        if not fields:
            if rows:
                break
            continue
        if len(fields) < 5:
            if rows:
                break
            continue
        rows.append((int(fields[0]), int(fields[1]), *map(float, fields[2:5])))
    atoms = np.asarray(rows, dtype=float)
    return atoms[:, 0].astype(np.int64), atoms[:, 1].astype(np.int8), atoms[:, 2:5].copy()


def signed_distance_and_normal(mask_w: np.ndarray, sigma_vox: float) -> tuple[np.ndarray, np.ndarray]:
    """Signed distance (Ag positive, W negative) with xy periodicity."""
    if mask_w.shape != (VOXELS_PER_SIDE,) * 3:
        raise ValueError(f"Expected original 80^3 QSGS morphology, got {mask_w.shape}")
    tiled = np.tile(mask_w, (3, 3, 1))
    pad_z = 24
    padded = np.pad(tiled, ((0, 0), (0, 0), (pad_z, pad_z)), mode="edge")
    dist_inside_w = distance_transform_edt(padded, sampling=(VOXEL_A,) * 3)
    dist_inside_ag = distance_transform_edt(~padded, sampling=(VOXEL_A,) * 3)
    phi_padded = (dist_inside_ag - dist_inside_w).astype(np.float32)
    zslice = slice(pad_z, pad_z + VOXELS_PER_SIDE)
    phi = phi_padded[VOXELS_PER_SIDE:2 * VOXELS_PER_SIDE,
                     VOXELS_PER_SIDE:2 * VOXELS_PER_SIDE, zslice].copy()
    del dist_inside_w, dist_inside_ag, padded, phi_padded
    normal_field = gaussian_filter(phi, sigma=sigma_vox, mode=("wrap", "wrap", "nearest"))
    gx = (np.roll(normal_field, -1, axis=0) - np.roll(normal_field, 1, axis=0)) / (2.0 * VOXEL_A)
    gy = (np.roll(normal_field, -1, axis=1) - np.roll(normal_field, 1, axis=1)) / (2.0 * VOXEL_A)
    gz = np.gradient(normal_field, VOXEL_A, axis=2, edge_order=2)
    normal = np.stack((gx, gy, gz), axis=-1)
    norms = np.linalg.norm(normal, axis=-1)
    good = norms > 1.0e-8
    normal[good] /= norms[good, None]
    normal[~good] = 0.0
    return phi, normal.astype(np.float32)


def interpolate_grid(field: np.ndarray, xyz: np.ndarray) -> np.ndarray:
    """Trilinear interpolation at voxel centers; xy wraps, z clamps."""
    u = xyz / VOXEL_A - 0.5
    lo = np.floor(u).astype(np.int64)
    frac = u - lo
    nx, ny, nz = field.shape[:3]
    out_shape = (len(xyz),) + field.shape[3:]
    out = np.zeros(out_shape, dtype=np.float64)
    for bx in (0, 1):
        ix = (lo[:, 0] + bx) % nx
        wx = frac[:, 0] if bx else 1.0 - frac[:, 0]
        for by in (0, 1):
            iy = (lo[:, 1] + by) % ny
            wy = frac[:, 1] if by else 1.0 - frac[:, 1]
            for bz in (0, 1):
                iz = np.clip(lo[:, 2] + bz, 0, nz - 1)
                wz = frac[:, 2] if bz else 1.0 - frac[:, 2]
                weight = (wx * wy * wz).reshape((-1,) + (1,) * (field.ndim - 3))
                out += weight * field[ix, iy, iz]
    return out


def pair_metrics(pos: np.ndarray, types: np.ndarray) -> dict:
    w_idx = np.flatnonzero(types == 1)
    ag_idx = np.flatnonzero(types == 2)
    wtree = cKDTree(pos[w_idx], boxsize=[BOX_XY_A, BOX_XY_A, 1000.0])
    agtree = cKDTree(pos[ag_idx], boxsize=[BOX_XY_A, BOX_XY_A, 1000.0])
    min_cross = float(min(wtree.query(pos[ag_idx], workers=-1)[0].min(),
                          agtree.query(pos[w_idx], workers=-1)[0].min()))
    min_ww = float(wtree.query(pos[w_idx], k=2, workers=-1)[0][:, 1].min())
    min_ag = float(agtree.query(pos[ag_idx], k=2, workers=-1)[0][:, 1].min())
    n_cross = int(wtree.sparse_distance_matrix(agtree, DISTANCE_FLOORS_A["Ag_W"], output_type="coo_matrix").nnz)
    n_ww = int(len(wtree.query_pairs(DISTANCE_FLOORS_A["W_W"], output_type="ndarray")))
    n_ag = int(len(agtree.query_pairs(DISTANCE_FLOORS_A["Ag_Ag"], output_type="ndarray")))
    return {
        "minimum_distances_A": {"Ag_W": min_cross, "W_W": min_ww, "Ag_Ag": min_ag},
        "pairs_below_frozen_screen": {"Ag_W": n_cross, "W_W": n_ww, "Ag_Ag": n_ag},
    }


def write_lammps_atomic(path: Path, atom_ids: np.ndarray, types: np.ndarray, pos: np.ndarray) -> None:
    with gzip.open(path, "wt", encoding="utf-8", compresslevel=6) as stream:
        stream.write("Ag-W QSGS interface reconstruction v5 geometry precondition\n\n")
        stream.write(f"{len(pos)} atoms\n\n2 atom types\n\n")
        stream.write("0 200 xlo xhi\n0 200 ylo yhi\n0 300 zlo zhi\n\n")
        stream.write("Masses\n\n1 183.84 # W\n2 107.8682 # Ag\n\nAtoms # atomic\n\n")
        for atom_id, atom_type, xyz in zip(atom_ids, types, pos):
            stream.write(f"{int(atom_id)} {int(atom_type)} {xyz[0]:.12f} {xyz[1]:.12f} {xyz[2]:.12f}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUT_ROOT / "geometry_preconditioned")
    parser.add_argument("--half-gap-A", type=float, default=1.35,
                        help="signed normal shift amplitude per phase at the interface")
    parser.add_argument("--decay-A", type=float, default=10.0,
                        help="Gaussian decay length for the lattice-preserving displacement")
    parser.add_argument("--mask-smoothing-vox", type=float, default=2.0,
                        help="Gaussian regularization before thresholding to the original QSGS W volume")
    parser.add_argument("--normal-smoothing-vox", type=float, default=1.0)
    args = parser.parse_args()
    if args.half_gap_A <= 0 or args.decay_A <= 0 or args.normal_smoothing_vox < 0 or args.mask_smoothing_vox < 0:
        parser.error("normal shift/decay must be positive and smoothing non-negative")
    out_dir = args.output.resolve()
    if out_dir.exists() and any(out_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output directory {out_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)

    source_ids, source_types, source_pos = read_lammps_atomic(SOURCE_DATA)
    morphology = np.load(SOURCE_MASK)
    source_mask = morphology["W_backbone"].astype(bool)
    if len(source_ids) != 482_471 or set(np.unique(source_types)) != {1, 2}:
        raise RuntimeError("Unexpected original QSGS atom count or type convention")
    if abs(float((source_types == 2).sum() * MASS_AG / ((source_types == 1).sum() * MASS_W)) - 1.0) > 1e-3:
        raise RuntimeError("Original Ag/W source mass ratio unexpectedly fails the 0.1% gate")

    sys.path.insert(0, str(PACK / "scripts"))
    import build_ag_w_qsgs as builder
    if args.mask_smoothing_vox > 0.0:
        occupancy = gaussian_filter(source_mask.astype(np.float32), sigma=args.mask_smoothing_vox,
                                    mode=("wrap", "wrap", "nearest"))
        n_w_voxels = int(source_mask.sum())
        take = np.argsort(occupancy.ravel(), kind="stable")[-n_w_voxels:]
        w_mask = np.zeros(occupancy.size, dtype=bool)
        w_mask[take] = True
        w_mask = w_mask.reshape(source_mask.shape)
    else:
        w_mask = source_mask.copy()
    w_pool, w_strain = builder.make_lattice_sites("W")
    ag_pool, ag_strain = builder.make_lattice_sites("Ag")
    w_positions = w_pool[builder.in_mask(w_mask, w_pool)]
    ag_positions = ag_pool[~builder.in_mask(w_mask, ag_pool)]
    types = np.concatenate((np.ones(len(w_positions), dtype=np.int8),
                            np.full(len(ag_positions), 2, dtype=np.int8)))
    original = np.vstack((w_positions, ag_positions))
    atom_ids = np.arange(1, len(types) + 1, dtype=np.int64)

    phi, normal_field = signed_distance_and_normal(w_mask, args.normal_smoothing_vox)
    phi_atom = interpolate_grid(phi, original)
    normal_atom = interpolate_grid(normal_field, original)
    norm = np.linalg.norm(normal_atom, axis=1)
    active = (np.abs(phi_atom) < 4.0 * args.decay_A) & (norm > 0.25)
    normal_atom[active] /= norm[active, None]
    normal_atom[~active] = 0.0
    amplitude = args.half_gap_A * np.exp(-0.5 * (np.abs(phi_atom) / args.decay_A) ** 2)
    phase_sign = np.where(types == 1, -1.0, 1.0)
    displacement = normal_atom * (phase_sign * amplitude)[:, None]
    proposed = original + displacement
    # Enforce the requested free-surface extent without moving atoms along a
    # direction unrelated to the local interface normal.
    scale = np.ones(len(original))
    positive = (displacement[:, 2] > 0.0) & (proposed[:, 2] > MATERIAL_Z_A)
    negative = (displacement[:, 2] < 0.0) & (proposed[:, 2] < 0.0)
    scale[positive] = np.maximum(0.0, (MATERIAL_Z_A - original[positive, 2]) / displacement[positive, 2])
    scale[negative] = np.maximum(0.0, -original[negative, 2] / displacement[negative, 2])
    displacement *= scale[:, None]
    pos = original + displacement
    pos[:, :2] %= BOX_XY_A

    # Use the full, strain-adjusted lattice populations to define phase
    # number-density targets for the regularized QSGS mask.
    w_volume = float(w_mask.mean()) * BOX_XY_A**2 * MATERIAL_Z_A
    ag_volume = float((~w_mask).mean()) * BOX_XY_A**2 * MATERIAL_Z_A
    n_w, n_ag = int(np.count_nonzero(types == 1)), int(np.count_nonzero(types == 2))
    density_w = n_w / w_volume
    density_ag = n_ag / ag_volume
    target_density_w = len(w_pool) / (BOX_XY_A**2 * MATERIAL_Z_A)
    target_density_ag = len(ag_pool) / (BOX_XY_A**2 * MATERIAL_Z_A)
    density_error = {
        "W_relative": density_w / target_density_w - 1.0,
        "Ag_relative": density_ag / target_density_ag - 1.0,
    }

    input_distances = pair_metrics(original, types)
    candidate_distances = pair_metrics(pos, types)
    topology = {
        "W": builder.topology(pos[types == 1], "W"),
        "Ag": builder.topology(pos[types == 2], "Ag"),
    }
    topology_pass = all(v["largest_component_share"] >= 0.98 and v["xy_winding_rank"] == 2
                        and v["largest_component_spans_z"] for v in topology.values())
    box_pass = bool(np.all(pos[:, 2] >= -1.0e-8) and np.all(pos[:, 2] <= MATERIAL_Z_A + 1.0e-8))
    distance_pass = all(candidate_distances["pairs_below_frozen_screen"][k] == 0
                        for k in ("Ag_W", "W_W", "Ag_Ag"))
    density_pass = max(abs(v) for v in density_error.values()) <= 0.01
    mass_ratio = n_ag * MASS_AG / (n_w * MASS_W)
    mass_pass = abs(mass_ratio - 1.0) <= 1.0e-3
    geometry_preflight = mass_pass and density_pass and topology_pass and box_pass and distance_pass

    candidate_path = out_dir / "ag_w_qsgs_v5_interface_precondition.data.gz"
    write_lammps_atomic(candidate_path, atom_ids, types, pos)
    digest = hashlib.sha256(candidate_path.read_bytes()).hexdigest()
    report = {
        "title": "Ag-W QSGS interface reconstruction v5 signed-normal geometry precondition",
        "status": "GEOMETRY_PREFLIGHT_PASS_STATIC_EFS_NOT_RUN" if geometry_preflight else "GEOMETRY_PREFLIGHT_FAIL_STATIC_EFS_NOT_RUN",
        "md_start_gate": "FAIL",
        "static_efs_status": "NOT_RUN",
        "method": "Regularize the 2.5 A QSGS interface at fixed W voxel volume, repopulate both phases with their own strained ideal lattices, then shift interface-band atoms along smoothed signed-distance normals away from the unlike phase. No atom deletion is used; the regenerated lattice population may differ slightly from the source count.",
        "source": {
            "data_file": str(SOURCE_DATA.relative_to(PACK)),
            "data_sha256": hashlib.sha256(SOURCE_DATA.read_bytes()).hexdigest(),
            "morphology_file": str(SOURCE_MASK.relative_to(PACK)),
            "qsgs_voxels": list(w_mask.shape),
            "voxel_size_A": VOXEL_A,
            "original_W_voxel_fraction": float(source_mask.mean()),
            "regularized_W_voxel_fraction": float(w_mask.mean()),
            "morphology_voxels_reassigned": int(np.count_nonzero(w_mask != source_mask)),
            "morphology_voxel_fraction_changed": float(np.mean(w_mask != source_mask)),
            "source_atom_counts": {"W": int(np.count_nonzero(source_types == 1)), "Ag": int(np.count_nonzero(source_types == 2)), "total": int(len(source_types))},
        },
        "parameters": {"mask_smoothing_vox": args.mask_smoothing_vox,
                       "half_gap_A": args.half_gap_A, "decay_A": args.decay_A,
                       "normal_smoothing_vox": args.normal_smoothing_vox,
                       "frozen_screen_distances_A": DISTANCE_FLOORS_A},
        "box_A": [BOX_XY_A, BOX_XY_A, BOX_Z_A],
        "material_A": [BOX_XY_A, BOX_XY_A, MATERIAL_Z_A],
        "boundary": "p p f",
        "counts": {"W": n_w, "Ag": n_ag, "total": int(len(types)),
                   "atom_count_change_from_source": int(len(types) - len(source_types)),
                   "atom_count_change_fraction": float(len(types) / len(source_types) - 1.0)},
        "mass_ratio_Ag_over_W": mass_ratio,
        "mass_ratio_error_fraction": abs(mass_ratio - 1.0),
        "phase_number_density_atoms_A3": {"W": density_w, "Ag": density_ag},
        "target_phase_number_density_atoms_A3": {"W": target_density_w, "Ag": target_density_ag},
        "phase_density_error_fraction": density_error,
        "density_gate_pass": density_pass,
        "minimum_distances_in_original_source_A": pair_metrics(source_pos, source_types)["minimum_distances_A"],
        "minimum_distances_before_normal_precondition_A": input_distances["minimum_distances_A"],
        "minimum_distances_after_geometry_precondition_A": candidate_distances["minimum_distances_A"],
        "pairs_below_frozen_screen_before": input_distances["pairs_below_frozen_screen"],
        "pairs_below_frozen_screen_after": candidate_distances["pairs_below_frozen_screen"],
        "distance_screen_pass": distance_pass,
        "maximum_displacement_A": float(np.linalg.norm(displacement, axis=1).max()),
        "displacement_percentiles_A": {str(q): float(np.percentile(np.linalg.norm(displacement, axis=1), q)) for q in (50, 90, 95, 99, 99.9)},
        "atoms_moved": int(np.count_nonzero(np.linalg.norm(displacement, axis=1) > 1.0e-8)),
        "atoms_at_material_z_limits_after_surface_guard": int(np.count_nonzero((pos[:, 2] < 1.0e-8) | (pos[:, 2] > MATERIAL_Z_A - 1.0e-8))),
        "topology": topology,
        "topology_gate_pass": topology_pass,
        "material_bounds_gate_pass": box_pass,
        "mass_gate_pass": mass_pass,
        "candidate": {"file": candidate_path.name, "sha256": digest},
        "warning": "Geometry-only preflight is never an MD-start pass. Bonded same-phase elastic preconditioning and staged local EFS relaxation must be followed by fully spring-free EFS convergence and a separate static audit.",
        "lattice_reference": {"W_xy_strain": w_strain, "Ag_xy_strain": ag_strain},
    }
    (out_dir / "geometry_precondition_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: report[k] for k in (
        "status", "md_start_gate", "counts", "mass_ratio_Ag_over_W", "phase_density_error_fraction",
        "minimum_distances_in_original_source_A", "minimum_distances_before_normal_precondition_A",
        "minimum_distances_after_geometry_precondition_A",
        "pairs_below_frozen_screen_after", "maximum_displacement_A", "topology_gate_pass",
    )}, indent=2))


if __name__ == "__main__":
    main()
