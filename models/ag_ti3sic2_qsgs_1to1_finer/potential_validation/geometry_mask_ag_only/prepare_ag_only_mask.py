#!/usr/bin/env python3
"""Export the existing QSGS Ag network and build a signed geometric mask.

No Ti/Si/C atoms or Ag--skeleton cross interactions are exported. The mask is
an external geometric constraint, not an atomistic interface potential.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import distance_transform_edt

MODEL = Path(__file__).resolve().parents[2]
SOURCE_DATA = MODEL / "potential_validation/relaxed_backbone_candidate/ag_ti3sic2_qsgs.data.gz"
SOURCE_GEOMETRY = MODEL / "potential_validation/relaxed_backbone_candidate/geometry_report.json"
SOURCE_MORPHOLOGY = MODEL / "potential_validation/relaxed_backbone_candidate/qsgs_morphology.npz"
OUT = Path(__file__).resolve().parent
AG_MASS = 107.8682
WALL_CLEARANCE_A = 0.0
WALL_STIFFNESS_EV_A2 = 5.0
AG_EAM_A0_A = 4.09
FCC_CELLS_XYZ = np.array([73, 74, 39], dtype=int)


def sha256(path: Path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_ag(path: Path):
    with gzip.open(path, "rt") as f:
        lines = f.readlines()
    natoms = next(int(line.split()[0]) for line in lines if len(line.split()) >= 2 and line.split()[1] == "atoms")
    bounds = {}
    for line in lines:
        fields = line.split()
        if len(fields) >= 4 and fields[-2:] == ["xlo", "xhi"]:
            bounds[0] = tuple(map(float, fields[:2]))
        elif len(fields) >= 4 and fields[-2:] == ["ylo", "yhi"]:
            bounds[1] = tuple(map(float, fields[:2]))
        elif len(fields) >= 4 and fields[-2:] == ["zlo", "zhi"]:
            bounds[2] = tuple(map(float, fields[:2]))
    section = next(i for i, line in enumerate(lines) if line.startswith("Atoms"))
    ids, types, xyz = [], [], []
    for line in lines[section + 1 :]:
        if not line.strip():
            if ids:
                break
            continue
        f = line.split()
        if len(f) < 5:
            break
        if int(f[1]) == 4:
            ids.append(int(f[0]))
            types.append(int(f[1]))
            xyz.append([float(v) for v in f[2:5]])
    if len(bounds) != 3 or natoms != 916534:
        raise ValueError(f"Unexpected parent data header: atoms={natoms}, bounds={bounds}")
    return np.asarray(ids, np.int64), np.asarray(xyz, float), bounds


def trilinear_periodic_xy(phi: np.ndarray, xyz: np.ndarray, spacing: np.ndarray, lz: float):
    nx, ny, nz = phi.shape
    ijk = xyz / spacing[None, :] - 0.5
    # Clamp the free-surface coordinate itself before interpolation. Clamping
    # only the integer index while retaining its old fractional part blends
    # the first/last z voxels and overcounts wall penetration at z boundaries.
    ijk[:, 2] = np.clip(ijk[:, 2], 0.0, nz - 1.0)
    base = np.floor(ijk).astype(np.int64)
    frac = ijk - base
    base[:, 0] %= nx
    base[:, 1] %= ny
    outside_z = (xyz[:, 2] < 0.0) | (xyz[:, 2] >= lz)
    base[:, 2] = np.clip(base[:, 2], 0, nz - 1)
    out = np.zeros(len(xyz), dtype=np.float64)
    for ox in (0, 1):
        ix = (base[:, 0] + ox) % nx
        wx = frac[:, 0] if ox else 1.0 - frac[:, 0]
        for oy in (0, 1):
            iy = (base[:, 1] + oy) % ny
            wy = frac[:, 1] if oy else 1.0 - frac[:, 1]
            for oz in (0, 1):
                iz = np.clip(base[:, 2] + oz, 0, nz - 1)
                wz = frac[:, 2] if oz else 1.0 - frac[:, 2]
                out += wx * wy * wz * phi[ix, iy, iz]
    out[outside_z] = 1.0e6
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    geom = json.loads(SOURCE_GEOMETRY.read_text())
    maskdoc = np.load(SOURCE_MORPHOLOGY)
    skeleton = maskdoc["skeleton"].astype(bool)
    size = maskdoc["size_A"].astype(float)
    if not np.allclose(size, geom["material_dimensions_A"], rtol=0.0, atol=1e-8):
        raise ValueError("Mask size does not match the geometry/data parent")
    ids, parent_xyz, bounds = parse_ag(SOURCE_DATA)
    parent_lattice_a = np.array([size[0] / FCC_CELLS_XYZ[0], size[1] / FCC_CELLS_XYZ[1], AG_EAM_A0_A])
    half_parent_a = parent_lattice_a / 2.0
    lattice_indices = np.rint(parent_xyz / half_parent_a[None, :]).astype(np.int64)
    lattice_error = parent_xyz - lattice_indices * half_parent_a[None, :]
    if np.max(np.abs(lattice_error)) > 1.0e-7:
        raise ValueError(f"Ag atoms do not map to the expected parent FCC sites: max error={np.max(np.abs(lattice_error))}")
    parity_violations = int(np.count_nonzero(np.mod(lattice_indices.sum(axis=1), 2)))
    if parity_violations:
        raise ValueError(f"Ag sites violate the FCC basis parity: {parity_violations}")

    phase_path = OUT / "fcc_lattice_phase_search.json"
    phase_result = json.loads(phase_path.read_text()) if phase_path.exists() else {
        "selected": {"orientation": "identity", "dx_A": 0.0, "dy_A": 0.0}
    }
    phase = phase_result["selected"]
    orientation = phase["orientation"]
    half_a = AG_EAM_A0_A / 2.0
    cells = FCC_CELLS_XYZ.copy()
    if orientation == "identity":
        output_indices = lattice_indices.copy()
    elif orientation == "rotate_z_90":
        output_indices = np.column_stack((lattice_indices[:, 1], 2 * FCC_CELLS_XYZ[0] - lattice_indices[:, 0], lattice_indices[:, 2]))
        cells = FCC_CELLS_XYZ[[1, 0, 2]]
    elif orientation == "rotate_z_180":
        output_indices = np.column_stack((2 * FCC_CELLS_XYZ[0] - lattice_indices[:, 0], 2 * FCC_CELLS_XYZ[1] - lattice_indices[:, 1], lattice_indices[:, 2]))
    elif orientation == "rotate_z_270":
        output_indices = np.column_stack((2 * FCC_CELLS_XYZ[1] - lattice_indices[:, 1], lattice_indices[:, 0], lattice_indices[:, 2]))
        cells = FCC_CELLS_XYZ[[1, 0, 2]]
    else:
        raise ValueError(f"Unknown FCC phase-search orientation: {orientation}")
    transformed_parity_violations = int(np.count_nonzero(np.mod(output_indices.sum(axis=1), 2)))
    if transformed_parity_violations:
        raise ValueError(f"FCC transformation changed lattice parity for {transformed_parity_violations} sites")

    size = cells.astype(float) * AG_EAM_A0_A
    xyz = output_indices * half_a
    xyz[:, 0] = np.mod(xyz[:, 0] + float(phase["dx_A"]), size[0])
    xyz[:, 1] = np.mod(xyz[:, 1] + float(phase["dy_A"]), size[1])
    vacuum_z = bounds[2][1] - geom["material_dimensions_A"][2]
    bounds = {0: (0.0, float(size[0])), 1: (0.0, float(size[1])), 2: (0.0, float(size[2] + vacuum_z))}
    lx, ly, lz = size
    grid = np.asarray(skeleton.shape, dtype=int)
    spacing = size / grid

    # Tile the periodic xy mask, while padding z with accessible space so the
    # geometric wall does not close the free top/bottom surfaces.
    tiled = np.tile(skeleton, (3, 3, 1))
    pad_z = 12
    padded = np.pad(tiled, ((0, 0), (0, 0), (pad_z, pad_z)), constant_values=False)
    outside_dist = distance_transform_edt(~padded, sampling=spacing)
    inside_dist = distance_transform_edt(padded, sampling=spacing)
    half_voxel = 0.5 * float(np.min(spacing))
    signed = np.where(padded, -(inside_dist - half_voxel), outside_dist - half_voxel)
    signed = signed[grid[0] : 2 * grid[0], grid[1] : 2 * grid[1], pad_z : pad_z + grid[2]].astype(np.float32)
    wall_sdf = OUT / "qsgs_signed_distance_A.npz"
    np.savez_compressed(wall_sdf, phi_A=signed, spacing_A=spacing, size_A=size, clearance_A=WALL_CLEARANCE_A, source_seed=int(maskdoc["seed"]))

    phi_atoms = trilinear_periodic_xy(signed, xyz, spacing, lz)
    active = phi_atoms < WALL_CLEARANCE_A
    inside = phi_atoms <= 0.0

    ag_data = OUT / "Ag_only_qsgs.data.gz"
    with gzip.open(ag_data, "wt") as f:
        f.write("QSGS Ag continuous phase only; QSGS geometry mask is external\n\n")
        f.write(f"{len(ids)} atoms\n1 atom types\n\n")
        for axis, name in enumerate(("x", "y", "z")):
            lo, hi = bounds[axis]
            f.write(f"{lo:.10f} {hi:.10f} {name}lo {name}hi\n")
        f.write(f"\nMasses\n\n1 {AG_MASS:.8f} # Ag\n\nAtoms # atomic\n\n")
        for new_id, pos in enumerate(xyz, 1):
            f.write(f"{new_id} 1 {pos[0]:.10f} {pos[1]:.10f} {pos[2]:.10f}\n")

    id_map = OUT / "Ag_only_id_map.csv.gz"
    with gzip.open(id_map, "wt", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(("ag_only_id", "source_model_id", "x_A", "y_A", "z_A", "mask_signed_distance_A",
                         "fcc_i", "fcc_j", "fcc_k", "parent_x_A", "parent_y_A", "parent_z_A"))
        for new_id, (old_id, pos, phi, index, parent_pos) in enumerate(zip(ids, xyz, phi_atoms, lattice_indices, parent_xyz), 1):
            writer.writerow((new_id, int(old_id), *[f"{v:.10f}" for v in pos], f"{phi:.8f}",
                             *[int(v) for v in index], *[f"{v:.10f}" for v in parent_pos]))

    report = {
        "status": "AG_ONLY_QSGS_MASK_REBUILT_ON_EAM_FCC_LATTICE_NOT_COLD_VALIDATED",
        "source_data": str(SOURCE_DATA.relative_to(MODEL)),
        "source_data_sha256": sha256(SOURCE_DATA),
        "source_geometry_sha256": sha256(SOURCE_GEOMETRY),
        "source_morphology_sha256": sha256(SOURCE_MORPHOLOGY),
        "source_geometry_status": geom.get("status"),
        "qsgs_seed": int(maskdoc["seed"]),
        "qsgs_mask_grid": [int(x) for x in grid],
        "qsgs_mask_size_A": [float(x) for x in size],
        "qsgs_mask_voxel_size_A": [float(x) for x in spacing],
        "lattice_reconstruction": {
            "method": "preserve all Ag sites, remap to cubic Ag_u3.eam lattice, then apply the selected bounded in-plane orientation/translation",
            "Ag_u3_eam_equilibrium_a0_A": AG_EAM_A0_A,
            "fcc_cells_xyz": [int(x) for x in cells],
            "selected_orientation": orientation,
            "selected_translation_xy_A": [float(phase["dx_A"]), float(phase["dy_A"])],
            "phase_search_result": phase_path.name if phase_path.exists() else None,
            "parent_effective_a_xyz_A": [float(x) for x in parent_lattice_a],
            "parent_to_new_box_scale_xyz": ([float(AG_EAM_A0_A / parent_lattice_a[1]), float(AG_EAM_A0_A / parent_lattice_a[0]), 1.0]
                                             if orientation in ("rotate_z_90", "rotate_z_270")
                                             else [float(AG_EAM_A0_A / parent_lattice_a[0]), float(AG_EAM_A0_A / parent_lattice_a[1]), 1.0]),
            "fcc_parity_violations": parity_violations,
            "transformed_fcc_parity_violations": transformed_parity_violations,
            "max_parent_site_mapping_error_A": float(np.max(np.abs(lattice_error))),
            "atomic_occupancy_preserved": True,
            "Ag_atom_count_preserved": int(len(ids)),
        },
        "source_total_atoms": 916534,
        "exported_Ag_atoms": int(len(ids)),
        "excluded_TiSiC_atoms": 703734,
        "boundary_conditions_for_future_test": "p p f; all Ag above z=20 A mobile; bottom Ag z<=20 A fixed during the cold gate",
        "geometry_constraint": {
            "method": "smooth signed-distance mask derived only from the existing QSGS ceramic occupancy voxels",
            "wall_clearance_A": WALL_CLEARANCE_A,
            "wall_stiffness_eV_A2": WALL_STIFFNESS_EV_A2,
            "force_form": "F_geom = k * max(0, clearance - phi) * grad(phi)",
            "energy_form": "U_geom = 0.5 * k * max(0, clearance - phi)^2",
            "initial_atoms_with_wall_force": int(active.sum()),
            "initial_atoms_inside_skeleton_voxels": int(inside.sum()),
            "initial_phi_min_A": float(phi_atoms.min()),
            "initial_phi_p01_A": float(np.percentile(phi_atoms, 1)),
            "initial_phi_median_A": float(np.median(phi_atoms)),
            "initial_phi_p99_A": float(np.percentile(phi_atoms, 99)),
        },
        "interpretation_limits": [
            "Ti/Si/C atoms and all Ag-Ti/Si/C cross-potential forces are absent.",
            "The external signed-distance wall is an explicit geometric boundary condition; its stiffness is numerical and requires a bounded sensitivity check.",
            "The QSGS mask does not provide physical adhesion, wetting, interfacial heat transfer, or ceramic thermal response.",
            "No minimization, thermostat, MD, TTM, heating, or phase claim was run by this preparation script.",
        ],
        "artifacts": {"ag_data": ag_data.name, "id_map": id_map.name, "signed_distance_grid": wall_sdf.name},
    }
    (OUT / "geometry_mask_metrics.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
