#!/usr/bin/env python3
"""Create an isolated 1% XY-strain/1% count-reduction mapping candidate."""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
TTM = HERE.parent
MODEL = TTM.parent
REPO_ROOT = MODEL.parents[3]
V6 = MODEL / "interface_reconstruction_v6_exploratory"
SOURCE_DATA = MODEL / "interface_reconstruction_v5/efs_static_v1_fixed_bottom/efs_minimized.data.gz"
SOURCE_DUMP = V6 / "nve_fixedids_after_nvt_final.dump.gz"
OUT = HERE / "candidate_xy202A_atoms_minus1pct"
OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(TTM))
sys.path.insert(0, str(HERE))
from agw_ttm_reference import AG_ATOMIC_VOLUME_A3, W_ATOMIC_VOLUME_A3, _raw_cic_volume_power  # noqa: E402
from diagnose_agw_ttm_mapping_v7 import read_dump, read_lammps_data_atoms  # noqa: E402

SCALE = 1.01
GRID = (20, 20, 16)
SPACING = (10.0, 10.0, 10.0)
OFFSET = (0.0, 2.5)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    ids, types, base_xyz = read_lammps_data_atoms(SOURCE_DATA)
    dump_ids, dump_types, xyz = read_dump(SOURCE_DUMP)
    if not (np.array_equal(ids, dump_ids) and np.array_equal(types, dump_types)):
        raise RuntimeError("Source data and frozen restart ID/type columns do not match")
    mapped = ((base_xyz[:, 2] >= 40.0 - 1e-9) & (base_xyz[:, 2] <= 300.0 + 1e-9)
              & (xyz[:, 2] >= 40.0 - 1e-9) & (xyz[:, 2] <= 200.0 + 1e-9))
    i_ag = np.flatnonzero(mapped & (types == 2))
    i_w = np.flatnonzero(mapped & (types == 1))
    xy = xyz[:, :2].copy()
    xy[:, 0] = (xy[:, 0] - OFFSET[0]) % 200.0
    xy[:, 1] = (xy[:, 1] - OFFSET[1]) % 200.0
    a = np.column_stack((xy[i_ag], xyz[i_ag, 2] - 40.0))
    w = np.column_stack((xy[i_w], xyz[i_w, 2] - 40.0))
    map_a, vol_a, _ = _raw_cic_volume_power(a, np.zeros(len(a)), GRID, SPACING, AG_ATOMIC_VOLUME_A3)
    map_w, vol_w, _ = _raw_cic_volume_power(w, np.zeros(len(w)), GRID, SPACING, W_ATOMIC_VOLUME_A3)
    total_fill = (vol_a + vol_w).ravel() / np.prod(SPACING)
    score_a = np.sum(map_a.weights * total_fill[map_a.cell_ids], axis=1)
    score_w = np.sum(map_w.weights * total_fill[map_w.cell_ids], axis=1)
    all_ag = int(np.count_nonzero(types == 2))
    all_w = int(np.count_nonzero(types == 1))
    n_remove_ag = int(round(all_ag * 0.01))
    n_remove_w = int(round(all_w * 0.01))
    keep_a = np.ones(len(i_ag), dtype=bool)
    keep_w = np.ones(len(i_w), dtype=bool)
    keep_a[np.argsort(score_a, kind="stable")[-n_remove_ag:]] = False
    keep_w[np.argsort(score_w, kind="stable")[-n_remove_w:]] = False
    remove_indices = np.concatenate((i_ag[~keep_a], i_w[~keep_w]))
    remove_ids = np.sort(ids[remove_indices])
    remove_ids_path = OUT / "removed_atom_ids.csv"
    with remove_ids_path.open("w", encoding="ascii", newline="\n") as f:
        f.write("atom_id,type\n")
        for atom_index in remove_indices:
            f.write(f"{int(ids[atom_index])},{int(types[atom_index])}\n")
    remove_set = set(map(int, remove_ids))
    new_ag = all_ag - n_remove_ag
    new_w = all_w - n_remove_w
    mass_ratio = new_ag * 107.8682 / (new_w * 183.84)
    if abs(mass_ratio - 1.0) > 1.0e-3:
        raise RuntimeError(f"Candidate mass ratio {mass_ratio:.9g} exceeds ±0.1%")

    base_scaled = base_xyz.copy()
    base_scaled[:, :2] *= SCALE
    data_out = OUT / "candidate.data.gz"
    with gzip.open(data_out, "wt", encoding="ascii", newline="\n") as f:
        f.write("Ag-W v8 isolated mapping sensitivity candidate; XY scale 1.01, 1% count removal\n\n")
        f.write(f"{len(ids) - len(remove_ids)} atoms\n2 atom types\n\n")
        f.write("0 202 xlo xhi\n0 202 ylo yhi\n0 300 zlo zhi\n\nMasses\n\n1 183.84\n2 107.8682\n\nAtoms # atomic\n\n")
        for atom_id, atom_type, pos in zip(ids, types, base_scaled):
            if int(atom_id) in remove_set:
                continue
            f.write(f"{int(atom_id)} {int(atom_type)} {pos[0]:.12g} {pos[1]:.12g} {pos[2]:.12g} 0 0 0\n")

    restart_out = OUT / "candidate.restart.dump.gz"
    atom_by_id = {}
    with gzip.open(SOURCE_DUMP, "rt", encoding="ascii") as f:
        lines = f.readlines()
    atom_header = next(i for i, line in enumerate(lines) if line.startswith("ITEM: ATOMS"))
    atom_cols = lines[atom_header].split()[2:]
    if atom_cols[:8] != ["id", "type", "x", "y", "z", "vx", "vy", "vz"]:
        raise RuntimeError(f"Unexpected restart columns: {atom_cols}")
    for line in lines[atom_header + 1:]:
        row = line.split()
        atom_by_id[int(row[0])] = row
    with gzip.open(restart_out, "wt", encoding="ascii", newline="\n") as f:
        f.write("ITEM: TIMESTEP\n4200\nITEM: NUMBER OF ATOMS\n")
        f.write(f"{len(ids) - len(remove_ids)}\nITEM: BOX BOUNDS pp pp ff\n")
        f.write("0 202\n0 202\n0 300\nITEM: ATOMS id type x y z vx vy vz\n")
        for atom_id, atom_type, _ in zip(ids, types, base_xyz):
            atom_id_int = int(atom_id)
            if atom_id_int in remove_set:
                continue
            row = atom_by_id[atom_id_int]
            x = float(row[2]) * SCALE
            y = float(row[3]) * SCALE
            f.write(f"{atom_id_int} {int(atom_type)} {x:.12g} {y:.12g} {row[4]} {row[5]} {row[6]} {row[7]}\n")

    manifest_out = OUT / "candidate_manifest.json"
    manifest = {
        "title": "Ag-W v8 isolated exploratory mapping candidate",
        "md_start_gate": "EXPLORATORY_ONLY",
        "static_audit_status": "STATIC_EFS_FAIL",
        "source_template_data_sha256": sha256(SOURCE_DATA),
        "source_restart_sha256": sha256(SOURCE_DUMP),
        "data_file": data_out.name,
        "data_sha256": sha256(data_out),
        "restart_dump": restart_out.name,
        "restart_dump_sha256": sha256(restart_out),
        "restart_dump_timestep": 4200,
        "candidate_box_xy_A": [202.0, 202.0],
        "candidate_atom_count": int(len(ids) - len(remove_ids)),
        "removed_atoms": {"Ag": n_remove_ag, "W": n_remove_w},
        "mass_ratio_Ag_over_W": float(mass_ratio),
        "mass_ratio_target_tolerance_fraction": 0.001,
        "position_edit": "uniform 1% xy affine dilation; z unchanged",
        "atom_removal_rule": "remove 1% of each full-system species from highest frozen-CIC local density score inside mapped Newton layer",
        "removed_atom_ids_file": remove_ids_path.name,
        "removed_atom_ids_sha256": sha256(remove_ids_path),
        "scope": "mapping-only preflight candidate; not a physical initial-state acceptance",
        "restrictions": [
            "Never overwrite the frozen v6 restart or v5 template.",
            "No NVT/NVE or laser source is authorized by this candidate manifest.",
            "Requires force, density, and 3D connectivity review before any MD interpretation.",
        ],
    }
    manifest_out.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "candidate_data": str(data_out), "candidate_restart": str(restart_out),
        "candidate_manifest": str(manifest_out), "removed_Ag": n_remove_ag,
        "removed_W": n_remove_w, "remaining_atoms": int(len(ids) - len(remove_ids)),
        "mass_ratio_Ag_over_W": float(mass_ratio), "data_sha256": manifest["data_sha256"],
        "restart_sha256": manifest["restart_dump_sha256"],
    }, indent=2))


if __name__ == "__main__":
    main()
