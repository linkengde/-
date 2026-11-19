#!/usr/bin/env python3
"""Build center-force DFT crop candidates from the current 1M-atom MACE shortlist."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from ase import Atoms
from ase.io import write

ROOT = Path(__file__).resolve().parent
SOURCE = Path("/workspace/-/models/ag_ti3sic2_qsgs_1to1/ag_ti3sic2_qsgs.data.gz")
SHORTLIST = Path("/workspace/-/archive/mace-fit-ace-exploratory-20261119/followup_20261119/acquisition/label_design_center_shortlist.csv")
SCREEN = ROOT / "local_environment_screen.csv"
OUT_XYZ = ROOT / "qsgs_hotspot_candidates.extxyz"
OUT_JSON = ROOT / "qsgs_hotspot_candidates.json"
BOX = np.array([299.635, 301.857484616085, 258.6115], dtype=float)
TYPE_TO_SYMBOL = {1: "Ti", 2: "Si", 3: "C", 4: "Ag"}
RADIUS = 7.0
CELL = 26.0
SELECTED = [
    {"atom_id": 819467, "role": "training", "reason": "highest-force Ag center in the 81-center shortlist"},
    {"atom_id": 35127, "role": "training", "reason": "highest-force Ti center in the 81-center shortlist"},
    {"atom_id": 783134, "role": "training", "reason": "high-force Si center with Ti nearest neighbor"},
    {"atom_id": 558025, "role": "development_holdout", "reason": "highest-force Si center; held out to score the registered Si-center metric"},
]


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_lammps_data(path: Path):
    atom_count = None
    with gzip.open(path, "rt") as f:
        for line in f:
            fields = line.split()
            if len(fields) == 2 and fields[1] == "atoms":
                try:
                    atom_count = int(fields[0])
                    break
                except ValueError:
                    continue
    if atom_count is None:
        raise ValueError("could not read atom count")
    coords = np.zeros((atom_count + 1, 3), dtype=float)
    types = np.zeros(atom_count + 1, dtype=np.int8)
    maxid = 0
    in_atoms = False
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("Atoms #"):
                in_atoms = True
                continue
            if not in_atoms:
                continue
            fields = line.split()
            if len(fields) < 5:
                if maxid:
                    continue
                continue
            try:
                atom_id, typ = int(fields[0]), int(fields[1])
                position = [float(x) for x in fields[2:5]]
            except ValueError:
                if maxid:
                    break
                continue
            coords[atom_id] = position
            types[atom_id] = typ
            maxid = max(maxid, atom_id)
    valid = types > 0
    ids = np.flatnonzero(valid)
    return ids, coords[:maxid + 1], types[:maxid + 1]


def main() -> None:
    with SCREEN.open(newline="") as f:
        screen = {int(row["atom_id"]): row for row in csv.DictReader(f)}
    with SHORTLIST.open(encoding="utf-8-sig", newline="") as f:
        shortlist = {int(row["atom_id"]): row for row in csv.DictReader(f)}
    selected_ids = [int(row["atom_id"]) for row in SELECTED]
    if len(set(selected_ids)) != len(selected_ids) or any(i not in screen or i not in shortlist for i in selected_ids):
        raise SystemExit("a selected atom ID is duplicated or absent from the registered screening inputs")

    ids, coords, types = read_lammps_data(SOURCE)
    typ = types[ids]
    xyz = coords[ids]
    id_to_row = np.full(len(types), -1, dtype=np.int64)
    id_to_row[ids] = np.arange(len(ids))
    query_xyz = xyz.copy()
    query_xyz[:, 2] += BOX[2]
    tree = cKDTree(query_xyz, boxsize=[BOX[0], BOX[1], 3 * BOX[2]], compact_nodes=True, balanced_tree=True)

    frames = []
    records = []
    for item in SELECTED:
        center_id = int(item["atom_id"])
        row = int(id_to_row[center_id])
        if row < 0:
            raise SystemExit(f"center atom {center_id} is absent from the source geometry")
        center = xyz[row]
        dcoord = np.abs(center - np.array([float(shortlist[center_id][k]) for k in ("x_A", "y_A", "z_A")]))
        dcoord[0] = min(dcoord[0], BOX[0] - dcoord[0])
        dcoord[1] = min(dcoord[1], BOX[1] - dcoord[1])
        if float(dcoord.max()) >= 1e-3:
            raise SystemExit(f"center coordinate mismatch for {center_id}")
        qcenter = center.copy()
        qcenter[2] += BOX[2]
        query_rows = np.unique(ids[np.asarray(tree.query_ball_point(qcenter, r=RADIUS), dtype=np.int64)])
        delta = coords[query_rows] - center
        delta[:, 0] -= BOX[0] * np.rint(delta[:, 0] / BOX[0])
        delta[:, 1] -= BOX[1] * np.rint(delta[:, 1] / BOX[1])
        distance = np.linalg.norm(delta, axis=1)
        keep = (query_rows == center_id) | (distance <= RADIUS + 1e-8)
        local_ids = query_rows[keep]
        local_delta = delta[keep]
        local_distance = distance[keep]
        center_rows = np.flatnonzero(local_ids == center_id)
        if len(center_rows) != 1:
            raise SystemExit(f"selected center {center_id} does not occur exactly once")
        center_index = int(center_rows[0])
        if types[center_id] <= 0:
            raise SystemExit(f"center {center_id} has an invalid atom type")
        center_shell_count = int(np.count_nonzero((local_distance <= 6.0 + 1e-8) & (local_ids != center_id)))
        if center_shell_count < 1:
            raise SystemExit(f"center {center_id} has no neighbors in its 6 Å model cutoff")
        local_symbols = [TYPE_TO_SYMBOL[int(types[i])] for i in local_ids]
        atoms = Atoms(local_symbols, positions=local_delta + CELL / 2.0,
                      cell=np.eye(3) * CELL, pbc=True)
        atoms.new_array("source_atom_id", local_ids.astype(np.int64))
        center_mask = np.zeros(len(atoms), dtype=np.int8)
        center_mask[center_index] = 1
        atoms.new_array("is_center", center_mask)
        atoms.info.update({
            "config_type": f"MACE_1M_QSGS_hotspot_{center_id}",
            "source_model": "Ag/Ti3SiC2 1,021,793-atom QSGS MACE candidate",
            "source_data_sha256": sha(SOURCE),
            "center_atom_id": center_id,
            "center_species": local_symbols[center_index],
            "center_model_force_eV_A": float(screen[center_id]["model_force_eV_A"]),
            "crop_radius_A": RADIUS,
            "dft_target_policy": "center atom force only; crop energy and edge forces are excluded",
            "source_boundary_policy": "periodic x/y, open z",
            "crop_cell_A": CELL,
            "training_role": item["role"],
        })
        frames.append(atoms)
        records.append({
            **item,
            "center_species": local_symbols[center_index],
            "model_force_eV_A": float(screen[center_id]["model_force_eV_A"]),
            "nearest_neighbor_A": float(screen[center_id]["nearest_neighbor_A"]),
            "nearest_neighbor_element": screen[center_id]["nearest_neighbor_element"],
            "neighbors_within_6A": center_shell_count,
            "crop_atoms": len(atoms),
            "crop_formula": atoms.get_chemical_formula(),
            "center_image_gap_A": CELL - 2 * RADIUS,
            "source_coordinate_match": True,
            "label_role": item["role"],
        })
    write(OUT_XYZ, frames, format="extxyz")
    manifest = {
        "status": "geometry_candidate_crops_built_DFT_labels_pending",
        "source_data": str(SOURCE),
        "source_data_sha256": sha(SOURCE),
        "source_atoms": int(len(ids)),
        "source_boundary_policy": "periodic x/y, open z",
        "crop_radius_A": RADIUS,
        "model_cutoff_A": 6.0,
        "crop_cell_A": CELL,
        "center_image_gap_A": CELL - 2 * RADIUS,
        "DFT_candidate_method": "PBE/PW350/Gamma, center-force-only, finite crop energy and edge forces excluded",
        "training_centers": [x["atom_id"] for x in SELECTED if x["role"] == "training"],
        "development_centers": [x["atom_id"] for x in SELECTED if x["role"] == "development_holdout"],
        "development_is_blind_test": False,
        "records": records,
        "output_sha256": sha(OUT_XYZ),
        "center_species_counts": dict(Counter(x["center_species"] for x in records)),
    }
    OUT_JSON.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
