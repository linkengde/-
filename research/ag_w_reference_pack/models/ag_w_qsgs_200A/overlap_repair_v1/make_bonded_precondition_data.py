#!/usr/bin/env python3
"""Create a temporary harmonic-neighbor cage for Ag-W overlap preconditioning.

This file is an exploratory geometry tool. Its harmonic bonds are not a physical
Ag-W potential and must not be used for production MD or reported energies.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "ag_w_qsgs.data.gz"
OUTPUT = Path(__file__).with_name("bonded_precondition.data.gz")
BOX = np.array([200.0, 200.0, 1200.0])
Z_SHIFT = 500.0
NEIGHBOR_CUTOFF = {1: 2.95, 2: 3.15}
R0 = {1: 2.745, 2: 2.886}


def read_data(path: Path) -> np.ndarray:
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
    data = np.asarray(rows, dtype=float)
    return data[np.argsort(data[:, 0])]


def main() -> None:
    atoms = read_data(SOURCE)
    coords = atoms[:, 2:5].copy()
    tree_coords = coords.copy()
    tree_coords[:, 2] += Z_SHIFT
    bond_groups = []
    counts = {}
    for typ, cutoff in NEIGHBOR_CUTOFF.items():
        indices = np.flatnonzero(atoms[:, 1] == typ)
        pairs = cKDTree(tree_coords[indices], boxsize=BOX).query_pairs(
            cutoff, output_type="ndarray"
        )
        left, right = indices[pairs[:, 0]], indices[pairs[:, 1]]
        # Do not write cross-box bonds: LAMMPS image flags in the starting
        # atomic data do not encode a globally consistent unwrapped graph.
        delta = coords[right] - coords[left]
        interior = (np.abs(delta[:, 0]) < 100.0) & (np.abs(delta[:, 1]) < 100.0)
        left, right = left[interior], right[interior]
        counts[typ] = int(len(left))
        group = np.empty((len(left), 3), dtype=np.int64)
        group[:, 0] = atoms[left, 0].astype(np.int64)
        group[:, 1] = atoms[right, 0].astype(np.int64)
        group[:, 2] = typ
        bond_groups.append(group)
    bonds = np.concatenate(bond_groups, axis=0)

    with gzip.open(OUTPUT, "wt", compresslevel=2) as stream:
        stream.write("Temporary Ag-W overlap-preconditioning structure; NOT a physical model\n\n")
        stream.write(f"{len(atoms)} atoms\n{len(bonds)} bonds\n\n2 atom types\n2 bond types\n\n")
        stream.write("0.0 200.0 xlo xhi\n0.0 200.0 ylo yhi\n0.0 300.0 zlo zhi\n\n")
        stream.write("Masses\n\n1 183.83999995270821\n2 107.86819997225152\n\n")
        stream.write("Bond Coeffs\n\n1 10.0 2.745\n2 10.0 2.886\n\n")
        stream.write("Atoms # bond\n\n")
        for atom_id, typ, x, y, z in atoms:
            stream.write(f"{int(atom_id)} 1 {int(typ)} {x:.15g} {y:.15g} {z:.15g}\n")
        stream.write("\nBonds\n\n")
        for bond_id, (atom_i, atom_j, typ) in enumerate(bonds, 1):
            stream.write(f"{bond_id} {typ} {atom_i} {atom_j}\n")
    print({"atoms": len(atoms), "bonds_by_type": counts, "total_bonds": len(bonds), "output": str(OUTPUT)})


if __name__ == "__main__":
    main()
