#!/usr/bin/env python3
"""Read a GPAW checkpoint without evaluating energy, forces, or starting SCF."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read
from gpaw import GPAW
from gpaw.mpi import world


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def geometry_sha(atoms) -> str:
    h = hashlib.sha256()
    h.update(" ".join(atoms.get_chemical_symbols()).encode())
    h.update(np.asarray(atoms.cell, dtype="<f8").tobytes(order="C"))
    h.update(np.asarray(np.round(atoms.positions, 8), dtype="<f8").tobytes(order="C"))
    h.update(np.asarray(atoms.pbc, dtype=np.uint8).tobytes())
    h.update(np.asarray(atoms.arrays["source_atom_id"], dtype="<i8").tobytes())
    return h.hexdigest()


parser = argparse.ArgumentParser()
parser.add_argument("checkpoint", type=Path)
parser.add_argument("input", type=Path)
parser.add_argument("expected_geometry_sha256")
args = parser.parse_args()
checkpoint = args.checkpoint.resolve()
input_path = args.input.resolve()
expected = read(input_path, format="extxyz")
calc = GPAW(str(checkpoint), txt=None)
atoms = calc.atoms.copy()
# Force loading of saved wavefunctions. This does not run SCF or evaluate labels.
calc.dft.ibzwfs.make_sure_wfs_are_read_from_gpw_file()
checks = {
    "atom_count": len(atoms) == len(expected),
    "symbols_order": atoms.get_chemical_symbols() == expected.get_chemical_symbols(),
    "positions_exact": np.array_equal(np.asarray(atoms.positions), np.asarray(expected.positions)),
    "cell_exact": np.array_equal(np.asarray(atoms.cell), np.asarray(expected.cell)),
    "pbc_exact": np.array_equal(np.asarray(atoms.pbc), np.asarray(expected.pbc)),
}
if not all(checks.values()):
    raise RuntimeError(f"checkpoint geometry mismatch: {checks}")
ids_missing = "source_atom_id" not in atoms.arrays
if ids_missing:
    atoms.set_array("source_atom_id", expected.arrays["source_atom_id"].copy())
else:
    checks["source_atom_ids_exact"] = np.array_equal(
        atoms.arrays["source_atom_id"], expected.arrays["source_atom_id"])
if geometry_sha(atoms) != args.expected_geometry_sha256:
    raise RuntimeError("checkpoint geometry and atom-ID mapping SHA mismatch")
checks["geometry_sha_matches"] = True
if world.rank == 0:
    print(json.dumps({
        "status": "READBACK_PASS",
        "checkpoint_sha256": sha256(checkpoint),
        "checkpoint_bytes": checkpoint.stat().st_size,
        "input_sha256": sha256(input_path),
        "geometry_sha256": args.expected_geometry_sha256,
        "atom_count": len(atoms),
        "checks": checks,
        "source_atom_ids_restored_from_frozen_input": ids_missing,
        "wavefunctions": "loaded from GPW",
        "scf_started": False,
        "energy_forces_stress_read": False,
    }, indent=2, sort_keys=True))
world.barrier()
