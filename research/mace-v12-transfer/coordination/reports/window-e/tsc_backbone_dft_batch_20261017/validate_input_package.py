#!/usr/bin/env python3
"""Validate frozen structure candidates and metadata only; performs no scientific calculation."""
from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
from ase.io import read


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coordinate_sha(atoms) -> str:
    h = hashlib.sha256()
    h.update(" ".join(atoms.get_chemical_symbols()).encode())
    h.update(np.asarray(atoms.cell, dtype="<f8").tobytes(order="C"))
    # Match the 8-decimal Cartesian precision serialized by ASE extxyz.
    h.update(np.asarray(np.round(atoms.positions, 8), dtype="<f8").tobytes(order="C"))
    h.update(np.asarray(atoms.pbc, dtype=np.uint8).tobytes())
    h.update(np.asarray(atoms.arrays["source_atom_id"], dtype="<i8").tobytes())
    return h.hexdigest()


def min_distance(atoms) -> float:
    d = atoms.get_all_distances(mic=True)
    np.fill_diagonal(d, np.inf)
    return float(np.min(d))


def main() -> None:
    base = Path(__file__).resolve().parent
    manifest = json.loads((base / "SURFACE_CANDIDATE_MANIFEST.json").read_text())
    tasks = json.loads((base / "PUBLIC_TASK_MANIFEST.json").read_text())
    kpoints = json.loads((base / "KPOINT_REVIEW.json").read_text())
    entries = manifest["entries"]
    checks = []
    atoms_by_id = {}
    for entry in entries:
        path = base / entry["path"]
        atoms = read(path, format="extxyz")
        symbols = Counter(atoms.get_chemical_symbols())
        source_ids = atoms.arrays.get("source_atom_id")
        valid = {
            "file_sha_matches": sha256(path) == entry["file_sha256"],
            "coordinate_sha_matches": coordinate_sha(atoms) == entry["coordinate_sha256"],
            "atom_count_matches": len(atoms) == entry["atom_count"] == 100,
            "composition_matches": symbols == Counter(entry["composition"]),
            "pbc_is_ttt": atoms.pbc.tolist() == [True, True, True],
            "atom_ids_are_stable_1_to_N": source_ids is not None and np.array_equal(source_ids, np.arange(1, 101)),
            "no_dft_labels": "energy" not in atoms.info and "free_energy" not in atoms.info and "forces" not in atoms.arrays and "stress" not in atoms.info,
            "minimum_distance_ge_1p5A": min_distance(atoms) >= 1.5,
        }
        if not all(valid.values()):
            raise RuntimeError(f"Candidate failed validation {entry['structure_id']}: {valid}")
        atoms_by_id[entry["structure_id"]] = atoms
        checks.append({"structure_id": entry["structure_id"], "input_sha256": sha256(path), "checks": valid})

    pair_checks = []
    for termination in ("TI", "SI"):
        prefix = f"TSC100_{termination}0001_COD9009647_2x2x2"
        ref = atoms_by_id[prefix + "_REF"]
        resp = atoms_by_id[prefix + "_RESP_Z020"]
        delta = resp.positions - ref.positions
        changed = np.flatnonzero(np.linalg.norm(delta, axis=1) > 1e-9)
        expected = np.asarray([1, 25] if termination == "TI" else [5, 24]) - 1
        if not np.array_equal(changed, expected):
            raise RuntimeError(f"Unexpected changed atom IDs in {termination} response: {(changed+1).tolist()}")
        expected_vectors = np.zeros((2, 3))
        expected_vectors[:, 2] = [-0.02, 0.02]
        if not np.allclose(delta[expected], expected_vectors, atol=1e-12, rtol=0):
            raise RuntimeError(f"Unexpected response displacement for {termination}")
        if not np.array_equal(ref.arrays["source_atom_id"], resp.arrays["source_atom_id"]):
            raise RuntimeError("Atom IDs changed within paired family")
        if not np.allclose(ref.cell, resp.cell, atol=0, rtol=0) or ref.pbc.tolist() != resp.pbc.tolist():
            raise RuntimeError("Cell/PBC changed within paired family")
        pair_checks.append({
            "family": prefix,
            "only_changed_atom_ids": (changed + 1).tolist(),
            "displacement_vectors_angstrom": delta[changed].tolist(),
            "atom_order_cell_pbc_unchanged": True,
        })
    if not kpoints["surface_candidate"]["in_plane_grid_exactly_matches_parent"]:
        raise RuntimeError("Surface MP grid does not match parent in-plane mesh")
    if len(tasks["tasks"]) != 6 or tasks["tasks"][4]["status"] != "HOLD_NO_VERIFIED_INDEPENDENT_NONTRAIN_SOURCE":
        raise RuntimeError("Public task manifest scope/status mismatch")

    result = {
        "status": "PASS_GEOMETRY_AND_HASH_CHECKS; NO_DFT_RUN; DEV_HOLD",
        "surface_candidate_count": len(entries),
        "task_count_in_public_manifest": len(tasks["tasks"]),
        "C_input_review_required": True,
        "execution_authorized": False,
        "checks": checks,
        "paired_geometry_checks": pair_checks,
        "kpoint_status": kpoints["status"],
        "dev_status": "HOLD_NO_VERIFIED_INDEPENDENT_NONTRAIN_SOURCE",
    }
    (base / "INPUT_VALIDATION.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "candidate_count": len(entries), "pairs": len(pair_checks)}, indent=2))


if __name__ == "__main__":
    main()
