#!/usr/bin/env python3
"""Read-only checks for the provisional V17 training and development set."""
import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

import numpy as np
from ase.io import read

ENTRY = Path("research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v17_interface_energy")


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_rank(matrix):
    basis = []
    for source in matrix.tolist():
        row = [Fraction(int(x)) for x in source]
        for pivot, existing in basis:
            if row[pivot]:
                factor = row[pivot]
                row = [x - factor * y for x, y in zip(row, existing)]
        pivot = next((i for i, x in enumerate(row) if x), None)
        if pivot is not None:
            scale = row[pivot]
            basis.append((pivot, [x / scale for x in row]))
            basis.sort(key=lambda x: x[0])
    return len(basis)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", required=True, type=Path)
    args = parser.parse_args()
    repo = args.repo_root.resolve()
    data = repo / ENTRY / "data"
    manifest = json.loads((data / "dataset_manifest.json").read_text())
    require(manifest.get("status", "").startswith("PROVISIONAL"), "Dataset must remain marked provisional")
    require(manifest["roles"] == {"train": 43, "development_validation": 3, "test": 0}, "Wrong V17 split roles/counts")
    require(not (data / "test.extxyz").exists(), "V17 test split is not permitted in this screening run")
    require(sha(data / "train.extxyz") == manifest["train_input_sha256"], "Train hash mismatch")
    require(sha(data / "valid.extxyz") == manifest["development_validation_input_sha256"], "Dev hash mismatch")
    for line in (data / "SHA256SUMS.txt").read_text().splitlines():
        expected, name = line.split(maxsplit=1)
        require(sha(data / name.strip()) == expected, f"Dataset checksum mismatch: {name}")
    for relative, expected in manifest["source_files_sha256"].items():
        path = (repo / relative).resolve()
        require(path.is_relative_to(repo) and path.is_file() and sha(path) == expected, f"Source snapshot mismatch: {relative}")
    foundation = repo / "research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v12_interface_energy/foundation_models/mace-mp-0b3-medium.model"
    require(sha(foundation) == manifest["foundation_model_sha256"], "Pinned foundation model changed")
    train = read(data / "train.extxyz", index=":")
    valid = read(data / "valid.extxyz", index=":")
    require(len(train) == 43 and len(valid) == 3, "Serialized data count mismatch")
    elements = manifest["composition_matrix_elements"]
    matrix = np.asarray([[a.get_chemical_symbols().count(e) for e in elements] for a in train], dtype=int)
    require(matrix.tolist() == manifest["composition_matrix"] and exact_rank(matrix) == 4 == manifest["composition_matrix_rank_exact"], "Composition matrix/rank mismatch")
    for split, frames in (("train", train), ("valid", valid)):
        for i, atoms in enumerate(frames):
            require("REF_energy" in atoms.info and "REF_forces" in atoms.arrays, f"Missing labels: {split}[{i}]")
            energy = float(atoms.info["REF_energy"])
            forces = np.asarray(atoms.arrays["REF_forces"], dtype=float)
            require(np.isfinite(energy) and forces.shape == (len(atoms), 3) and np.isfinite(forces).all(), f"Invalid labels: {split}[{i}]")
            require(np.isfinite(atoms.positions).all() and np.isfinite(atoms.cell.array).all(), f"Invalid geometry: {split}[{i}]")
            require("v17_role" in atoms.info, f"Missing V17 role: {split}[{i}]")
            if split == "valid":
                require(len(atoms.info.get("v17_reference_sha256", "")) == 64, f"Missing dev reference hash: {i}")
    status = {}
    for atoms in train:
        key = str(atoms.info.get("v17_provenance_status", "unknown"))
        status[key] = status.get(key, 0) + 1
    require(status == manifest["train_provenance_status_counts"], "Per-frame provenance counts mismatch")
    require(len(manifest["verified_V17_pair_labels"]) == 8 and len(manifest["development_validation_labels"]) == 3, "Role label inventory mismatch")
    print(json.dumps({"status": "PROVISIONAL_PREFLIGHT_PASS", "train": len(train), "development_validation": len(valid), "test": 0, "composition_rank_exact": 4, "provenance_status_counts": status, "withheld_labels_opened": False}, indent=2))


if __name__ == "__main__":
    main()
