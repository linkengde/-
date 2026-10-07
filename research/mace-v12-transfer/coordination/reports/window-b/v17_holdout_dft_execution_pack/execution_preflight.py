#!/usr/bin/env python3
"""Geometry-only and authorization guards for a future, explicitly assigned run."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from ase.io import read


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def geometry_identity(atoms) -> str:
    payload = {
        "symbols_in_file_order": atoms.get_chemical_symbols(),
        "lammps_ids_in_file_order": np.asarray(atoms.arrays["lammps_id"], dtype=np.int64).tolist(),
        "positions_1e-8_A_integer": np.rint(atoms.positions * 1e8).astype(np.int64).tolist(),
        "cell_1e-8_A_integer": np.rint(atoms.cell.array * 1e8).astype(np.int64).tolist(),
        "pbc": np.asarray(atoms.pbc, dtype=bool).tolist(),
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return sha256(canonical).hexdigest()


def no_input_labels(atoms, label: str) -> None:
    info_keys = {"energy", "free_energy", "ref_energy", "pw_pbe_energy_ev"}
    array_keys = {"forces", "ref_forces", "pw_pbe_forces"}
    bad_info = sorted(key for key in atoms.info if key.lower() in info_keys)
    bad_arrays = sorted(key for key in atoms.arrays if key.lower() in array_keys)
    if bad_info or bad_arrays:
        raise SystemExit(f"{label}: geometry input unexpectedly contains labels")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--workspace-root", required=True, type=Path)
    p.add_argument("--pack-root", required=True, type=Path)
    p.add_argument("--input-root", required=True, type=Path)
    p.add_argument("--run-root", required=True, type=Path)
    p.add_argument("--archive-root", required=True, type=Path)
    p.add_argument("--authorization", required=True, type=Path)
    p.add_argument("--owner", required=True, choices=("window-a", "window-b"))
    p.add_argument("--labels", required=True, nargs="+")
    args = p.parse_args()

    workspace = args.workspace_root.resolve()
    pack = args.pack_root.resolve()
    input_root = args.input_root.resolve()
    run_root = args.run_root.resolve()
    archive_root = args.archive_root.resolve()
    auth_path = args.authorization.resolve()
    manifest = json.loads((pack / "input_manifest.json").read_text())
    records = {row["label"]: row for row in manifest["records"]}
    if len(set(args.labels)) != len(args.labels):
        raise SystemExit("Duplicate label requested")
    unknown = sorted(set(args.labels) - set(records))
    if unknown:
        raise SystemExit(f"Unknown labels: {unknown}")
    if not auth_path.is_file():
        raise SystemExit("BLOCKED: A authorization file is absent; do not create one in this B package")
    authorization = json.loads(auth_path.read_text())
    if authorization.get("status") != "authorized":
        raise SystemExit("BLOCKED: A authorization status is not authorized")
    if authorization.get("reviewed_by") != "window-a":
        raise SystemExit("BLOCKED: roles have not been reviewed by A")
    expected_role_hash = manifest["source_geometry_bundle"]["role_manifest_sha256"]
    if authorization.get("role_manifest_sha256") != expected_role_hash:
        raise SystemExit("BLOCKED: authorization references a different role manifest")

    for name in ("role_manifest", "screening_manifest", "source_inventory", "audit_report", "candidate_bundle_sha256_list"):
        item = manifest["source_geometry_bundle"]
        path = workspace / item[name]
        expected = item[name + "_sha256"]
        if not path.is_file() or digest(path) != expected:
            raise SystemExit(f"Source bundle hash mismatch: {name}")

    assignment = authorization.get("owner_assignments", {})
    selected_model = authorization.get("frozen_selection", {})
    for label in args.labels:
        rec = records[label]
        if assignment.get(label) != args.owner:
            raise SystemExit(f"BLOCKED: A has not assigned {label} to {args.owner}")
        if authorization.get("role_assignments", {}).get(label) != rec["role"]:
            raise SystemExit(f"BLOCKED: A has not accepted the role for {label}")
        expected_input_hash = authorization.get("input_sha256", {}).get(label)
        if expected_input_hash != rec["input_sha256"]:
            raise SystemExit(f"BLOCKED: authorization input hash mismatch for {label}")
        if rec["role"] == "withheld_test":
            if selected_model.get("frozen_before_test_label_read") is not True:
                raise SystemExit("BLOCKED: A has not recorded frozen checkpoint selection before test-label access")
            for key in ("selected_model_sha256", "selection_record_sha256", "frozen_at_utc"):
                if not selected_model.get(key):
                    raise SystemExit(f"BLOCKED: frozen_selection.{key} is required for withheld tests")

        geometry = input_root / Path(rec["input"]).name
        if not geometry.is_file() or digest(geometry) != rec["input_sha256"]:
            raise SystemExit(f"Input file hash mismatch for {label}")
        atoms = read(geometry, format="extxyz")
        no_input_labels(atoms, label)
        if geometry_identity(atoms) != rec["pack_geometry_identity_sha256"]:
            raise SystemExit(f"Ordered geometry identity mismatch for {label}")
        if len(atoms) != rec["atom_count"] or atoms.get_chemical_formula() != rec["formula"]:
            raise SystemExit(f"Atom-count/formula mismatch for {label}")
        if atoms.get_chemical_symbols() != rec["symbols_in_file_order"]:
            raise SystemExit(f"Element order mismatch for {label}")
        if np.asarray(atoms.arrays["lammps_id"], dtype=int).tolist() != rec["atom_ids_in_file_order"]:
            raise SystemExit(f"Atom-ID order mismatch for {label}")
        if not np.allclose(atoms.cell.array, rec["cell_A"], atol=1e-10, rtol=0):
            raise SystemExit(f"Cell mismatch for {label}")
        if np.asarray(atoms.pbc, dtype=bool).tolist() != rec["pbc"]:
            raise SystemExit(f"PBC mismatch for {label}")
        for source_key in ("direct_source_path", "parent_input_path"):
            source = workspace / rec[source_key]
            hash_key = "direct_source_sha256" if source_key == "direct_source_path" else "parent_input_sha256"
            if not source.is_file() or digest(source) != rec[hash_key]:
                raise SystemExit(f"Geometry provenance source hash mismatch for {label}: {source_key}")

        # Never overwrite or resume a prior directory implicitly.
        for root, kind in ((run_root, "run"), (archive_root, "archive")):
            target = root / label
            if target.exists():
                raise SystemExit(f"Refusing existing {kind} path; inspect and assign recovery explicitly: {target}")
            if list(root.glob(f".{label}.staging-*")) if root.exists() else []:
                raise SystemExit(f"Refusing stale {kind} staging path for {label}")

    print("PRE-FLIGHT PASS: A authorization, source hashes and geometry identities verified.")


if __name__ == "__main__":
    main()
