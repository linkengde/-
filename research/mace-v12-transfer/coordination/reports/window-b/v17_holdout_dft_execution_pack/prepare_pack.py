#!/usr/bin/env python3
"""Build a geometry-only manifest for the frozen V17 candidate input files."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from ase.io import read


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def geometry_identity(atoms) -> tuple[str, str, list[int]]:
    ids = np.asarray(atoms.arrays["lammps_id"], dtype=np.int64)
    payload = {
        "symbols_in_file_order": atoms.get_chemical_symbols(),
        "lammps_ids_in_file_order": ids.tolist(),
        # 1e-8 A bins make this identity stable across extxyz round trips.
        "positions_1e-8_A_integer": np.rint(atoms.positions * 1e8).astype(np.int64).tolist(),
        "cell_1e-8_A_integer": np.rint(atoms.cell.array * 1e8).astype(np.int64).tolist(),
        "pbc": np.asarray(atoms.pbc, dtype=bool).tolist(),
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ids_canonical = json.dumps(ids.tolist(), separators=(",", ":")).encode()
    return sha256(canonical).hexdigest(), sha256(ids_canonical).hexdigest(), ids.tolist()


def reject_labels(atoms, label: str) -> None:
    forbidden_info = {"energy", "free_energy", "ref_energy", "pw_pbe_energy_ev"}
    forbidden_arrays = {"forces", "ref_forces", "pw_pbe_forces"}
    found_info = sorted(key for key in atoms.info if key.lower() in forbidden_info)
    found_arrays = sorted(key for key in atoms.arrays if key.lower() in forbidden_arrays)
    if found_info or found_arrays:
        raise SystemExit(f"Candidate {label} contains label fields: {found_info} {found_arrays}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace-root", required=True, type=Path)
    parser.add_argument("--pack-root", required=True, type=Path)
    args = parser.parse_args()
    workspace = args.workspace_root.resolve()
    pack = args.pack_root.resolve()
    split_root = workspace / "research/mace-v12-transfer/coordination/reports/window-b/v17_split_geometry_design/final_set"
    role_path = split_root / "role_manifest.json"
    screening_path = split_root / "screening_manifest.json"
    source_inventory_path = split_root / "source_inventory.json"
    roles = json.loads(role_path.read_text())

    records = []
    geometry_sources: dict[str, dict[str, str]] = {}
    for candidate in roles["candidates"]:
        label = candidate["candidate_id"]
        source = split_root / candidate["file"]
        target = pack / "inputs" / source.name
        if not source.is_file() or not target.is_file():
            raise SystemExit(f"Missing candidate input: {label}")
        observed_hash = file_hash(source)
        if observed_hash != candidate["file_sha256"] or file_hash(target) != observed_hash:
            raise SystemExit(f"Candidate byte/hash mismatch: {label}")
        for path_key, hash_key in (("direct_source_path", "direct_source_sha256"), ("parent_input_path", "parent_input_sha256")):
            source_path = workspace / candidate[path_key]
            source_hash = file_hash(source_path) if source_path.is_file() else "missing"
            if source_hash != candidate[hash_key]:
                raise SystemExit(f"Geometry provenance hash mismatch: {label} {path_key}")
            geometry_sources[str(source_path.relative_to(workspace))] = {
                "sha256": source_hash,
                "role": "geometry_only_parent_input; no DFT result file read",
            }
        atoms = read(target, format="extxyz")
        reject_labels(atoms, label)
        if len(atoms) != candidate["atom_count"]:
            raise SystemExit(f"Atom-count mismatch: {label}")
        if atoms.get_chemical_formula() != candidate["formula"]:
            raise SystemExit(f"Formula mismatch: {label}")
        if np.asarray(atoms.pbc, dtype=bool).tolist() != candidate["pbc"]:
            raise SystemExit(f"PBC mismatch: {label}")
        if not np.allclose(atoms.cell.array, candidate["cell_A"], atol=1e-10, rtol=0):
            raise SystemExit(f"Cell mismatch: {label}")
        identity_hash, ids_hash, atom_ids = geometry_identity(atoms)
        role = candidate["proposed_role"]
        proposed_owner = "window-a" if role == "development_validation" else "window-b"
        screening = candidate["screening"]
        record = {
            "label": label,
            "interface": candidate["interface"],
            "role": role,
            "role_status": "proposed_pending_A_role_review_and_explicit_owner_assignment",
            "proposed_owner": proposed_owner,
            "owner_assignment_status": "proposal_only_not_claimed",
            "launch_status": "blocked_pending_A_role_review_and_owner_assignment",
            "input": f"inputs/{target.name}",
            "input_sha256": observed_hash,
            "split_geometry_file_sha256": candidate["file_sha256"],
            "split_geometry_identity_sha256": candidate["geometry_sha256"],
            "pack_geometry_identity_sha256": identity_hash,
            "ordered_atom_ids_sha256": ids_hash,
            "atom_ids_in_file_order": atom_ids,
            "atom_count": len(atoms),
            "formula": atoms.get_chemical_formula(),
            "symbols_in_file_order": atoms.get_chemical_symbols(),
            "cell_A": np.asarray(atoms.cell.array, dtype=float).tolist(),
            "pbc": np.asarray(atoms.pbc, dtype=bool).tolist(),
            "marked_pair": candidate["marked_pair"],
            "direct_source_path": candidate["direct_source_path"],
            "direct_source_sha256": candidate["direct_source_sha256"],
            "parent_input_path": candidate["parent_input_path"],
            "parent_input_sha256": candidate["parent_input_sha256"],
            "lineage": candidate["lineage"],
            "species_pair_minima_A": candidate["species_pair_minima_A"],
            "geometry_screen": {
                "decision": candidate["decision"],
                "exact_or_near_hits": screening["exact_or_near_hits"],
                "closest_inventory_geometry": screening["nearest_inventory_geometry"],
                "contact_floor_rule": screening["contact_floor_rule"],
                "species_contact_checks": screening["species_contact_checks"],
            },
            "future_test_access_policy": (
                "A records a frozen model/checkpoint-selection hash and timestamp before reading this label; "
                "B may perform only explicitly assigned DFT/archive verification and must not inspect or score test labels."
                if role == "withheld_test" else
                "Use for development validation only after A accepts the role and assigns execution."
            ),
        }
        records.append(record)

    if len(records) != 9:
        raise SystemExit(f"Expected 9 candidate inputs, found {len(records)}")
    if sum(row["role"] == "development_validation" for row in records) != 3:
        raise SystemExit("Expected 3 development-validation candidates")
    if sum(row["role"] == "withheld_test" for row in records) != 6:
        raise SystemExit("Expected 6 withheld-test candidates")

    manifest = {
        "task": "v17_holdout_dft_execution_pack",
        "prepared_by_owner_instance": "c5035b48-f43b-4da0-b8e4-e2862817f86a",
        "status": "preparation_complete_no_DFT_started",
        "authorization": "none; all launches blocked_pending_A_role_review_and_owner_assignment",
        "source_geometry_bundle": {
            "role_manifest": str(role_path.relative_to(workspace)),
            "role_manifest_sha256": file_hash(role_path),
            "screening_manifest": str(screening_path.relative_to(workspace)),
            "screening_manifest_sha256": file_hash(screening_path),
            "source_inventory": str(source_inventory_path.relative_to(workspace)),
            "source_inventory_sha256": file_hash(source_inventory_path),
            "audit_report": str((split_root / "audit_report.md").relative_to(workspace)),
            "audit_report_sha256": file_hash(split_root / "audit_report.md"),
            "candidate_bundle_sha256_list": str((split_root / "SHA256SUMS.txt").relative_to(workspace)),
            "candidate_bundle_sha256_list_sha256": file_hash(split_root / "SHA256SUMS.txt"),
            "source_dataset_coordinate_hashes_are_provenance_only_not_reread": True,
        },
        "proposed_owner_split": {
            "window-a": [row["label"] for row in records if row["role"] == "development_validation"],
            "window-b": [row["label"] for row in records if row["role"] == "withheld_test"],
            "status": "proposal_only_A_review_and_assignment_required",
        },
        "method": {
            "calculator": "GPAW 26.7.0",
            "xc": "PBE",
            "basis": "plane wave",
            "cutoff_eV": 500,
            "kpts": [1, 1, 1],
            "smearing": "FermiDirac",
            "smearing_eV": 0.1,
            "mixer": {"type": "Pulay", "beta": 0.05, "nmaxold": 8, "weight": 100},
            "scf_convergence": {"energy": 1e-5, "density": 1e-5, "eigenstates": 1e-5},
            "maxiter": 160,
            "mpi_ranks": 4,
            "threads_per_rank": 1,
            "energy_convention": "GPAW native extrapolated energy (force_consistent=False); free energy recorded separately when supported",
            "geometry": "fixed-geometry single point; no relaxation",
        },
        "blind_policy": {
            "withheld_test_never_selects_model_or_parameters": True,
            "A_freezes_selected_checkpoint_before_reading_withheld_test_labels": True,
            "B_must_not_analyze_or_model_score_withheld_test_labels": True,
        },
        "records": records,
    }
    (pack / "input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    inventory = {
        "task": "v17_holdout_dft_execution_pack",
        "scope": "candidate extxyz and published geometry-only parent input files; no calculated result paths",
        "source_geometry_bundle": manifest["source_geometry_bundle"],
        "candidate_inputs": [
            {"label": row["label"], "path": row["input"], "sha256": row["input_sha256"], "role": row["role"]}
            for row in records
        ],
        "geometry_parent_inputs": [
            {"path": path, **metadata} for path, metadata in sorted(geometry_sources.items())
        ],
        "provenance_only_hashes_not_reread": [
            {"path": row["lineage"].get("source_dataset_path"), "sha256": row["lineage"].get("source_dataset_sha256")}
            for row in records
        ],
    }
    (pack / "source_inventory.json").write_text(json.dumps(inventory, indent=2) + "\n")
    print(f"Prepared {len(records)} geometry-only manifest records; no GPAW import or calculation.")


if __name__ == "__main__":
    main()
