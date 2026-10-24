#!/usr/bin/env python3
"""Independently audit the six COD-derived Ti3SiC2 geometry proposals.

This checks source/file hashes, exact parent-supercell geometry, composition,
periodicity, finite coordinates, minimum distances, label-free status, and
that displaced variants preserve their parent cell and ordered atom identities.
It does not assign dataset roles or enable DFT.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.build import make_supercell
from ase.geometry import find_mic
from ase.io import read
from ase.neighborlist import neighbor_list


PACKAGE = Path(__file__).resolve().parent
REPO = PACKAGE.parents[5]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_geometry_sha256(atoms) -> str:
    payload = {
        "numbers": atoms.numbers.tolist(),
        "cell_A_rounded_8dp": np.round(atoms.cell.array, 8).tolist(),
        "scaled_positions_wrapped_rounded_8dp": np.round(
            atoms.get_scaled_positions(wrap=True), 8
        ).tolist(),
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def minimum_distinct_atom_distance(atoms) -> float:
    i, j, distances = neighbor_list("ijd", atoms, 4.0)
    values = [float(d) for d, a, b in zip(distances, i, j) if a != b]
    if not values:
        raise ValueError("No distinct atom pairs found inside 4 Å cutoff")
    return min(values)


def audit() -> dict:
    manifest = json.loads((PACKAGE / "candidate_manifest.json").read_text())
    source = REPO / manifest["source_path"]
    source_hash = sha256(source)
    if source_hash != manifest["source_cif_sha256"]:
        raise ValueError("Source CIF hash differs from the reviewed manifest")
    source_atoms = read(source)
    expected_source_formula = "C4Si2Ti6"
    if len(source_atoms) != 12 or source_atoms.get_chemical_formula() != expected_source_formula:
        raise ValueError("Unexpected COD source CIF contents")

    parents = {}
    records = []
    for expected in manifest["records"]:
        path = PACKAGE / expected["file"]
        file_hash = sha256(path)
        if file_hash != expected["sha256"]:
            raise ValueError(f"Candidate hash mismatch: {expected['file']}")
        atoms = read(path)
        if len(atoms) != expected["atoms"]:
            raise ValueError(f"Atom-count mismatch: {expected['file']}")
        if atoms.get_chemical_formula() != expected["formula"]:
            raise ValueError(f"Formula mismatch: {expected['file']}")
        if not np.asarray(atoms.pbc, dtype=bool).all():
            raise ValueError(f"Candidate is not fully periodic: {expected['file']}")
        if not np.isfinite(atoms.positions).all() or not np.isfinite(atoms.cell.array).all():
            raise ValueError(f"Non-finite geometry: {expected['file']}")
        if any(k.lower() in {"energy", "free_energy", "forces", "force"} for k in atoms.info):
            raise ValueError(f"Label-like info field found: {expected['file']}")
        if any(k.lower() in {"energy", "free_energy", "forces", "force"} for k in atoms.arrays):
            raise ValueError(f"Label-like array found: {expected['file']}")
        if str(atoms.info.get("dft_launch_enabled", "")).lower() not in {"false", "0"}:
            raise ValueError(f"DFT launch is not explicitly disabled: {expected['file']}")
        role = str(atoms.info.get("dataset_role", ""))
        if "geometry_only" not in role or "no_DFT_label" not in role:
            raise ValueError(f"Unexpected dataset-role metadata: {expected['file']}")
        if atoms.info.get("proposal_only_no_DFT_label") is not True:
            raise ValueError(f"Missing no-label proposal flag: {expected['file']}")
        if atoms.info.get("source_cif_sha256") != manifest["source_cif_sha256"]:
            raise ValueError(f"Candidate source hash metadata mismatch: {expected['file']}")
        if atoms.info.get("source_accession") != manifest["source_accession"]:
            raise ValueError(f"Candidate source accession mismatch: {expected['file']}")
        if atoms.info.get("source_doi") != manifest["source_doi"]:
            raise ValueError(f"Candidate source DOI mismatch: {expected['file']}")
        if atoms.info.get("parent_family_id") != expected["parent_family_id"]:
            raise ValueError(f"Parent-family metadata mismatch: {expected['file']}")
        if atoms.info.get("proposal_kind") != expected["proposal_kind"]:
            raise ValueError(f"Proposal-kind metadata mismatch: {expected['file']}")
        expected_seed = -1 if expected["seed"] is None else expected["seed"]
        if int(atoms.info.get("proposal_seed", -999)) != expected_seed:
            raise ValueError(f"Proposal-seed metadata mismatch: {expected['file']}")
        if not np.isclose(float(atoms.info.get("displacement_sigma_A", np.nan)), expected["displacement_sigma_A"]):
            raise ValueError(f"Displacement-scale metadata mismatch: {expected['file']}")

        matrix = np.asarray(expected["supercell_matrix"], dtype=int)
        expected_atoms = make_supercell(source_atoms, matrix, wrap=True)
        minimum_distance = minimum_distinct_atom_distance(atoms)
        if minimum_distance <= 1.8:
            raise ValueError(f"Minimum-distance screen failed: {expected['file']}")

        record = {
            "file": expected["file"],
            "sha256": file_hash,
            "canonical_geometry_sha256": canonical_geometry_sha256(atoms),
            "formula": atoms.get_chemical_formula(),
            "atoms": len(atoms),
            "cell_volume_A3": float(atoms.get_volume()),
            "minimum_distinct_atom_distance_A": minimum_distance,
            "pbc": atoms.pbc.tolist(),
            "info_keys": sorted(atoms.info),
            "array_keys": sorted(atoms.arrays),
            "proposal_kind": expected["proposal_kind"],
            "seed": expected["seed"],
            "parent_family_id": expected["parent_family_id"],
            "lineage_group_id": manifest["lineage_group_id"],
            "dft_launch_enabled": False,
            "energy_or_force_labels_present": False,
        }

        if expected["proposal_kind"] == "unperturbed_supercell":
            if atoms.get_chemical_symbols() != expected_atoms.get_chemical_symbols():
                raise ValueError(f"Parent atom-order mismatch: {expected['file']}")
            mic, _ = find_mic(atoms.positions - expected_atoms.positions, atoms.cell, pbc=True)
            source_delta = np.linalg.norm(mic, axis=1)
            record["source_supercell_max_site_delta_A"] = float(source_delta.max())
            if record["source_supercell_max_site_delta_A"] > 1.0e-7:
                raise ValueError(f"Parent is not the declared COD supercell: {expected['file']}")
            parents[expected["parent_family_id"]] = atoms.copy()
        else:
            parent = parents.get(expected["parent_family_id"])
            if parent is None:
                raise ValueError(f"Parent must precede displacement: {expected['file']}")
            if atoms.get_chemical_symbols() != parent.get_chemical_symbols():
                raise ValueError(f"Displaced atom order differs from parent: {expected['file']}")
            if not np.allclose(atoms.cell.array, parent.cell.array, atol=1e-10, rtol=0.0):
                raise ValueError(f"Displaced cell differs from parent: {expected['file']}")
            mic, _ = find_mic(atoms.positions - parent.positions, parent.cell, pbc=True)
            displacement_norms = np.linalg.norm(mic, axis=1)
            record.update(
                {
                    "displacement_component_sigma_A": expected["displacement_sigma_A"],
                    "displacement_rms_vector_A": float(np.sqrt(np.mean(displacement_norms**2))),
                    "displacement_max_vector_A": float(displacement_norms.max()),
                }
            )
        records.append(record)

    if len(records) != 6:
        raise ValueError(f"Expected six reviewed candidates, found {len(records)}")
    if len({r["sha256"] for r in records}) != len(records):
        raise ValueError("Candidate file hashes are not unique")
    if len({r["canonical_geometry_sha256"] for r in records}) != len(records):
        raise ValueError("Candidate geometry hashes are not unique")
    if {r["lineage_group_id"] for r in records} != {"COD9009647_Ti3SiC2"}:
        raise ValueError("Candidate source-lineage grouping changed unexpectedly")

    return {
        "status": "INDEPENDENT_GEOMETRY_CHECK_PASS_ROLE_AND_DFT_HOLD",
        "scope": "Geometry/file/source-lineage audit only; not source-independent validation or dataset approval.",
        "source_accession": manifest["source_accession"],
        "source_doi": manifest["source_doi"],
        "source_cif_sha256": source_hash,
        "candidate_count": len(records),
        "unique_file_hashes": True,
        "unique_geometry_hashes": True,
        "all_candidates_share_one_source_lineage": True,
        "checks": {
            "manifest_file_hashes": "PASS",
            "COD_CIF_reconstruction_of_both_unperturbed_supercells": "PASS",
            "formula_atom_count_and_full_PBC": "PASS",
            "finite_coordinates_and_cells": "PASS",
            "minimum_distances_above_1p8_A": "PASS",
            "no_energy_or_force_labels": "PASS",
            "DFT_explicitly_disabled": "PASS",
            "displacements_preserve_parent_cell_and_atom_order": "PASS",
        },
        "role_decision": "Keep all six unassigned and in one COD9009647_Ti3SiC2 lineage group. The 72- and 96-atom cells and random displacement seeds do not create independent source families. Do not use these as independent validation or blind test structures.",
        "DFT_decision": "Do not launch a B1 label pilot from this review. The source/role gate remains unresolved, the environment task card restricts the active owner to registered V20 numerical controls, and these six geometries are not registered DFT inputs.",
        "limitations": [
            "The source DOI identifies the COD experimental structure lineage; this audit does not reproduce the original paper's experiment.",
            "Minimum-distance screening is only a geometry sanity check, not a stability or relaxed-structure test.",
            "The separately prepared Bai2010 DFT-derived candidate package is a different source lineage; its exact paper parameters remain pending full-text verification.",
        ],
        "records": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=PACKAGE / "independent_review.json")
    args = parser.parse_args()
    result = audit()
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"status": result["status"], "candidate_count": result["candidate_count"], "output": str(args.output)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
