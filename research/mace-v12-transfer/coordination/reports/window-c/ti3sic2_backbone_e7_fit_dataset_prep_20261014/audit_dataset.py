#!/usr/bin/env python3
"""Verify the Stage C geometry import and emit a read-only dataset audit.

This is a geometry/provenance audit only. It never instantiates a calculator,
reads or writes a potential, performs DFT, or optimizes parameters.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
from ase.io import read


ROOT = Path(__file__).resolve().parent
IMPORT = ROOT / "stage_c_import"
GEOMS = IMPORT / "geometry_only_inputs"
SOURCE_PACKAGE = ROOT.parent / "ti3sic2_DFT_extension_20261014"
SOURCE_GEOMETRY_VALIDATION = SOURCE_PACKAGE / "geometry_validation.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def formula(symbols: list[str]) -> str:
    counts = Counter(symbols)
    return "".join(f"{element}{counts[element]}" for element in ("Ag", "Ti", "Si", "C") if counts[element])


def main() -> None:
    # Verify the complete original Stage C source manifest, then check that
    # the imported manifest/validation and each geometry are exact copies.
    source_manifest = SOURCE_PACKAGE / "SHA256SUMS.txt"
    manifest_entries = []
    for line in source_manifest.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        expected, relative = line.split(maxsplit=1)
        relative = relative.lstrip("* ")
        source_path = SOURCE_PACKAGE / relative
        if not source_path.is_file() or sha256(source_path) != expected:
            raise SystemExit(f"Original Stage C manifest failed for {relative}")
        manifest_entries.append(relative)
    if len(manifest_entries) != 20:
        raise SystemExit(f"Expected 20 entries in original Stage C manifest, found {len(manifest_entries)}")

    with (IMPORT / "DFT_STRUCTURE_LIST.csv").open(newline="", encoding="utf-8") as stream:
        inventory = list(csv.DictReader(stream))
    source_validation = json.loads(SOURCE_GEOMETRY_VALIDATION.read_text())
    if sha256(IMPORT / "DFT_STRUCTURE_LIST.csv") != "0ee7bdfcdfec9a63c7b5b67e4c81b7aa05b0d7aa1bd720c60a19e46b22ac9fb2":
        raise SystemExit("Imported DFT_STRUCTURE_LIST.csv differs from the frozen source copy")
    if sha256(IMPORT / "geometry_validation.json") != "41dbdfcdf7479f437fee536027cfb42b5e87ee850b84416110ad7bdeae43e81d":
        raise SystemExit("Imported geometry_validation.json differs from the frozen source copy")
    if "START_DFT=NO" not in (SOURCE_PACKAGE / "FINAL_DECISION.md").read_text(encoding="utf-8"):
        raise SystemExit("Source Stage C package no longer declares START_DFT=NO")
    source_by_id = {row["id"]: row for row in source_validation["structures"]}

    audit_rows: list[dict[str, object]] = []
    for row in inventory:
        structure_id = row["structure_id"]
        status = row["structure_status"]
        if status != "PASS_GEOMETRY_ONLY":
            continue
        path = GEOMS / Path(row["structure_file"]).name
        actual_hash = sha256(path)
        expected_hash = row["structure_sha256"]
        atoms = read(path, format="extxyz")
        symbols = atoms.get_chemical_symbols()
        ids = atoms.arrays.get("local_new_id")
        coords = atoms.positions
        cell = atoms.cell.array
        expected_cell = np.asarray(json.loads(row["cell_A"]), dtype=float)
        distances = atoms.get_all_distances(mic=True)
        np.fill_diagonal(distances, np.inf)
        minimum = float(np.min(distances))
        source_record = source_by_id[structure_id]
        expected_formula = row["composition"]
        actual_formula = formula(symbols)
        pbc = "".join("T" if item else "F" for item in atoms.pbc)
        finite = bool(np.isfinite(coords).all() and np.isfinite(cell).all())
        ids_unique = bool(ids is not None and len(ids) == len(atoms) and len(set(map(int, ids))) == len(atoms))
        ids_one_based = bool(ids_unique and set(map(int, ids)) == set(range(1, len(atoms) + 1)))
        calc_results = getattr(getattr(atoms, "calc", None), "results", {}) or {}
        no_dft_labels = not (
            {"forces", "force"} & set(atoms.arrays)
            or {"energy", "free_energy", "stress"} & set(atoms.info)
            or {"forces", "force"} & set(calc_results)
            or {"energy", "free_energy", "stress"} & set(calc_results)
        )
        checks = {
            "sha_match": actual_hash == expected_hash == source_record["file_sha256"],
            "natoms_match": len(atoms) == int(row["natoms"]) == int(source_record["natoms"]),
            "composition_match": actual_formula == expected_formula,
            "pbc_match": pbc == row["pbc"] == "TTT",
            "cell_match": bool(np.allclose(cell, expected_cell, rtol=0.0, atol=1e-9)),
            "coordinates_and_cell_finite": finite,
            "atom_ids_unique_and_1based": ids_one_based,
            "minimum_distance_match": abs(minimum - float(row["minimum_distance_A"])) <= 2e-5,
            "geometry_validation_source_pass": bool(source_record["geometry_pass"]),
            "no_dft_labels_in_extxyz": no_dft_labels,
        }
        audit_rows.append({
            "structure_id": structure_id,
            "category": row["category"],
            "parent_family": row["parent_family"],
            "role_proposal": row["role_proposal"],
            "role_status": row["role_status"],
            "structure_sha256": actual_hash,
            "natoms": len(atoms),
            "composition": actual_formula,
            "pbc": pbc,
            "cell_matches_manifest": bool(np.allclose(cell, expected_cell, rtol=0.0, atol=1e-9)),
            "atom_id_field": "local_new_id",
            "atom_id_unique_1based": ids_one_based,
            "coordinates_cell_finite": finite,
            "minimum_distance_A": f"{minimum:.9f}",
            "dft_energy_force_stress_labels_present": not no_dft_labels,
            "audit_status": "PASS_GEOMETRY_ONLY" if all(checks.values()) else "FAIL",
            "checks_json": json.dumps(checks, sort_keys=True),
        })

    if len(audit_rows) != 10:
        raise SystemExit(f"Expected exactly 10 source-backed geometries, found {len(audit_rows)}")
    if any(row["audit_status"] != "PASS_GEOMETRY_ONLY" for row in audit_rows):
        raise SystemExit("At least one imported geometry failed an identity or geometry check")

    out = ROOT / "STRUCTURE_AUDIT.csv"
    with out.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(audit_rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(audit_rows)

    # Count label-bearing records only from the E5 source inventory. Stage C's
    # ten extxyz files are explicitly geometry-only; no labels are synthesized.
    result = {
        "scope": "offline data identity and geometry audit; no calculator or fitting invoked",
        "stage_c_source_commit": "0a8585e03522b88075e5f62770a8b5409031bf69",
        "stage_c_structure_list_sha256": sha256(IMPORT / "DFT_STRUCTURE_LIST.csv"),
        "stage_c_geometry_validation_sha256": sha256(IMPORT / "geometry_validation.json"),
        "geometry_count": len(audit_rows),
        "geometry_audit_pass_count": sum(r["audit_status"] == "PASS_GEOMETRY_ONLY" for r in audit_rows),
        "new_stage_c_dft_label_count": 0,
        "new_stage_c_energy_force_stress_records": 0,
        "source_families_in_stage_c": sorted({str(r["parent_family"]) for r in audit_rows}),
        "all_stage_c_geometry_candidates_share_one_parent_family": len({str(r["parent_family"]) for r in audit_rows}) == 1,
        "independent_dev_structures_present": False,
        "independent_test_structures_present_or_inspected": False,
        "sealed_test_data_accessed": False,
        "approved_frozen_train_count": 0,
        "approved_frozen_dev_count": 0,
        "approved_frozen_test_count": 0,
        "fit_ready_labelled_structure_count": 0,
        "existing_e5_dft_force_label_count_for_future_role_review": 6,
        "existing_e5_label_role_note": "Five 48-atom baseline-family bulk structures and one interface-derived 114-atom TSC-only slab were used in E5 diagnostics; retain current roles; no independent DEV/TEST use.",
        "stage_c_status": "GEOMETRIES_VALIDATED_LABELS_ABSENT",
    }
    (ROOT / "audit_summary.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"verified {len(audit_rows)} geometries; wrote STRUCTURE_AUDIT.csv and audit_summary.json")


if __name__ == "__main__":
    main()
