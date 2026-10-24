#!/usr/bin/env python3
"""Independently verify the serialized Bai2010 label-free supercell proposals."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "Ti3SiC2_Bai2010_DFT_candidate.cif"
OUTDIR = ROOT / "supercell_candidates"
OUTPUT = ROOT / "supercell_integrity_audit.json"
SOURCE_DOI = "10.1016/j.solidstatesciences.2010.03.007"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def geometry_sha(atoms) -> str:
    digest = hashlib.sha256()
    digest.update(" ".join(atoms.get_chemical_symbols()).encode())
    digest.update(np.asarray(atoms.cell).round(10).astype("<f8").tobytes())
    digest.update(np.asarray(atoms.positions).round(10).astype("<f8").tobytes())
    digest.update(np.asarray(atoms.pbc, dtype=np.uint8).tobytes())
    return digest.hexdigest()


def main() -> None:
    source_hash = sha256(SOURCE)
    expected = {
        "Bai2010_TSC72_SX3_SY2_SZ1": (72, "C24Si12Ti36"),
        "Bai2010_TSC96_SX2_SY2_SZ2": (96, "C32Si16Ti48"),
    }
    paths = sorted(OUTDIR.glob("*.extxyz"))
    rows = []
    for path in paths:
        atoms = read(path, format="extxyz")
        family = atoms.info.get("parent_family_id", "")
        expected_n, expected_formula = expected.get(family, (-1, ""))
        distances = atoms.get_all_distances(mic=True)
        np.fill_diagonal(distances, np.inf)
        minimum = float(np.min(distances))
        finite = bool(np.isfinite(atoms.positions).all() and np.isfinite(np.asarray(atoms.cell)).all())
        checks = {
            "expected_family": family in expected,
            "expected_atom_count": len(atoms) == expected_n,
            "expected_formula": atoms.get_chemical_formula() == expected_formula,
            "fully_periodic": bool(atoms.pbc.all()),
            "finite_geometry": finite,
            "minimum_pair_distance_above_1p8_A": minimum > 1.8,
            "correct_source_cif_hash": atoms.info.get("source_cif_sha256") == source_hash,
            "correct_source_doi": atoms.info.get("source_doi") == SOURCE_DOI,
            "no_dft_label": atoms.info.get("proposal_only_no_DFT_label") is True
            and atoms.info.get("dft_label_present") is False,
            "dft_launch_disabled": atoms.info.get("dft_launch_enabled") is False,
            "role_unfrozen": atoms.info.get("role_freeze_status") == "DRAFT_NOT_FROZEN_B_REVIEW_REQUIRED",
        }
        rows.append(
            {
                "file": path.name,
                "file_sha256": sha256(path),
                "geometry_sha256": geometry_sha(atoms),
                "formula": atoms.get_chemical_formula(),
                "natoms": len(atoms),
                "minimum_pair_distance_A": minimum,
                "checks": checks,
                "status": "PASS" if all(checks.values()) else "FAIL",
            }
        )

    unique_files = len({row["file_sha256"] for row in rows}) == len(rows)
    unique_geometries = len({row["geometry_sha256"] for row in rows}) == len(rows)
    payload = {
        "status": "PASS_GEOMETRY_ONLY" if len(rows) == 6 and unique_files and unique_geometries and all(r["status"] == "PASS" for r in rows) else "FAIL",
        "source_cif_sha256": source_hash,
        "record_count": len(rows),
        "unique_file_hashes": unique_files,
        "unique_geometry_hashes": unique_geometries,
        "records": rows,
        "limitations": [
            "The source paper's exact numerical parameters remain pending full-text verification.",
            "These are unlabeled geometry proposals; no DFT, training, validation, or production MD was run.",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": payload["status"], "records": len(rows), "unique_files": unique_files, "unique_geometries": unique_geometries}, ensure_ascii=False))
    if payload["status"] != "PASS_GEOMETRY_ONLY":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
