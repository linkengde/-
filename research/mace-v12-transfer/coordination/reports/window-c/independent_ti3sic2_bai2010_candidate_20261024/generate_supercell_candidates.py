#!/usr/bin/env python3
"""Create label-free supercell geometry proposals from the audited Bai2010 CIF."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from ase.build import make_supercell
from ase.io import read, write


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "Ti3SiC2_Bai2010_DFT_candidate.cif"
OUTDIR = ROOT / "supercell_candidates"
MANIFEST = ROOT / "supercell_geometry_manifest.json"
SOURCE_DOI = "10.1016/j.solidstatesciences.2010.03.007"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def geometry_sha(atoms) -> str:
    digest = hashlib.sha256()
    digest.update(" ".join(atoms.get_chemical_symbols()).encode())
    digest.update(np.asarray(atoms.cell).round(10).astype("<f8").tobytes())
    digest.update(np.asarray(atoms.positions).round(10).astype("<f8").tobytes())
    digest.update(np.asarray(atoms.pbc, dtype=np.uint8).tobytes())
    return digest.hexdigest()


def min_pair_distance(atoms) -> float:
    distances = atoms.get_all_distances(mic=True)
    np.fill_diagonal(distances, np.inf)
    return float(np.min(distances))


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    source_hash = sha256_bytes(SOURCE.read_bytes())
    source = read(SOURCE, format="cif")
    if source.get_chemical_formula() != "C4Si2Ti6" or len(source) != 12:
        raise ValueError("Unexpected source CIF expansion; refusing to generate candidates")

    families = {
        "Bai2010_TSC72_SX3_SY2_SZ1": np.diag([3, 2, 1]),
        "Bai2010_TSC96_SX2_SY2_SZ2": np.diag([2, 2, 2]),
    }
    rows = []
    for family, matrix in families.items():
        parent = make_supercell(source, matrix, wrap=True)
        if len(parent) not in (72, 96):
            raise ValueError(f"Unexpected atom count for {family}: {len(parent)}")
        variants = [("parent", None)] + [(f"displacement_seed_{seed}", seed) for seed in (3101, 9073)]
        for variant, seed in variants:
            atoms = parent.copy()
            if seed is not None:
                rng = np.random.default_rng(seed)
                displacements = rng.normal(0.0, 0.025, size=(len(atoms), 3))
                displacements -= displacements.mean(axis=0, keepdims=True)
                atoms.positions += displacements

            # Drop parser-specific CIF objects and attach only portable provenance.
            atoms.info.clear()
            atoms.info.update(
                {
                    "config_type": f"{family}_{variant}",
                    "candidate_role": "B1_TSC_Bai2010_geometry_candidate_pending_source_and_B_review",
                    "dataset_role": "geometry_only_no_DFT_labels_no_split_assignment",
                    "proposal_only_no_DFT_label": True,
                    "dft_label_present": False,
                    "dft_launch_enabled": False,
                    "source_paper": "Bai et al., Solid State Sciences 12 (2010), 1220-1225",
                    "source_doi": SOURCE_DOI,
                    "source_cif_sha256": source_hash,
                    "source_status": "independent_computational_source_candidate_exact_paper_parameters_pending_fulltext_check",
                    "source_scope": "independent_source_lineage_vs_COD_9009647_same_alpha_Ti3SiC2_structure_type",
                    "parent_family_id": family,
                    "supercell_matrix": json.dumps(matrix.tolist(), separators=(",", ":")),
                    "proposal_kind": "unperturbed_supercell" if seed is None else "small_random_displacement",
                    "proposal_seed": -1 if seed is None else seed,
                    "displacement_sigma_A": 0.0 if seed is None else 0.025,
                    "role_freeze_status": "DRAFT_NOT_FROZEN_B_REVIEW_REQUIRED",
                }
            )
            path = OUTDIR / f"{family}_{variant}.extxyz"
            write(path, atoms, format="extxyz")
            min_dist = min_pair_distance(atoms)
            file_hash = sha256_bytes(path.read_bytes())
            geom_hash = geometry_sha(atoms)
            rows.append(
                {
                    "file": path.name,
                    "file_sha256": file_hash,
                    "geometry_sha256": geom_hash,
                    "formula": atoms.get_chemical_formula(),
                    "natoms": len(atoms),
                    "minimum_pair_distance_A": min_dist,
                    "parent_family_id": family,
                    "source_doi": SOURCE_DOI,
                    "source_cif_sha256": source_hash,
                    "dft_label_present": False,
                    "dft_launch_enabled": False,
                    "proposed_role": "geometry-only; source and B review pending",
                }
            )

    unique_geometry = len({row["geometry_sha256"] for row in rows}) == len(rows)
    unique_files = len({row["file_sha256"] for row in rows}) == len(rows)
    min_distance_gate = all(row["minimum_pair_distance_A"] > 1.8 for row in rows)
    manifest = {
        "status": "PASS_GEOMETRY_ONLY" if unique_geometry and unique_files and min_distance_gate else "FAIL",
        "source_file": SOURCE.name,
        "source_cif_sha256": source_hash,
        "source_doi": SOURCE_DOI,
        "source_novelty_scope": "Independent computational-paper lineage relative to COD 9009647; same alpha-Ti3SiC2 structure type.",
        "records": rows,
        "integrity": {
            "record_count": len(rows),
            "unique_geometry_hashes": unique_geometry,
            "unique_file_hashes": unique_files,
            "all_minimum_pair_distances_above_1p8_A": min_distance_gate,
            "all_label_free_and_DFT_disabled": all(
                not row["dft_label_present"] and not row["dft_launch_enabled"] for row in rows
            ),
        },
        "limitations": [
            "Source paper parameter values remain pending independent full-text table verification.",
            "All structures are geometry proposals only; no relaxation, DFT labels, training, or validation role was assigned.",
            "The two supercell families are structures from the same source CIF and are not independent validation samples from separate publications.",
        ],
    }
    MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    if manifest["status"] != "PASS_GEOMETRY_ONLY":
        raise RuntimeError("Geometry integrity checks failed; inspect the manifest")
    print(
        json.dumps(
            {
                "status": manifest["status"],
                "records": len(rows),
                "atom_counts": sorted({row["natoms"] for row in rows}),
                "unique_geometry_hashes": unique_geometry,
                "min_distance_gate": min_distance_gate,
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
