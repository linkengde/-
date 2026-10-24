#!/usr/bin/env python3
"""Build reproducible, label-free fcc Ag candidates from COD 9013045 lattice data."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from ase.build import bulk
from ase.io import read, write


HERE = Path(__file__).resolve().parent
OUT = HERE / "candidates"
LATTICE_A = 4.086
SOURCE_COD = "9013045"
SOURCE_CIF_SHA256 = "6dbb7ae83a23963204724bb2da599fbe6304a6c875d5cc3f28429639a582fe31"
SOURCE_DOI = "10.1007/BF01174717"
DISPLACEMENT_SIGMA_A = 0.025
DISPLACEMENT_SEEDS = (3101, 9073)


def decorate(atoms, *, config: str, family: str, variant: str, seed: int | None) -> None:
    atoms.info.update({
        "config_type": config,
        "candidate_role": "B1_Ag_training_candidate_pending_role_audit",
        "dataset_role": "geometry_only_proposal_no_DFT_label",
        "proposal_only_no_DFT_label": True,
        "dft_launch_enabled": False,
        "source_database": "Crystallography Open Database",
        "source_accession": f"COD {SOURCE_COD}",
        "source_cif_url": f"https://www.crystallography.net/cod/{SOURCE_COD}.cif",
        "source_cif_sha256": SOURCE_CIF_SHA256,
        "source_doi": SOURCE_DOI,
        "source_citation": "I.-K. Suh, H. Ohta, and Y. Waseda, J. Mater. Sci. 23 (1988) 757-760",
        "source_access_note": "COD record identifies the cited paper and experimental 293 K structure; cite the paper when reusing the record.",
        "source_lattice_parameter_A": LATTICE_A,
        "source_measurement_temperature_K": 293,
        "source_crystal_system": "cubic",
        "source_space_group": "Fm-3m (No. 225)",
        "source_family_id": family,
        "parent_family_id": family,
        "candidate_variant": variant,
        "displacement_sigma_A": DISPLACEMENT_SIGMA_A if seed is not None else 0.0,
        "displacement_seed": int(seed) if seed is not None else -1,
        "generator": "ASE bulk('Ag', 'fcc', a=4.086 Å); cubic conventional or primitive cell; seeded Gaussian displacement variants",
        "training_label_policy": "B source/role review and separate batch authorization required",
    })


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    records = []
    for family_name, cubic, repeat in (
        ("Ag32_COD9013045_cubic_2x2x2", True, (2, 2, 2)),
        ("Ag64_COD9013045_primitive_4x4x4", False, (4, 4, 4)),
    ):
        parent = bulk("Ag", "fcc", a=LATTICE_A, cubic=cubic).repeat(repeat)
        variants: list[tuple[str, int | None, object]] = [("parent", None, parent.copy())]
        for seed in DISPLACEMENT_SEEDS:
            variant = parent.copy()
            rng = np.random.default_rng(seed)
            variant.positions += rng.normal(0.0, DISPLACEMENT_SIGMA_A, variant.positions.shape)
            variant.wrap()
            variants.append((f"displacement_seed_{seed}", seed, variant))

        for variant_name, seed, atoms in variants:
            config = f"{family_name}_{variant_name}"
            decorate(atoms, config=config, family=family_name, variant=variant_name, seed=seed)
            distances = atoms.get_all_distances(mic=True)
            np.fill_diagonal(distances, np.inf)
            min_dist = float(np.min(distances))
            if (
                atoms.get_chemical_formula() not in {"Ag32", "Ag64"}
                or not np.isfinite(atoms.positions).all()
                or not np.isfinite(np.asarray(atoms.cell)).all()
                or not atoms.pbc.all()
                or min_dist < 2.5
            ):
                raise ValueError(f"Invalid candidate geometry: {config}")

            path = OUT / f"{config}.extxyz"
            write(path, atoms, format="extxyz")
            loaded = read(path)
            if (
                loaded.get_chemical_formula() != atoms.get_chemical_formula()
                or loaded.info["source_cif_sha256"] != SOURCE_CIF_SHA256
                or loaded.info["proposal_only_no_DFT_label"] is not True
                or loaded.info["dft_launch_enabled"] is not False
            ):
                raise ValueError(f"Serialized candidate failed read-back checks: {path}")

            records.append({
                "config_type": config,
                "file": str(path.relative_to(HERE)),
                "file_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "formula": loaded.get_chemical_formula(),
                "atoms": len(loaded),
                "source_family_id": family_name,
                "variant": variant_name,
                "seed": seed,
                "cell_lengths_A": [float(x) for x in loaded.cell.lengths()],
                "volume_A3": float(loaded.get_volume()),
                "minimum_pair_distance_A": min_dist,
                "no_DFT_label": loaded.info["proposal_only_no_DFT_label"],
                "DFT_launch_enabled": loaded.info["dft_launch_enabled"],
            })

    if len(records) != 6 or len({r["file_sha256"] for r in records}) != 6:
        raise ValueError("Expected six unique solid-Ag geometry proposals")
    manifest = {
        "status": "PASS_GEOMETRY_ONLY_SOURCE_LINKED_CANDIDATES",
        "candidate_count": len(records),
        "source": {
            "cod_accession": SOURCE_COD,
            "cif_url": f"https://www.crystallography.net/cod/{SOURCE_COD}.cif",
            "downloaded_cif_bytes": 4769,
            "downloaded_cif_sha256": SOURCE_CIF_SHA256,
            "citation": "I.-K. Suh, H. Ohta, and Y. Waseda, High-temperature thermal expansion of six metallic elements measured by dilatation method and X-ray diffraction, Journal of Materials Science 23 (1988) 757–760.",
            "doi": SOURCE_DOI,
            "temperature_K": 293,
            "lattice_parameter_A": LATTICE_A,
            "space_group": "Fm-3m (No. 225)",
            "formula_units_per_conventional_cell": 4,
            "source_note": "COD 9013045 record links the 293 K structure and article citation; the publisher full text was not retrieved in this audit.",
        },
        "geometry_method": "ASE-built fcc cells at the COD lattice parameter; two seeded 0.025 Å displacement variants per supercell family.",
        "role_limit": "One source family for training-candidate review only; no dev/test role, no labels, no dataset integration, no DFT authorization.",
        "records": records,
    }
    (HERE / "candidate_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({
        "status": manifest["status"],
        "candidate_count": len(records),
        "minimum_pair_distance_A": min(r["minimum_pair_distance_A"] for r in records),
        "unique_hashes": len({r["file_sha256"] for r in records}),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
