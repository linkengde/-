#!/usr/bin/env python3
"""Audit the supplied reconstructed Ti3SiC2 CIF without rewriting it."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import numpy as np
from ase.io import read
from ase.neighborlist import neighbor_list


ROOT = Path(__file__).resolve().parent
CIF = ROOT / "Ti3SiC2_Bai2010_DFT_candidate.cif"
SOURCE_AUDIT = ROOT / "source_audit_user_provided.md"
OUTPUT = ROOT / "geometry_audit.json"
EXPECTED_HASH = "cb361a88ab35948e9fc8e0e8f86b0ec7030ae0934c3bbf7c33ef68b674d3d3a5"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    atoms = read(CIF, format="cif")
    scaled = atoms.get_scaled_positions(wrap=True)
    ti_2a = [
        i
        for i, atom in enumerate(atoms)
        if atom.symbol == "Ti" and np.linalg.norm(scaled[i, :2]) < 1e-7
    ]
    ti_4f = [i for i, atom in enumerate(atoms) if atom.symbol == "Ti" and i not in ti_2a]
    carbon = [i for i, atom in enumerate(atoms) if atom.symbol == "C"]
    silicon = [i for i, atom in enumerate(atoms) if atom.symbol == "Si"]
    all_ti = [i for i, atom in enumerate(atoms) if atom.symbol == "Ti"]
    nearest_2a_c = min(atoms.get_distance(i, j, mic=True) for i in ti_2a for j in carbon)
    nearest_4f_c = min(atoms.get_distance(i, j, mic=True) for i in ti_4f for j in carbon)
    nearest_ti_si = min(atoms.get_distance(i, j, mic=True) for i in all_ti for j in silicon)
    neighbor_i, neighbor_j, neighbor_d = neighbor_list("ijd", atoms, cutoff=0.5)
    overlap_pairs = [
        {"i": int(i), "j": int(j), "distance_angstrom": float(d)}
        for i, j, d in zip(neighbor_i, neighbor_j, neighbor_d)
        if i < j
    ]

    cif_hash = sha256(CIF)
    claimed_hashes = re.findall(r"[0-9a-f]{64}", SOURCE_AUDIT.read_text(encoding="utf-8"))
    sg = atoms.info.get("spacegroup")
    frac = atoms.get_scaled_positions(wrap=True)
    ti_4f_z = [
        float(frac[i, 2]) for i, atom in enumerate(atoms)
        if atom.symbol == "Ti" and i not in ti_2a and 0.0 < frac[i, 2] < 0.5
    ]
    carbon_z = [float(frac[i, 2]) for i in carbon if frac[i, 2] > 0.5]
    source_parameters = {
        "a_angstrom": float(atoms.cell.lengths()[0]),
        "c_angstrom": float(atoms.cell.lengths()[2]),
        "zTi2": min(ti_4f_z) if ti_4f_z else None,
        "zC": min(carbon_z) if carbon_z else None,
        "matches_Bai2010_Table1": bool(
            np.isclose(atoms.cell.lengths()[0], 3.0592, atol=5e-5)
            and np.isclose(atoms.cell.lengths()[2], 17.6360, atol=5e-5)
            and bool(ti_4f_z)
            and np.isclose(min(ti_4f_z), 0.1351, atol=5e-5)
            and bool(carbon_z)
            and np.isclose(min(carbon_z), 0.5726, atol=5e-5)
        ),
    }
    checks = {
        "cif_sha256": cif_hash,
        "matches_hash_in_user_audit": cif_hash in claimed_hashes,
        "matches_expected_hash": cif_hash == EXPECTED_HASH,
        "atom_count": len(atoms),
        "formula": atoms.get_chemical_formula(),
        "formula_units_per_cell": 2 if len(atoms) == 12 else None,
        "site_counts_after_symmetry_expansion": {
            "Ti_2a": len(ti_2a),
            "Ti_4f": len(ti_4f),
            "Si_2b": len(silicon),
            "C_4f": len(carbon),
        },
        "space_group_number_as_parsed_by_ASE": getattr(sg, "no", None),
        "cell_lengths_angstrom": [float(x) for x in atoms.cell.lengths()],
        "cell_angles_degrees": [float(x) for x in atoms.cell.angles()],
        "volume_angstrom3": float(atoms.get_volume()),
        "minimum_Ti2a_C_distance_angstrom": float(nearest_2a_c),
        "minimum_Ti4f_C_distance_angstrom": float(nearest_4f_c),
        "minimum_Ti_Si_distance_angstrom": float(nearest_ti_si),
        "overlap_cutoff_angstrom": 0.5,
        "overlap_pairs": overlap_pairs,
        "all_expanded_site_occupancies_full": all(
            abs(float(occ) - 1.0) < 1e-8
            for site in atoms.info.get("occupancy", {}).values()
            for occ in site.values()
        ),
        "Bai2010_Table1_parameters": source_parameters,
    }

    geometry_checks_passed = (
        checks["cif_sha256"] == EXPECTED_HASH
        and checks["matches_hash_in_user_audit"]
        and checks["atom_count"] == 12
        and checks["space_group_number_as_parsed_by_ASE"] == 194
        and checks["site_counts_after_symmetry_expansion"] == {
            "Ti_2a": 2,
            "Ti_4f": 4,
            "Si_2b": 2,
            "C_4f": 4,
        }
        and checks["all_expanded_site_occupancies_full"]
        and source_parameters["matches_Bai2010_Table1"]
        and not checks["overlap_pairs"]
    )
    result = {
        "status": "geometry_only_candidate_pass" if geometry_checks_passed else "geometry_audit_failed",
        "geometry_checks_passed": geometry_checks_passed,
        "source": {
            "paper": "Y. Bai et al., Chemical bonding and elastic properties of Ti3AC2 phases (A = Si, Ge, and Sn): A first-principle study, Solid State Sciences 12 (2010), 1220-1225.",
            "doi": "10.1016/j.solidstatesciences.2010.03.007",
            "bibliographic_metadata_check": "Crossref and OpenAlex metadata agree on title, authorship, journal, and 2010 publication year.",
            "parameter_evidence": "User-supplied accepted manuscript, Table 1 (printed page 18; PDF page 20), visually checked. SHA256 is recorded in ../provided_primary_papers_audit_20261024/README.md.",
            "parameter_check": "Candidate CIF parameters are checked against the values transcribed from the supplied Table 1 in the checks object.",
            "parameter_check_limit": "The accepted manuscript PDF is user-supplied; no publisher-hosted source-original CIF was obtained.",
            "independence_scope": "Independent computational paper relative to COD 9009647; same alpha-Ti3SiC2 structure type, not a distinct experimental record or polymorph.",
            "cif_origin": "Reconstructed derivative supplied by the user; not a downloaded source-original CIF.",
        },
        "checks": checks,
        "limitations": [
            "No original source CIF was supplied or downloaded; this check verifies the supplied reconstructed derivative only.",
            "The CIF is a reconstructed derivative; source-paper values are checked, but no publisher-hosted source-original CIF was obtained.",
            "No DFT calculation, relaxation, model training, or independent validation was performed.",
        ],
    }
    OUTPUT.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(checks, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
