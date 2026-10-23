#!/usr/bin/env python3
"""Read-only follow-up on V19 source and numerical-method gaps.

Reads approved training rows and exact DFT control archives only. It does not
read development, test, or sealed data and does not launch GPAW/MPI.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read

REPO = Path(__file__).resolve().parents[6]
BASE = REPO / "research/mace-v12-transfer/periodic_interface_v4"
OUT = Path(__file__).resolve().parent / "review.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def find_frame(frames, label):
    matches = [(i, atoms) for i, atoms in enumerate(frames) if str(atoms.info.get("config_type", "")).startswith(label)]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one frame for {label}, got {len(matches)}")
    return matches[0]


def main() -> None:
    v16_train_path = BASE / "mace_periodic_v16_interface_energy/data/train.extxyz"
    v19_train_path = BASE / "mace_periodic_v19_residual_response/data/train.extxyz"
    v16_old = read(v16_train_path, index=0)
    v19_train = read(v19_train_path, index=":")
    v19_index, v19_row = find_frame(v19_train, "Ti_gap_2p3_rank_recovery_v17_01")

    rank_root = BASE / "pbe_interface_v17_rank_recovery_parallel"
    label = "Ti_gap_2p3_rank_recovery_v17_01"
    label_dir = rank_root / "calculations" / label
    output_path = label_dir / f"{label}_PW_PBE.extxyz"
    summary_path = label_dir / "summary.json"
    verification_path = label_dir / "verification.json"
    output = read(output_path)
    summary = json.loads(summary_path.read_text())
    verification = json.loads(verification_path.read_text())
    force_delta = v16_old.arrays["REF_forces"] - output.arrays["PW_PBE_forces"]
    geometry_same = {
        "ordered_numbers": bool(np.array_equal(v16_old.numbers, output.numbers)),
        "positions_exact": bool(np.array_equal(v16_old.positions, output.positions)),
        "cell_exact": bool(np.array_equal(v16_old.cell.array, output.cell.array)),
        "pbc_equal": bool(np.array_equal(v16_old.pbc, output.pbc)),
    }
    current_matches_archive = {
        "ordered_numbers": bool(np.array_equal(v19_row.numbers, output.numbers)),
        "positions_exact": bool(np.array_equal(v19_row.positions, output.positions)),
        "cell_exact": bool(np.array_equal(v19_row.cell.array, output.cell.array)),
        "forces_exact": bool(np.array_equal(v19_row.arrays["REF_forces"], output.arrays["PW_PBE_forces"])),
        "energy_exact": float(v19_row.info["REF_energy"]) == float(output.info["PW_PBE_energy_eV"]),
    }

    dipole_report_path = REPO / "research/mace-v12-transfer/coordination/reports/window-a/v19_boundary_dipole_comparison/report.json"
    dipole_report = json.loads(dipole_report_path.read_text())
    controls = dipole_report["matched_parent_dipole_and_vacuum_controls"]

    result = {
        "status": "READ_ONLY_METHOD_AND_SOURCE_FOLLOWUP",
        "scope": "Approved training data, exact source archives, and existing fixed-geometry control report only. No valid/dev/test/sealed file read; no DFT/MPI run.",
        "legacy_rank_recovery_row": {
            "historical_v16_train_frame": 0,
            "current_v19_train_frame": v19_index,
            "historical_method_string": v16_old.info.get("dft_method"),
            "current_verified_method": summary["method"],
            "legacy_geometry_matches_current_output": geometry_same,
            "current_v19_row_matches_verified_gamma_archive": current_matches_archive,
            "historical_energy_eV_cell": float(v16_old.info["REF_energy"]),
            "verified_gamma_energy_eV_cell": float(output.info["PW_PBE_energy_eV"]),
            "gamma_minus_historical_energy_eV_cell": float(output.info["PW_PBE_energy_eV"] - v16_old.info["REF_energy"]),
            "energy_difference_meV_atom": float((output.info["PW_PBE_energy_eV"] - v16_old.info["REF_energy"]) * 1000 / len(output)),
            "force_component_RMSE_eV_A": float(np.sqrt(np.mean(force_delta**2))),
            "force_vector_RMSE_eV_A": float(np.sqrt(np.mean(np.sum(force_delta**2, axis=1)))),
            "maximum_atom_vector_difference_eV_A": float(np.linalg.norm(force_delta, axis=1).max()),
            "verified_gamma_archive_status": verification.get("status"),
            "verified_gamma_output_sha256": sha(output_path),
            "verified_gamma_summary_sha256": sha(summary_path),
            "interpretation": "V19 uses the separately archived, verified Gamma label. It exactly matches the current V19 training row. The older V16 value has the same geometry but no original SCF log, verified output identity, or getter evidence; the large difference cannot be assigned to k-point sampling alone and the legacy value must not be substituted into train36.",
        },
        "fixed_ag_si_parent_boundary_controls": {
            "source_report": dipole_report_path.relative_to(REPO).as_posix(),
            "source_report_sha256": sha(dipole_report_path),
            "status": controls["status"],
            "dataset_role": controls["dataset_role"],
            "comparisons": controls["pairs"],
            "interpretation": "The 26 A dipole ON/OFF and 26 A to 36 A dipole ON comparisons pass the existing fixed-geometry energy and force-vector budgets. This does not establish nonperiodic-z slab k-mesh or smearing convergence.",
        },
        "independent_parent_search": {
            "qualified_second_Ti3SiC2_parent_found": False,
            "local_checkout": "Only CIF copies identifying COD 9009647 were found; existing local audit normalized both copies to the same CIF bytes and publication/source family.",
            "hidden_workspace_cif_inventory": {
                "scope": "/workspace including hidden directories and installed environments",
                "unique_ASE_mineral_fixture_CIFs": 8,
                "environment_copies": 3,
                "additional_project_CIFs": 2,
                "other_CIF_or_MCIF_files_found": 0,
                "reduced_Ti3SiC2_matches": 0,
                "note": "The eight ASE fixture CIFs were parsed; their formulas are Ni2Sb2, CSi, H4Al2F6O4Pb2, CoAsS, Co2As2S2, Ni3S2, MoS3, and SiO2. The two project CIFs both identify COD 9009647."
            },
            "external_requests": {
                "materials_project_optimade": {"url": "https://optimade.materialsproject.org/v1/structures?filter=chemical_formula_reduced%3D%22Ti3SiC2%22", "result": "Tunnel connection failed: HTTP 403 Forbidden"},
                "oqmd_optimade": {"url": "https://oqmd.org/optimade/v1/structures?filter=chemical_formula_reduced%3D%22Ti3SiC2%22", "result": "Tunnel connection failed: HTTP 403 Forbidden"},
                "nomad_optimade": {"url": "https://nomad-lab.eu/prod/v1/api/v1/optimade/structures?filter=chemical_formula_reduced%3D%22Ti3SiC2%22", "result": "Tunnel connection failed: HTTP 403 Forbidden"},
                "cod_direct_cif": {"url": "https://www.crystallography.net/cod/cif/9009647.cif", "result": "Tunnel connection failed: HTTP 403 Forbidden"},
                "cod_ugr_mirror": {"url": "https://qiserver.ugr.es/cod/cif/9009647.cif", "result": "Tunnel connection failed: HTTP 403 Forbidden"},
            },
            "public_git_fixture_scans": [
                {"repository": "https://github.com/materialsproject/pymatgen.git", "commit": "c6507866179a6c37ee7d0df6cbb79896a2ebba52", "cif_mcif_files_scanned": 88, "Ti3SiC2_reduced_formula_matches": 0},
                {"repository": "https://github.com/materialsproject/pymatgen-test-files.git", "commit": "85005108e99e40cbd59e78f067459b0ca664efc6", "cif_mcif_files_scanned": 110, "Ti3SiC2_reduced_formula_matches": 0, "json_structure_matches": 0},
                {"repository": "https://github.com/cod-developers/cod-tools.git", "commit": "0060d4652a5a262fbda2276801f1244060d67128", "cif_mcif_files_scanned": 1516, "chemical_formula_fields_scanned": 2065, "Ti3SiC2_reduced_formula_matches": 0},
            ],
            "decision": "Do not construct or label an independent development parent from a perturbation of COD 9009647 or from an untraceable software fixture. Source access remains the hard blocker.",
        },
        "remaining_numerical_coverage": {
            "validated_for_this_frozen_parent": ["26 A dipole ON/OFF", "26 A to 36 A with dipole ON", "6x6x1, sigma 0.10 eV at the registered 26 A parent"],
            "not_yet_validated_for_nonperiodic_z_slab": ["matched 4x4x1/6x6x1/8x8x1 k-mesh ladder", "matched sigma 0.05/0.10/0.20 eV ladder"],
            "controls_are_not_training_or_independent_validation": True,
        },
    }
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": result["status"], "current_row_archive_match": current_matches_archive, "second_parent_found": False, "output": OUT.relative_to(REPO).as_posix()}, indent=2))


if __name__ == "__main__":
    main()
