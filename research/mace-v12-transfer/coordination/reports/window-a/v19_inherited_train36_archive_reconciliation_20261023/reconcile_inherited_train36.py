#!/usr/bin/env python3
"""Reconcile the 12 formerly unresolved V19 train36 rows to archived DFT outputs.

Reads only the approved V19 train.extxyz plus the exact historical input/output
archives named below. It does not read validation/test/sealed dataset files or
run GPAW/MPI.
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import numpy as np
from ase.io import read

REPO = Path(__file__).resolve().parents[6]
BASE = REPO / "research/mace-v12-transfer/periodic_interface_v4"
TRAIN = BASE / "mace_periodic_v19_residual_response/data/train.extxyz"
OUT = Path(__file__).resolve().parent / "reconciliation.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def record(manifest: Path, label: str) -> tuple[dict, str]:
    data = json.loads(manifest.read_text())
    rows = [row for row in data["records"] if row["label"] == label]
    if len(rows) != 1:
        raise RuntimeError(f"Expected one manifest record for {label} in {manifest}")
    return rows[0], data.get("owner_task", "historical archive")


SOURCES = [
    ("AgSi_d2p2", "pbe_interface_energy_additions_v12/calculations/AgSi_d2p2", "pbe_interface_energy_additions_v12/inputs/AgSi_d2p2.extxyz", "pbe_interface_energy_additions_v12/archive_manifest.json", "v12_contact"),
    ("AgSi_d2p8", "pbe_interface_energy_additions_v12/calculations/AgSi_d2p8", "pbe_interface_energy_additions_v12/inputs/AgSi_d2p8.extxyz", "pbe_interface_energy_additions_v12/archive_manifest.json", "v12_contact"),
    ("AgC_d2p0", "pbe_interface_energy_additions_v12/calculations/AgC_d2p0", "pbe_interface_energy_additions_v12/inputs/AgC_d2p0.extxyz", "pbe_interface_energy_additions_v12/archive_manifest.json", "v12_contact"),
    ("AgC_d2p4", "pbe_interface_energy_additions_v12/calculations/AgC_d2p4", "pbe_interface_energy_additions_v12/inputs/AgC_d2p4.extxyz", "pbe_interface_energy_additions_v12/archive_manifest.json", "v12_contact"),
    ("AgC_validation_d2p30", "pbe_interface_energy_additions_v12/calculations/validation_distance_scans/AgC_validation_d2p30", "pbe_interface_energy_additions_v12/inputs/validation_distance_scans/AgC_validation_d2p30.extxyz", "pbe_interface_energy_additions_v12/calculations/validation_distance_scans/archive_manifest.json", "v12_distance_scan"),
    ("AgSi_validation_d2p60", "pbe_interface_energy_additions_v12/calculations/validation_distance_scans/AgSi_validation_d2p60", "pbe_interface_energy_additions_v12/inputs/validation_distance_scans/AgSi_validation_d2p60.extxyz", "pbe_interface_energy_additions_v12/calculations/validation_distance_scans/archive_manifest.json", "v12_distance_scan"),
    ("AgTi_validation_d2p60", "pbe_interface_energy_additions_v12/calculations/validation_distance_scans/AgTi_validation_d2p60", "pbe_interface_energy_additions_v12/inputs/validation_distance_scans/AgTi_validation_d2p60.extxyz", "pbe_interface_energy_additions_v12/calculations/validation_distance_scans/archive_manifest.json", "v12_distance_scan"),
    ("AgC_lateral_registry_holdout_v13", "pbe_interface_v13_holdouts/calculations/AgC_lateral_registry_holdout_v13", "pbe_interface_v13_holdouts/inputs/AgC_lateral_registry_holdout_v13.extxyz", "pbe_interface_v13_holdouts/archive_manifest.json", "v13_lateral_holdout_now_train"),
    ("AgSi_lateral_registry_holdout_v13", "pbe_interface_v13_holdouts/calculations/AgSi_lateral_registry_holdout_v13", "pbe_interface_v13_holdouts/inputs/AgSi_lateral_registry_holdout_v13.extxyz", "pbe_interface_v13_holdouts/archive_manifest.json", "v13_lateral_holdout_now_train"),
    ("AgTi_registry_probe_v13_01", "pbe_interface_v14_parallel_acquisition/calculations/AgTi_registry_probe_v13_01", "pbe_interface_v14_parallel_acquisition/inputs/AgTi_registry_probe_v13_01.extxyz", "pbe_interface_v14_parallel_acquisition/archive_manifest.json", "v14_acquisition"),
    ("AgC_registry_strain_acq_v14_01", "pbe_interface_v14_parallel_acquisition/calculations/AgC_registry_strain_acq_v14_01", "pbe_interface_v14_parallel_acquisition/inputs/AgC_registry_strain_acq_v14_01.extxyz", "pbe_interface_v14_parallel_acquisition/archive_manifest.json", "v14_acquisition"),
    ("AgSi_registry_strain_acq_v14_01", "pbe_interface_v14_parallel_acquisition/calculations/AgSi_registry_strain_acq_v14_01", "pbe_interface_v14_parallel_acquisition/inputs/AgSi_registry_strain_acq_v14_01.extxyz", "pbe_interface_v14_parallel_acquisition/archive_manifest.json", "v14_acquisition"),
]


def main() -> None:
    train = read(TRAIN, index=":")
    train_hash = sha(TRAIN)
    rows = []
    for label, folder_rel, input_rel, manifest_rel, source_class in SOURCES:
        folder = BASE / folder_rel
        input_path = BASE / input_rel
        manifest_path = BASE / manifest_rel
        manifest_row, _ = record(manifest_path, label)
        input_atoms = read(input_path)
        output_path = folder / f"{label}_PW_PBE.extxyz"
        output_atoms = read(output_path)
        summary_path = folder / "summary.json"
        summary = json.loads(summary_path.read_text())
        progress_path = folder / "progress.json"
        progress = json.loads(progress_path.read_text())
        log_path = folder / "gpaw.log"
        log = log_path.read_text(errors="replace")
        backend_match = re.search(r"GPAW_MPI_BACKEND:\s*(\S+)", log)
        domain_match = re.search(r"^\s*domain:\s*(\d+)\s*$", log, re.M)
        matches = [(i, atoms) for i, atoms in enumerate(train) if str(atoms.info.get("config_type", "")).startswith(label + "_")]
        if len(matches) != 1:
            raise RuntimeError(f"Expected unique train36 row for {label}, got {len(matches)}")
        index, frame = matches[0]

        # Older compact manifests use source_sha256 for the input hash. Newer
        # manifests use input_sha256 and hold per-file output hashes in sha256.
        expected_input = manifest_row.get("input_sha256", manifest_row.get("source_sha256"))
        expected_output = manifest_row.get("extxyz_sha256", manifest_row.get("sha256", {}).get(output_path.name))
        expected_summary = manifest_row.get("summary_sha256", manifest_row.get("sha256", {}).get("summary.json"))
        iteration = int(summary["scf_iterations"])
        method = summary.get("method", {})
        geometry_checks = {
            "input_output_atomic_numbers": bool(np.array_equal(input_atoms.numbers, output_atoms.numbers)),
            "train_output_atomic_numbers": bool(np.array_equal(frame.numbers, output_atoms.numbers)),
            "input_output_lammps_ids": bool(np.array_equal(input_atoms.arrays["lammps_id"], output_atoms.arrays["lammps_id"])),
            "train_output_lammps_ids": bool(np.array_equal(frame.arrays["lammps_id"], output_atoms.arrays["lammps_id"])),
            "input_output_positions_exact": bool(np.array_equal(input_atoms.positions, output_atoms.positions)),
            "train_output_positions_exact": bool(np.array_equal(frame.positions, output_atoms.positions)),
            "input_output_cell_exact": bool(np.array_equal(input_atoms.cell.array, output_atoms.cell.array)),
            "train_output_cell_exact": bool(np.array_equal(frame.cell.array, output_atoms.cell.array)),
            "train_output_pbc_equal": bool(np.array_equal(frame.pbc, output_atoms.pbc)),
            "input_output_pbc_equal": bool(np.array_equal(input_atoms.pbc, output_atoms.pbc)),
            "train_output_force_max_abs_eV_A": float(np.max(np.abs(frame.arrays["REF_forces"] - output_atoms.arrays["PW_PBE_forces"]))),
            "train_output_force_exact": bool(np.array_equal(frame.arrays["REF_forces"], output_atoms.arrays["PW_PBE_forces"])),
        }
        energies = {
            "train_REF_energy_eV_cell": float(frame.info["REF_energy"]),
            "summary_energy_eV_cell": float(summary["energy_eV_cell"]),
            "output_PW_PBE_energy_eV_cell": float(output_atoms.info["PW_PBE_energy_eV"]),
            "train_summary_abs_difference_eV": abs(float(frame.info["REF_energy"]) - float(summary["energy_eV_cell"])),
            "output_summary_abs_difference_eV": abs(float(output_atoms.info["PW_PBE_energy_eV"]) - float(summary["energy_eV_cell"])),
            "free_energy_stored_in_output": any(key in output_atoms.info for key in ("free_energy", "free_energy_eV", "PW_PBE_free_energy_eV")),
        }
        hashes = {
            "input_sha256": sha(input_path),
            "manifest_input_sha256": expected_input,
            "summary_source_sha256": summary.get("source_sha256"),
            "output_sha256": sha(output_path),
            "manifest_output_sha256": expected_output,
            "summary_sha256": sha(summary_path),
            "manifest_summary_sha256": expected_summary,
            "input_match": sha(input_path) == expected_input == summary.get("source_sha256"),
            "output_match": expected_output is not None and sha(output_path) == expected_output,
            "summary_match": expected_summary is not None and sha(summary_path) == expected_summary,
        }
        checks = {
            "hashes": all(hashes[key] for key in ("input_match", "output_match", "summary_match")),
            "geometry": all(geometry_checks[key] for key in (
                "input_output_atomic_numbers", "train_output_atomic_numbers", "input_output_lammps_ids",
                "train_output_lammps_ids", "input_output_positions_exact", "train_output_positions_exact",
                "input_output_cell_exact",
                "train_output_cell_exact", "train_output_pbc_equal", "input_output_pbc_equal", "train_output_force_exact")),
            "energy": energies["train_summary_abs_difference_eV"] <= 1e-10 and energies["output_summary_abs_difference_eV"] <= 1e-10,
            "scf": summary.get("scf_converged") is True and progress.get("status") == "complete" and bool(re.search(rf"Converged in\s+{iteration}\s+(?:steps|iterations)", log, re.I)),
            "summary_geometry": summary.get("atoms") == len(output_atoms) and summary.get("formula") == output_atoms.get_chemical_formula(),
            "method": method.get("xc") == "PBE" and method.get("basis") == "plane wave" and method.get("cutoff_eV") == 500 and method.get("kpts") == [1, 1, 1] and method.get("smearing_eV") == 0.1 and method.get("GPAW_version") == "26.7.0",
        }
        rows.append({
            "label": label,
            "train36_index": index,
            "train_config_type": frame.info["config_type"],
            "historical_source_class": source_class,
            "current_role": "approved_train36",
            "input_path": input_path.relative_to(REPO).as_posix(),
            "output_path": output_path.relative_to(REPO).as_posix(),
            "summary_path": summary_path.relative_to(REPO).as_posix(),
            "source_method": method,
            "scf_iterations": iteration,
            "progress_status": progress.get("status"),
            "execution": {
                "summary_mpi_ranks": summary.get("mpi_ranks", summary.get("execution", {}).get("mpi_ranks")),
                "gpaw_mpi_backend_from_log": backend_match.group(1) if backend_match else None,
                "domain_parallelization_from_log": int(domain_match.group(1)) if domain_match else None,
            },
            "hashes": hashes,
            "geometry": geometry_checks,
            "energy": energies,
            "checks": checks,
            "status": "PASS" if all(checks.values()) else "FAIL",
        })

    result = {
        "status": "PASS" if all(row["status"] == "PASS" for row in rows) else "FAIL",
        "scope": "Reconciles exactly the 12 formerly unresolved approved V19 train36 rows. Reads no valid/dev/test/sealed files beyond the already approved train.extxyz and exact source archives. No GPAW/MPI run.",
        "train36_sha256": train_hash,
        "row_count": len(rows),
        "historical_source_class_counts": {key: sum(row["historical_source_class"] == key for row in rows) for key in sorted({row["historical_source_class"] for row in rows})},
        "method_summary": "All 12 source summaries record GPAW 26.7.0, PBE plane waves, 500 eV, Gamma (1x1x1), and Fermi smearing 0.10 eV. Archived logs show nine runs using the cgpaw backend with domain=4 and three runs using the serial backend with domain=1; this execution difference does not alter the recorded DFT method. The four run scripts call atoms.get_potential_energy() and store it as PW_PBE_energy_eV/energy_eV_cell; none calls get_free_energy(), and the archived extxyz files do not store a separate free-energy value.",
        "limitations": [
            "Exact source method/getter and numerical results are recovered for these 12 labels, resolving the earlier missing-archive lineage flag.",
            "Execution backend differs across the historical archives: three V12 contact calculations used serial GPAW; the remaining nine use four-domain cgpaw. MPI rank count was not separately recorded in the three serial summaries.",
            "This does not establish slab-method consistency: train36 is fully periodic Gamma-point data, while dev3 is TTF 6x6x1 slab data with sigma 0.10 eV and xy dipole correction.",
            "Historical V12 distance-scan and V13 holdout labels are now train36 members; they are not independent V19 development data.",
            "The single V17 rank-recovery row still has an unresolved source-versus-replacement k-point convergence issue; it is outside this 12-row reconciliation."
        ],
        "rows": rows,
    }
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": result["status"], "rows": len(rows), "passed": sum(row["status"] == "PASS" for row in rows), "output": OUT.relative_to(REPO).as_posix()}, indent=2))
    if result["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
