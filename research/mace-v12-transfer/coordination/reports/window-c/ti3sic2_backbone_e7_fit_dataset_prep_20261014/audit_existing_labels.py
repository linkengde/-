#!/usr/bin/env python3
"""Read-only hash/completeness audit of the five archived B bulk DFT labels."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read


ROOT = Path(__file__).resolve().parent
CALCS = ROOT.parents[3] / "periodic_interface_v4" / "pbe_interface_v19_pure_phase_controls" / "calculations"
B_RESULTS = ROOT.parents[1] / "window-b" / "v19_Ti3SiC2_backbone_force_response_execution_window_b" / "final_results.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(label: str, expected_hash: str, expected_n: int, output_dir: Path,
           output_file: str, summary: dict) -> dict:
    path = output_dir / output_file
    actual_hash = sha256(path)
    atoms = read(path, format="extxyz")
    properties = atoms.info
    forces = atoms.arrays.get("PW_PBE_forces")
    force_key = "PW_PBE_forces"
    energy_key = "PW_PBE_energy_eV"
    energy = float(properties[energy_key]) if energy_key in properties else float("nan")
    expected_energy = float(summary["energy_eV_cell"])
    checks = {
        "file_sha256_match": actual_hash == expected_hash,
        "natoms_match": len(atoms) == expected_n,
        "forces_present": force_key in atoms.arrays and atoms.arrays[force_key].shape == (expected_n, 3),
        "native_energy_present": energy_key in properties,
        "native_energy_matches_summary": bool(np.isfinite(energy) and abs(energy - expected_energy) <= 1e-10),
        "forces_finite": bool(forces is not None and np.isfinite(forces).all()),
        "free_energy_recorded": bool(np.isfinite(float(summary["free_energy_eV_cell"]))),
        "scf_converged": summary.get("scf_converged", summary.get("status") == "converged"),
    }
    if not all(checks.values()):
        raise SystemExit(f"Existing B label failed verification: {label}: {checks}")
    return {
        "label": label,
        "extxyz_sha256": actual_hash,
        "natoms": len(atoms),
        "force_property": force_key,
        "native_energy_property": energy_key,
        "force_finite_components": bool(np.isfinite(forces).all()),
        "energy_eV_cell": summary.get("energy_eV_cell"),
        "free_energy_eV_cell": summary.get("free_energy_eV_cell"),
        "stress_label": "NOT_RECORDED_IN_B_FINAL_RESULTS",
        "checks": checks,
    }


def main() -> None:
    results = json.loads(B_RESULTS.read_text(encoding="utf-8"))
    records = []

    baseline = CALCS / "Ti3SiC2_baseline_k10x10x4_v19"
    baseline_file = "Ti3SiC2_baseline_k10x10x4_v19_PW_PBE.extxyz"
    baseline_manifest = {
        row.split()[1]: row.split()[0]
        for row in (baseline / "SHA256SUMS.txt").read_text(encoding="utf-8").splitlines()
        if row.strip()
    }
    baseline_summary = json.loads((baseline / "summary.json").read_text(encoding="utf-8"))
    records.append(verify("Ti3SiC2_baseline_k10x10x4_v19", baseline_manifest[baseline_file], 48,
                          baseline, baseline_file, baseline_summary))

    for label, item in results["points"].items():
        output_file = f"{label}_PW_PBE.extxyz"
        summary = {
            "status": item.get("status"),
            "scf_converged": item.get("status") == "converged" and item.get("archive_verification") == "PASS",
            "energy_eV_cell": item.get("energy_eV_cell"),
            "free_energy_eV_cell": item.get("free_energy_eV_cell"),
        }
        records.append(verify(label, item["archive_file_sha256"][output_file], 48,
                              CALCS / label, output_file, summary))

    report = {
        "scope": "read-only verification of existing B DFT output files; no calculation invoked",
        "source_commit": "34b8be988015af96d3fd8791fc041f6cc811e2b7",
        "method_manifest_sha256": results["method"]["method_manifest_sha256"],
        "verified_record_count": len(records),
        "all_records_have_48_atoms_energy_and_forces": True,
        "stress_completeness": "No stress recorded for these five bulk records in the reviewed source; do not use as stress labels.",
        "records": records,
    }
    (ROOT / "B_EXISTING_LABEL_VERIFICATION.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"verified {len(records)} existing 48-atom B DFT energy/force records")


if __name__ == "__main__":
    main()
