#!/usr/bin/env python3
"""Diagnostic inference on existing pure-phase controls; never assigns data roles."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import torch
from ase.io import read
from mace.calculators import MACECalculator

REPO = Path(__file__).resolve().parents[6]
PURE = REPO / "research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_pure_phase_controls"
MANIFEST = PURE / "input_manifest.json"
MODEL = REPO / "research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v19_residual_response/runs/v19_train36_dev3_20261023_A/models/MACE_periodic_v19_residual_response.model"
OUT = Path(__file__).resolve().parent / "metrics.json"
SELECTED = (
    "Ag_baseline_k8x8x8_v19",
    "Ti3SiC2_baseline_k10x10x2_v19",
    "Ti3SiC2_baseline_k10x10x4_v19",
    "Ti3SiC2_Ti4f_id0003_z_m020A_k10x10x4_v19",
    "Ti3SiC2_Ti4f_id0003_z_p020A_k10x10x4_v19",
    "Ti3SiC2_C4f_id0009_z_m020A_k10x10x4_v19",
    "Ti3SiC2_C4f_id0009_z_p020A_k10x10x4_v19",
)
CONTROL_ROLE = "numerical_pure_phase_control_not_training_or_independent_validation"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_calc_dir(directory: Path) -> dict[str, object]:
    summary = json.loads((directory / "summary.json").read_text())
    verification = json.loads((directory / "verification.json").read_text())
    checksums = {}
    for line in (directory / "SHA256SUMS.txt").read_text().splitlines():
        if not line.strip():
            continue
        digest, filename = line.split(maxsplit=1)
        artifact = directory / filename.strip()
        checksums[artifact.name] = artifact.is_file() and sha256(artifact) == digest
    assert summary.get("scf_converged") is True, directory
    assert summary.get("dataset_role") == CONTROL_ROLE, directory
    assert verification.get("status") == "PASS", directory
    assert verification.get("checks") and all(verification["checks"].values()), directory
    assert checksums and all(checksums.values()), directory
    return {"scf_converged": True, "archive_verification": "PASS", "sha256_manifest_pass": True}


def main() -> None:
    manifest = json.loads(MANIFEST.read_text())
    records = {record["label"]: record for record in manifest["records"]}
    calc_dirs = sorted(path for path in (PURE / "calculations").iterdir() if path.is_dir())
    inventory = {path.name: verify_calc_dir(path) for path in calc_dirs}
    assert len(inventory) == 26, f"Unexpected archive count: {len(inventory)}"
    torch.set_default_dtype(torch.float64)
    calculator = MACECalculator(model_paths=str(MODEL), device="cpu", default_dtype="float64")
    results = []

    for label in SELECTED:
        record = records[label]
        assert record["role"] == CONTROL_ROLE
        input_path = PURE / record["input"]
        assert sha256(input_path) == record["input_sha256"]
        directory = PURE / "calculations" / label
        result_path = directory / f"{label}_PW_PBE.extxyz"
        source = read(input_path)
        atoms = read(result_path)
        assert np.array_equal(source.numbers, atoms.numbers)
        assert np.allclose(source.positions, atoms.positions, rtol=0.0, atol=1.0e-10)
        assert np.allclose(source.cell.array, atoms.cell.array, rtol=0.0, atol=1.0e-10)
        assert np.array_equal(source.pbc, atoms.pbc)
        for key in ("lammps_id", "atom_id"):
            if key in source.arrays or key in atoms.arrays:
                assert key in source.arrays and key in atoms.arrays
                assert np.array_equal(source.arrays[key], atoms.arrays[key])

        force_key = next(key for key in ("REF_forces", "PW_PBE_forces", "forces") if key in atoms.arrays)
        reference = np.asarray(atoms.arrays[force_key], dtype=float)
        atoms.calc = calculator
        predicted = np.asarray(atoms.get_forces(), dtype=float)
        delta = predicted - reference
        vector_error = np.linalg.norm(delta, axis=1)
        symbols = np.asarray(atoms.get_chemical_symbols())
        rowspecies = {
            element: float(np.sqrt(np.mean(vector_error[symbols == element] ** 2)))
            for element in sorted(set(symbols))
        }
        results.append(
            {
                "label": label,
                "formula": atoms.get_chemical_formula(),
                "n_atoms": len(atoms),
                "role": record["role"],
                "kpts": record["kpts"],
                "pbc": [bool(value) for value in atoms.pbc],
                "input_sha256": sha256(input_path),
                "reference_sha256": sha256(result_path),
                "force_key": force_key,
                "force_vector_RMSE_eV_A": float(np.sqrt(np.mean(vector_error**2))),
                "force_component_RMSE_eV_A": float(np.sqrt(np.mean(delta**2))),
                "force_vector_max_error_eV_A": float(vector_error.max()),
                "force_component_RMSE_by_axis_eV_A": {
                    axis: float(np.sqrt(np.mean(delta[:, index] ** 2)))
                    for index, axis in enumerate("xyz")
                },
                "per_species_force_vector_RMSE_eV_A": rowspecies,
                "geometry_matches_registered_input": True,
                "archive_verification": "PASS",
            }
        )

    result = {
        "status": "DIAGNOSTIC_ONLY",
        "scope": "V19 inference on existing, archive-verified pure-phase numerical controls; no test or sealed data",
        "model": str(MODEL.relative_to(REPO)),
        "model_sha256": sha256(MODEL),
        "inventory": {
            "calculation_archive_count": len(inventory),
            "all_scfs_converged": True,
            "all_archive_sha256_manifests_pass": True,
            "all_records_keep_control_role": True,
        },
        "interpretation_limit": "These controls are not independent interface validation and their role was not changed. Their host CIF source is COD 9009647, shared with the current dev3 family.",
        "references": results,
    }
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
