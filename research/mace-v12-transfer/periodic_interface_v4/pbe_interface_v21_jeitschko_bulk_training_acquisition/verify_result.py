#!/usr/bin/env python3
"""Geometry, finite-label, and calculation-metadata checks for one DFT result."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from ase.io import read


ROOT = Path(__file__).resolve().parent
if len(sys.argv) != 2:
    raise SystemExit("usage: verify_result.py LABEL")
label = sys.argv[1]
manifest = json.loads((ROOT / "input_manifest.json").read_text())
record = next(item for item in manifest["records"] if item["label"] == label)
folder = ROOT / "calculations" / label
input_path = ROOT / record["input"]
output_path = folder / f"{label}_PW_PBE.extxyz"
summary = json.loads((folder / "summary.json").read_text())
progress = json.loads((folder / "progress.json").read_text())
source = read(input_path)
result = read(output_path)
forces = result.arrays["PW_PBE_forces"]
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()

checks = {
    "registered_input_hash": sha(input_path) == record["input_sha256"] == summary["input_sha256"],
    "source_cif_hash": record["source_cif_sha256"] == summary["source_cif_sha256"],
    "label_identity": summary["label"] == progress["label"] == label,
    "atom_count_formula": len(source) == len(result) == record["atom_count"] and source.get_chemical_formula() == result.get_chemical_formula(),
    "ordered_symbols": source.get_chemical_symbols() == result.get_chemical_symbols(),
    "exact_positions": np.array_equal(source.positions, result.positions),
    "exact_cell": np.array_equal(source.cell.array, result.cell.array),
    "pbc": source.pbc.tolist() == result.pbc.tolist() == [True, True, True],
    "local_ids": np.array_equal(source.arrays["local_new_id"], result.arrays["local_new_id"]),
    "finite_energy_and_forces": np.isfinite(summary["native_extrapolated_energy_eV_cell"]) and np.isfinite(summary["free_energy_eV_cell"]) and forces.shape == (record["atom_count"], 3) and np.isfinite(forces).all(),
    "converged": summary["scf_converged"] is True and progress["scf_converged"] is True and progress["status"] == "complete" and summary["scf_iterations"] == progress["iteration"],
    "four_mpi_ranks": summary["mpi_ranks"] == progress["mpi_ranks"] == 4,
    "fixed_method": summary["method"] == record["method"],
    "role_not_training_integrated": summary["training_integration"].startswith("not performed"),
}
checks = {key: bool(value) for key, value in checks.items()}
if not all(checks.values()):
    raise SystemExit(json.dumps({"status": "FAIL", "checks": checks}, indent=2))
print(json.dumps({
    "status": "PASS",
    "label": label,
    "checks": checks,
    "input_sha256": sha(input_path),
    "source_cif_sha256": record["source_cif_sha256"],
    "native_extrapolated_energy_eV_cell": summary["native_extrapolated_energy_eV_cell"],
    "free_energy_eV_cell": summary["free_energy_eV_cell"],
    "scf_iterations": summary["scf_iterations"],
    "max_force_eV_A": summary["max_force_eV_A"],
}, indent=2))
