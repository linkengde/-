#!/usr/bin/env python3
"""Verify one registered result archive against its exact geometry and run record."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parent
if len(sys.argv) != 2:
    raise SystemExit("usage: verify_preflight.py EXACT_REGISTERED_LABEL")
label = sys.argv[1]
manifest = json.loads((ROOT / "input_manifest.json").read_text())
job = next(row for row in manifest["preflight_jobs"] if row["label"] == label)
folder = ROOT / "calculations" / label
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
source, output = read(ROOT / job["input"]), read(folder / f"{label}_PW_PBE.extxyz")
summary = json.loads((folder / "summary.json").read_text())
progress = json.loads((folder / "progress.json").read_text())
log = (folder / "gpaw.log").read_text()
forces = output.arrays["PW_PBE_forces"]
checks = {
    "input_hash": sha(ROOT / job["input"]) == job["input_sha256"] == summary["source_sha256"],
    "identity": summary["label"] == progress["label"] == label,
    "count_formula_symbols": len(source) == len(output) == 52 and source.get_chemical_formula() == output.get_chemical_formula() == "C16Ag4Si8Ti24" and source.get_chemical_symbols() == output.get_chemical_symbols(),
    "geometry_and_pbc": np.array_equal(source.positions, output.positions) and np.array_equal(source.cell.array, output.cell.array) and np.array_equal(source.pbc, output.pbc),
    "ids_and_contact_markers": np.array_equal(source.arrays["lammps_id"], output.arrays["lammps_id"]) and np.array_equal(source.arrays["central_pair"], output.arrays["central_pair"]),
    "convergence": summary["scf_converged"] is True and progress["scf_converged"] is True and progress["status"] == "complete" and f"Converged in {summary['scf_iterations']} steps" in log,
    "method": summary["method"]["xc"] == "PBE" and summary["method"]["cutoff_eV"] == 500 and summary["method"]["kpts"] == job["kpts"] and summary["method"]["smearing_eV"] == job["smearing_eV"] and summary["method"]["poissonsolver"] == {"dipolelayer": "xy"},
    "energy_force_consistency_recorded_separately": output.info["PW_PBE_energy_eV"] == summary["native_extrapolated_energy_eV_cell"] and output.info["PW_PBE_free_energy_eV"] == summary["free_energy_eV_cell"] and summary["training_integration"] == "not performed",
    "finite_forces_and_energies": forces.shape == (len(source), 3) and np.isfinite(forces).all() and np.isfinite(summary["native_extrapolated_energy_eV_cell"]) and np.isfinite(summary["free_energy_eV_cell"]),
    "force_summaries": np.isclose(summary["max_force_eV_A"], np.linalg.norm(forces, axis=1).max(), rtol=0, atol=1e-12) and np.isclose(summary["rms_force_vector_eV_A"], np.sqrt(np.mean(np.sum(forces * forces, axis=1))), rtol=0, atol=1e-12),
    "no_large_checkpoint": not list(folder.glob("*.gpw")),
}
checks = {k: bool(v) for k, v in checks.items()}
if not all(checks.values()):
    raise SystemExit(json.dumps(checks, indent=2))
inventory = {}
for line in (folder / "SHA256SUMS.txt").read_text().splitlines():
    digest, name = line.split("  ", 1)
    path = folder / name
    assert name not in inventory and path.is_file() and sha(path) == digest
    inventory[name] = digest
actual = {p.name for p in folder.iterdir() if p.is_file() and p.name != "SHA256SUMS.txt"}
assert actual == set(inventory)
result = {"status": "PASS", "label": label, "checks": checks, "input_sha256": job["input_sha256"], "archive_sha256": inventory}
print(json.dumps(result, indent=2))
