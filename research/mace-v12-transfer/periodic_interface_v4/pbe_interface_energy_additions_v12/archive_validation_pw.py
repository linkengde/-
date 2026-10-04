#!/usr/bin/env python3
"""Verify and manifest the three independent validation PW-PBE labels."""
from pathlib import Path
import hashlib
import json
import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parent
BASE = ROOT / "calculations/validation_distance_scans"
inputs = json.loads((ROOT / "validation_distance_input_manifest.json").read_text())
records = []
for input_record in inputs["records"]:
    label = input_record["label"]
    folder = BASE / label
    summary = json.loads((folder / "summary.json").read_text())
    progress = json.loads((folder / "progress.json").read_text())
    if summary.get("label") != label or summary.get("scf_converged") is not True:
        raise SystemExit(f"Non-converged validation result: {label}")
    if progress.get("status") != "complete" or progress.get("scf_converged") is not True:
        raise SystemExit(f"Incomplete progress record: {label}")
    source = ROOT / input_record["input"]
    result = read(folder / f"{label}_PW_PBE.extxyz")
    input_atoms = read(source)
    forces = np.asarray(result.arrays["PW_PBE_forces"], dtype=float)
    energy = float(result.info["PW_PBE_energy_eV"])
    central = np.flatnonzero(np.asarray(input_atoms.arrays["central_pair"], dtype=bool))
    saved_distance = float(np.linalg.norm(input_atoms.positions[central[1]] - input_atoms.positions[central[0]]))
    checks = {
        "input_hash": hashlib.sha256(source.read_bytes()).hexdigest() == input_record["input_sha256"] == summary["source_sha256"],
        "atoms": len(result) == len(input_atoms) == summary["atoms"],
        "elements": result.get_chemical_symbols() == input_atoms.get_chemical_symbols(),
        "ids": np.array_equal(result.arrays["lammps_id"], input_atoms.arrays["lammps_id"]),
        "positions": np.allclose(result.positions, input_atoms.positions, rtol=0, atol=1e-12),
        "cell": np.allclose(result.cell.array, input_atoms.cell.array, rtol=0, atol=1e-12),
        "pbc": bool(np.all(result.pbc)),
        "central_pair_distance": bool(np.isclose(summary["central_pair"]["distance_A"], saved_distance, rtol=0, atol=1e-10) and np.isclose(saved_distance, input_record["central_pair"]["input_distance_A"], rtol=0, atol=1e-10)),
        "finite": bool(np.isfinite(energy) and np.isfinite(forces).all()),
        "energy": bool(np.isclose(energy, summary["energy_eV_cell"], rtol=0, atol=1e-10)),
        "log": f"Converged in {summary['scf_iterations']} steps" in (folder / "gpaw.log").read_text(errors="replace"),
        "files": all((folder / name).is_file() and (folder / name).stat().st_size for name in (f"{label}_PW_PBE.extxyz", "summary.json", "gpaw.log", "progress.json")),
    }
    print(label, json.dumps(checks, sort_keys=True))
    if not all(checks.values()):
        raise SystemExit(f"Validation archive check failed: {label}")
    records.append({
        "label": label,
        "atoms": summary["atoms"],
        "central_pair": summary["central_pair"],
        "scf_iterations": summary["scf_iterations"],
        "energy_eV_cell": energy,
        "input_sha256": summary["source_sha256"],
        "extxyz_sha256": hashlib.sha256((folder / f"{label}_PW_PBE.extxyz").read_bytes()).hexdigest(),
        "summary_sha256": hashlib.sha256((folder / "summary.json").read_bytes()).hexdigest(),
        "state_gpw_archived": False,
    })
(BASE / "archive_manifest.json").write_text(json.dumps({"records": records}, indent=2) + "\n")
print(f"Independent validation archive PASS: {len(records)} labels")
