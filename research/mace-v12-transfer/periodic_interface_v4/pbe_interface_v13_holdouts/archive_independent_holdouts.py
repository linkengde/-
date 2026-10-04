#!/usr/bin/env python3
"""Verify independent v13 holdout labels and write the hash manifest."""
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parent
manifest = json.loads((ROOT / "input_manifest.json").read_text())
records = []
for source_record in manifest["records"]:
    label = source_record["label"]
    folder = ROOT / "calculations" / label
    input_path = ROOT / source_record["input"]
    input_atoms = read(input_path)
    output_path = folder / f"{label}_PW_PBE.extxyz"
    result = read(output_path)
    summary = json.loads((folder / "summary.json").read_text())
    progress = json.loads((folder / "progress.json").read_text())
    forces = np.asarray(result.arrays["PW_PBE_forces"], dtype=float)
    energy = float(result.info["PW_PBE_energy_eV"])
    pair = np.flatnonzero(np.asarray(input_atoms.arrays["central_pair"], dtype=bool))
    input_distance = float(input_atoms.get_distance(int(pair[0]), int(pair[1]), mic=True))
    # The manifest records the in-memory geometry before extxyz's 8-decimal
    # coordinate serialization. Keep that metadata check within the writer's
    # rounding error; the input hash and output/input position checks stay strict.
    manifest_distance_ok = np.isclose(
        input_distance, source_record["central_pair"]["input_distance_A"], atol=2e-8, rtol=0
    )
    summary_distance_ok = np.isclose(
        input_distance, summary["central_pair"]["distance_A"], atol=1e-10, rtol=0
    )
    checks = {
        "input_hash": sha256(input_path.read_bytes()).hexdigest() == source_record["input_sha256"] == summary["source_sha256"],
        "converged": summary.get("scf_converged") is True and progress.get("status") == "complete",
        "four_mpi_ranks": summary.get("mpi_ranks") == 4,
        "atoms": len(result) == len(input_atoms) == summary["atoms"],
        "elements": result.get_chemical_symbols() == input_atoms.get_chemical_symbols(),
        "ids": np.array_equal(result.arrays["lammps_id"], input_atoms.arrays["lammps_id"]),
        "positions": np.allclose(result.positions, input_atoms.positions, rtol=0, atol=1e-12),
        "cell": np.allclose(result.cell.array, input_atoms.cell.array, rtol=0, atol=1e-12),
        "pbc": bool(np.all(result.pbc)),
        "distance": bool(manifest_distance_ok and summary_distance_ok),
        "finite": bool(np.isfinite(energy) and forces.shape == (len(result), 3) and np.isfinite(forces).all()),
        "energy_consistent": bool(np.isclose(energy, summary["energy_eV_cell"], atol=1e-10, rtol=0)),
        "log": f"Converged in {summary['scf_iterations']} steps" in (folder / "gpaw.log").read_text(errors="replace"),
    }
    print(label, json.dumps(checks, sort_keys=True))
    if not all(checks.values()):
        raise SystemExit(f"Archive validation failed for {label}")
    verification = {
        "status": "PASS",
        "label": label,
        "role": "v13_independent_geometry_check; eligible for v14 training only after v13 scoring",
        "checks": checks,
        "input_sha256": source_record["input_sha256"],
        "sha256": {
            name: sha256((folder / name).read_bytes()).hexdigest()
            for name in (output_path.name, "summary.json", "gpaw.log", "progress.json")
        },
        "state_gpw_archived": False,
    }
    (folder / "verification.json").write_text(json.dumps(verification, indent=2) + "\n")
    records.append({
        "label": label,
        "input_sha256": source_record["input_sha256"],
        "extxyz_sha256": sha256(output_path.read_bytes()).hexdigest(),
        "summary_sha256": sha256((folder / "summary.json").read_bytes()).hexdigest(),
        "atoms": summary["atoms"],
        "central_pair": summary["central_pair"],
        "scf_iterations": summary["scf_iterations"],
        "energy_eV_cell": energy,
        "state_gpw_archived": False,
    })
(ROOT / "archive_manifest.json").write_text(json.dumps({"records": records}, indent=2) + "\n")
print(f"Independent archive PASS: {len(records)} holdouts")
