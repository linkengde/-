#!/usr/bin/env python3
"""Verify compact fixed-geometry v15 registry holdout results and hashes."""
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parent
manifest = json.loads((ROOT / "input_manifest.json").read_text())
required = set(sys.argv[1:])
known = {r["label"] for r in manifest["records"]}
if required - known:
    raise SystemExit("Unknown labels: " + str(sorted(required - known)))
records = []
missing = []
for rec in manifest["records"]:
    label = rec["label"]
    folder = ROOT / "calculations" / label
    if not (folder / "summary.json").exists():
        missing.append(label)
        continue
    source = ROOT / rec["input"]
    original = read(source)
    output = folder / f"{label}_PW_PBE.extxyz"
    atoms = read(output)
    summary = json.loads((folder / "summary.json").read_text())
    progress = json.loads((folder / "progress.json").read_text())
    energy = float(atoms.info["PW_PBE_energy_eV"])
    forces = np.asarray(atoms.arrays["PW_PBE_forces"], dtype=float)
    pair_indices = np.flatnonzero(np.asarray(original.arrays["central_pair"], dtype=bool))
    ia = int(next(i for i in pair_indices if original[i].symbol == "Ag"))
    ib = int(next(i for i in pair_indices if original[i].symbol != "Ag"))
    distance = float(original.get_distance(ia, ib, mic=True))
    checks = {
        "input_hash": sha256(source.read_bytes()).hexdigest() == rec["input_sha256"] == summary["source_sha256"],
        "convergence": summary.get("scf_converged") is True and progress.get("status") == "complete",
        "mpi_ranks": summary.get("mpi_ranks") == 4,
        "atoms": len(atoms) == len(original) == summary["atoms"] == rec["atoms"],
        "formula": atoms.get_chemical_formula() == original.get_chemical_formula() == rec["formula"],
        "elements": atoms.get_chemical_symbols() == original.get_chemical_symbols(),
        "ids": np.array_equal(atoms.arrays["lammps_id"], original.arrays["lammps_id"]),
        "positions": np.allclose(atoms.positions, original.positions, atol=1e-12, rtol=0),
        "cell": np.allclose(atoms.cell.array, original.cell.array, atol=1e-12, rtol=0),
        "pbc": np.array_equal(atoms.pbc, original.pbc) and np.all(atoms.pbc),
        "central_pair_ids": set(map(int, original.arrays["lammps_id"][pair_indices])) == set(rec["central_pair"]["persistent_ids"]),
        "distance": np.isclose(distance, rec["central_pair"]["input_distance_A"], atol=2e-8, rtol=0)
                    and np.isclose(distance, summary["central_pair"]["distance_A"], atol=1e-10, rtol=0),
        "finite": np.isfinite(energy) and forces.shape == (len(atoms), 3) and np.isfinite(forces).all(),
        "energy": np.isclose(energy, summary["energy_eV_cell"], atol=1e-10, rtol=0),
        "log": f"Converged in {summary['scf_iterations']} steps" in (folder / "gpaw.log").read_text(errors="replace"),
    }
    checks = {k: bool(v) for k, v in checks.items()}
    if not all(checks.values()):
        raise SystemExit("Archive validation failed: " + label + " " + json.dumps(checks))
    hashes = {name: sha256((folder / name).read_bytes()).hexdigest() for name in (output.name, "summary.json", "gpaw.log", "progress.json")}
    verification = {
        "status": "PASS", "label": label, "role": rec["role"], "checks": checks,
        "sha256": hashes, "input_sha256": rec["input_sha256"], "state_gpw_archived": False,
    }
    (folder / "verification.json").write_text(json.dumps(verification, indent=2) + "\n")
    records.append(verification)
    print(label, "archive PASS")
if required & set(missing):
    raise SystemExit("Required output missing: " + str(sorted(required & set(missing))))
if not records:
    raise SystemExit("No completed labels yet; do not report DFT PASS.")
report = {"owner_task": "window-a", "complete": not missing, "records": records, "missing_labels": missing}
dest = ROOT / ("partial_archive_manifest.json" if missing else "archive_manifest.json")
dest.write_text(json.dumps(report, indent=2) + "\n")
print("Verified", len(records), "labels; pending", len(missing))
