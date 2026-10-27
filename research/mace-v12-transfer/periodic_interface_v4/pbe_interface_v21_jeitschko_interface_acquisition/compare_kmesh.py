#!/usr/bin/env python3
"""Compare registered 6x6x1 and 8x8x1 same-geometry slab preflight results."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from ase.geometry import find_mic
from ase.io import read

ROOT = Path(__file__).resolve().parent
manifest = json.loads((ROOT / "input_manifest.json").read_text())
labels = {row["label"]: row for row in manifest["preflight_jobs"]}
fine = "AgSi_JN1967_registryB_k6x6x1_sigma0p10_preflight_v21"
control = "AgSi_JN1967_registryB_k8x8x1_sigma0p10_control_v21"
folders = {label: ROOT / "calculations" / label for label in (fine, control)}
for label in (fine, control):
    import subprocess, sys
    subprocess.run([sys.executable, str(ROOT / "verify_preflight.py"), label], check=True, stdout=subprocess.DEVNULL)
summaries = {label: json.loads((folder / "summary.json").read_text()) for label, folder in folders.items()}
atoms = {label: read(folder / f"{label}_PW_PBE.extxyz") for label, folder in folders.items()}
a, b = atoms[fine], atoms[control]
assert np.array_equal(a.numbers, b.numbers) and np.array_equal(a.positions, b.positions)
assert np.array_equal(a.cell.array, b.cell.array) and np.array_equal(a.pbc, b.pbc)
symbols = np.asarray(a.get_chemical_symbols())
marked = a.arrays["central_pair"].astype(bool)
ia = np.flatnonzero(marked & (symbols == "Ag"))
ix = np.flatnonzero(marked & (symbols == "Si"))
assert len(ia) == len(ix) == 1
ia, ix = int(ia[0]), int(ix[0])
unit = find_mic(a.positions[ix] - a.positions[ia], a.cell, a.pbc)[0]
unit /= np.linalg.norm(unit)
delta = b.arrays["PW_PBE_forces"] - a.arrays["PW_PBE_forces"]
metrics = {
    "native_energy_delta_meV_atom_8x8_minus_6x6": (summaries[control]["native_extrapolated_energy_eV_cell"] - summaries[fine]["native_extrapolated_energy_eV_cell"]) * 1000 / len(a),
    "free_energy_delta_meV_atom_8x8_minus_6x6": (summaries[control]["free_energy_eV_cell"] - summaries[fine]["free_energy_eV_cell"]) * 1000 / len(a),
    "force_vector_difference_RMS_eV_A": float(np.sqrt(np.mean(np.sum(delta * delta, axis=1)))),
    "maximum_atom_force_vector_difference_eV_A": float(np.linalg.norm(delta, axis=1).max()),
    "marked_Ag_Si_pair_projection_delta_eV_A": float(np.dot(delta[ix] - delta[ia], unit)),
}
budgets = {"native_energy_meV_per_atom": 2.0, "free_energy_meV_per_atom": 2.0, "force_vector_rms_eV_A": 0.01, "marked_pair_projection_eV_A": 0.02}
checks = {
    "native_energy": abs(metrics["native_energy_delta_meV_atom_8x8_minus_6x6"]) <= budgets["native_energy_meV_per_atom"],
    "free_energy": abs(metrics["free_energy_delta_meV_atom_8x8_minus_6x6"]) <= budgets["free_energy_meV_per_atom"],
    "force_vector_RMS": metrics["force_vector_difference_RMS_eV_A"] <= budgets["force_vector_rms_eV_A"],
    "marked_pair_projection": abs(metrics["marked_Ag_Si_pair_projection_delta_eV_A"]) <= budgets["marked_pair_projection_eV_A"],
}
result = {"status": "PASS_WITHIN_PROVISIONAL_BUDGET" if all(checks.values()) else "SENSITIVITY_EXCEEDS_BUDGET",
          "labels": [fine, control], "same_geometry_sha256": labels[fine]["input_sha256"],
          "metrics": metrics, "budgets": budgets, "checks": checks,
          "archive_sha256": {label: {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in folder.iterdir() if p.is_file() and p.name != "SHA256SUMS.txt"} for label, folder in folders.items()},
          "role": "8x8x1 is numerical control only; neither endpoint changes the V19 dataset"}
out = ROOT / "comparisons"
out.mkdir(exist_ok=True)
path = out / f"comparison_{fine}__{control}.json"
path.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
