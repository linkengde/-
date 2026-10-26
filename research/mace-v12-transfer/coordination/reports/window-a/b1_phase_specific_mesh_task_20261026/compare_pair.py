#!/usr/bin/env python3
"""Compare a same-geometry coarse/fine mesh pair using free energy and forces."""
import argparse
import json
from pathlib import Path

import numpy as np
from ase.io import read

HERE = Path(__file__).resolve().parent
TASK = json.loads((HERE / "task_manifest.json").read_text())
parser = argparse.ArgumentParser()
parser.add_argument("coarse_label")
parser.add_argument("fine_label")
parser.add_argument("--output", required=True)
args = parser.parse_args()
ROOT = HERE / TASK["outputs_directory"]
coarse_dir, fine_dir = ROOT / args.coarse_label, ROOT / args.fine_label
coarse = json.loads((coarse_dir / "summary.json").read_text())
fine = json.loads((fine_dir / "summary.json").read_text())
a, b = read(coarse_dir / "result.extxyz"), read(fine_dir / "result.extxyz")
assert coarse["source_sha256"] == fine["source_sha256"]
assert np.array_equal(a.numbers, b.numbers)
assert np.array_equal(a.positions, b.positions)
assert np.array_equal(a.cell.array, b.cell.array)
assert np.array_equal(a.pbc, b.pbc)
fa = np.asarray(a.arrays["PW_PBE_forces"], dtype=float)
fb = np.asarray(b.arrays["PW_PBE_forces"], dtype=float)
delta = fb - fa
per_atom = (fine["free_energy_eV_cell"] - coarse["free_energy_eV_cell"]) * 1000.0 / len(a)
native_per_atom = (fine["native_energy_eV_cell"] - coarse["native_energy_eV_cell"]) * 1000.0 / len(a)
force_rms = float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))
max_atom = float(np.linalg.norm(delta, axis=1).max())
limits = TASK["acceptance"]
checks = {
    "free_energy": abs(per_atom) <= limits["absolute_free_energy_difference_meV_per_atom_max"],
    "force_vector_rms": force_rms <= limits["force_vector_rms_difference_eV_per_A_max"],
    "max_atom_force_vector": max_atom <= limits["maximum_atom_force_vector_difference_eV_per_A_max"],
}
result = {
    "task_id": TASK["task_id"], "pair_id": f"{args.coarse_label}__vs__{args.fine_label}",
    "coarse_label": args.coarse_label, "fine_label": args.fine_label,
    "input_sha256": coarse["source_sha256"], "natoms": len(a),
    "coarse_kpts": coarse["method"]["kpts"], "fine_kpts": fine["method"]["kpts"],
    "free_energy_delta_meV_per_atom": per_atom,
    "native_energy_delta_meV_per_atom_diagnostic_only": native_per_atom,
    "force_vector_rms_delta_eV_per_A": force_rms,
    "maximum_atom_force_vector_delta_eV_per_A": max_atom,
    "checks": checks, "pass": all(checks.values()),
    "interpretation": "Numerical mesh control only; no dataset integration or validation role.",
}
output = HERE / "comparisons" / args.output
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
