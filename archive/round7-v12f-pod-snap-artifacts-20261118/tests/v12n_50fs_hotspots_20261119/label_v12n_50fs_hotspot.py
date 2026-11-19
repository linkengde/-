#!/workspace/venvs/dft-labels/bin/python
"""Single-point GPAW labels for the two severe local environments at POD v12n 50 fs."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from ase.io import read, write
from gpaw import FermiDirac, GPAW, Mixer, PW
from gpaw.mpi import world

HERE = Path(__file__).resolve().parent
PAIR = sys.argv[1] if len(sys.argv) > 1 else ""
KNOWN = {
    "tic_94075_99595": ([94075, 99595], ["Ti", "C"]),
    "cc_38346_33016": ([38346, 33016], ["C", "C"]),
}
if PAIR not in KNOWN:
    raise SystemExit(f"Usage: {Path(sys.argv[0]).name} <{'|'.join(KNOWN)}>")
TARGET_IDS, ELEMENTS = KNOWN[PAIR]
STEM = f"{PAIR}_t50fs_pw350_gamma"
SOURCE = HERE / "candidates" / f"{PAIR}_t50fs.extxyz"
OUTDIR = HERE / "labels"
OUTDIR.mkdir(exist_ok=True)
atoms = read(SOURCE)
ids = np.asarray(atoms.arrays["source_atom_id"], dtype=int)
target = np.asarray(atoms.arrays["is_DFT_force_target"], dtype=bool)
index = {int(atom_id): i for i, atom_id in enumerate(ids)}
if len(atoms) > 96 or int(target.sum()) != len(TARGET_IDS):
    raise RuntimeError("Unexpected crop size or target count")
if any(atom_id not in index for atom_id in TARGET_IDS):
    raise RuntimeError("Target atom missing from crop")
if [atoms[index[atom_id]].symbol for atom_id in TARGET_IDS] != ELEMENTS:
    raise RuntimeError("Target element mapping mismatch")

log = OUTDIR / f"{STEM}.txt"
start = time.time()
atoms.calc = GPAW(
    mode=PW(350), xc="PBE", kpts=(1, 1, 1), symmetry="off",
    occupations=FermiDirac(0.10),
    mixer=Mixer(beta=0.05, nmaxold=8, weight=100),
    convergence={"density": 1e-5, "eigenstates": 1e-5, "energy": 1e-5},
    maxiter=160, parallel={"domain": max(1, world.size), "band": 1},
    txt=str(log),
)
energy = float(atoms.get_potential_energy())
free_energy = float(atoms.get_potential_energy(force_consistent=True))
forces = np.asarray(atoms.get_forces(), float)
iterations = int(atoms.calc.get_number_of_iterations())
if not np.isfinite(forces).all() or not np.isfinite(energy):
    raise RuntimeError("Nonfinite DFT energy or force")

model_forces = np.asarray(atoms.arrays["POD_v12n_forces"], float)
center_forces = {str(atom_id): forces[index[atom_id]].tolist() for atom_id in TARGET_IDS}
center_errors = {
    str(atom_id): float(np.linalg.norm(forces[index[atom_id]] - model_forces[index[atom_id]]))
    for atom_id in TARGET_IDS
}
atoms.arrays["REF_forces"] = forces
atoms.info.update({
    "REF_energy": energy,
    "free_energy_eV": free_energy,
    "DFT_method": "GPAW PBE/PW350 Gamma, Fermi-Dirac sigma 0.10 eV",
    "DFT_center_forces_eV_A": json.dumps(center_forces),
    "DFT_POD_v12n_center_force_errors_eV_A": json.dumps(center_errors),
    "DFT_target_ids": json.dumps(TARGET_IDS),
    "label_status": "GPAW single-point; central atom force rows are targets; crop-edge rows are not fit targets",
    "model_data_role": "v12o training candidate; validation and Si01 holdout unchanged",
})
atoms.calc = None
output = OUTDIR / f"{STEM}.extxyz"
write(output, atoms, format="extxyz")
log_text = log.read_text(errors="replace")
status = "CONVERGED" if "Converged in " in log_text else "SCF_CONVERGENCE_REVIEW_REQUIRED"
record = {
    "status": status,
    "state": "POD_v12n_cold_50fs",
    "candidate": PAIR,
    "center_atom_ids": TARGET_IDS,
    "elements": ELEMENTS,
    "center_pair_distance_A": float(atoms.info["center_pair_distance_A"]),
    "natoms": len(atoms),
    "energy_eV": energy,
    "free_energy_eV": free_energy,
    "center_DFT_forces_eV_A": center_forces,
    "center_POD_v12n_force_errors_eV_A": center_errors,
    "max_DFT_force_norm_eV_A": float(np.linalg.norm(forces, axis=1).max()),
    "DFT_target_rows": len(TARGET_IDS),
    "iterations": iterations,
    "elapsed_s": time.time() - start,
    "output": str(output),
    "log": str(log),
    "method": "GPAW PBE/PW350 Gamma, Fermi-Dirac sigma 0.10 eV",
}
(OUTDIR / f"{STEM}.json").write_text(json.dumps(record, indent=2) + "\n")
if world.rank == 0:
    print(json.dumps(record, indent=2), flush=True)
