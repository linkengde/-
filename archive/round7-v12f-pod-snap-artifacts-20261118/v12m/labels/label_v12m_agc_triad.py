#!/workspace/venvs/dft-labels/bin/python
"""GPAW PBE/PW350 label for the v12l POD cold Ag-C-C collapse precursor."""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from ase.io import read, write
from gpaw import FermiDirac, GPAW, Mixer, PW
from gpaw.mpi import world

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "candidates/agc_triad_v12l_t28fs_step18.extxyz"
OUTDIR = HERE / "labels"
STEM = "agc_triad_v12l_t28fs_step18_pw350_gamma"
TARGET_IDS = [414900, 198354, 203868]
OUTDIR.mkdir(exist_ok=True)

atoms = read(SOURCE)
ids = np.asarray(atoms.arrays["source_atom_id"], dtype=int)
target = np.asarray(atoms.arrays["is_DFT_force_target"], dtype=bool)
index = {int(atom_id): i for i, atom_id in enumerate(ids)}
if len(atoms) > 64 or int(target.sum()) != 3 or any(atom_id not in index for atom_id in TARGET_IDS):
    raise RuntimeError("Unexpected Ag-C-C crop or target-row count")
if [atoms[index[i]].symbol for i in TARGET_IDS] != ["Ag", "C", "C"]:
    raise RuntimeError("Target IDs do not map to Ag-C-C")

log = OUTDIR / f"{STEM}.txt"
start = time.time()
atoms.calc = GPAW(
    mode=PW(350),
    xc="PBE",
    kpts=(1, 1, 1),
    symmetry="off",
    occupations=FermiDirac(0.10),
    mixer=Mixer(beta=0.05, nmaxold=8, weight=100),
    convergence={"density": 1e-5, "eigenstates": 1e-5, "energy": 1e-5},
    maxiter=160,
    parallel={"domain": 1, "band": max(1, world.size)},
    txt=str(log),
)
energy = float(atoms.get_potential_energy())
free_energy = float(atoms.get_potential_energy(force_consistent=True))
forces = np.asarray(atoms.get_forces(), dtype=float)
iterations = int(atoms.calc.get_number_of_iterations())
if not np.isfinite(forces).all() or not np.isfinite(energy):
    raise RuntimeError("DFT produced non-finite results")

model_forces = np.asarray(atoms.arrays["POD_v12l_forces"], dtype=float)
center_forces = {str(atom_id): forces[index[atom_id]].tolist() for atom_id in TARGET_IDS}
center_errors = {
    str(atom_id): float(np.linalg.norm(forces[index[atom_id]] - model_forces[index[atom_id]]))
    for atom_id in TARGET_IDS
}
atoms.arrays["REF_forces"] = forces
atoms.info.update(
    {
        "REF_energy": energy,
        "free_energy_eV": free_energy,
        "DFT_method": "GPAW PBE/PW350 Gamma, Fermi-Dirac sigma 0.10 eV",
        "DFT_center_forces_eV_A": json.dumps(center_forces),
        "DFT_POD_v12l_center_force_errors_eV_A": json.dumps(center_errors),
        "DFT_target_ids": json.dumps(TARGET_IDS),
        "label_status": "GPAW single-point; only three Ag-C-C center forces are DFT targets; edge rows excluded from DFT targets",
        "model_data_role": "v12m training candidate; Si01 remains development holdout",
    }
)
atoms.calc = None
output = OUTDIR / f"{STEM}.extxyz"
write(output, atoms, format="extxyz")
log_text = log.read_text(errors="replace")
record = {
    "status": "CONVERGED" if "Converged in " in log_text else "SCF_CONVERGENCE_REVIEW_REQUIRED",
    "state": "v12l_POD_cold_Ag-C-C_precursor_t28fs",
    "center_atom_ids": TARGET_IDS,
    "elements": ["Ag", "C", "C"],
    "center_pair_distances_A": json.loads(atoms.info["center_pair_distances_A"]),
    "natoms": len(atoms),
    "energy_eV": energy,
    "free_energy_eV": free_energy,
    "center_DFT_forces_eV_A": center_forces,
    "center_POD_v12l_force_errors_eV_A": center_errors,
    "max_DFT_force_norm_eV_A": float(np.linalg.norm(forces, axis=1).max()),
    "DFT_target_rows": int(target.sum()),
    "iterations": iterations,
    "elapsed_s": time.time() - start,
    "output": str(output),
    "log": str(log),
    "method": "GPAW PBE/PW350 Gamma, Fermi-Dirac sigma 0.10 eV",
}
(OUTDIR / f"{STEM}.json").write_text(json.dumps(record, indent=2) + "\n")
if world.rank == 0:
    print(json.dumps(record, indent=2), flush=True)
