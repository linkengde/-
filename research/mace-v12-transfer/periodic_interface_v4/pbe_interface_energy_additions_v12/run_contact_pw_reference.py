#!/usr/bin/env python3
"""Run a GPAW PW-PBE single point on an existing AgSi or AgC interface geometry."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import sys
import time
import numpy as np
from ase.io import read, write

if len(sys.argv) != 4:
    raise SystemExit("usage: run_contact_pw_reference.py LABEL INPUT.extxyz OUTPUT_DIR")
LABEL = sys.argv[1]
SOURCE = Path(sys.argv[2]).resolve()
OUT = Path(sys.argv[3]).resolve()
OUT.mkdir(parents=True, exist_ok=False)
if LABEL not in {"AgSi_d2p2", "AgSi_d2p8", "AgC_d2p0", "AgC_d2p4"}:
    raise SystemExit(f"Unexpected label: {LABEL}")

import gpaw_data
os.environ["GPAW_SETUP_PATH"] = str(gpaw_data.datapath())
from gpaw import GPAW, PW
from gpaw.mixer import Mixer
from gpaw.occupations import FermiDirac

atoms = read(SOURCE)
atoms.pbc = True
contact = "Si" if LABEL.startswith("AgSi") else "C"
ids = np.asarray(atoms.arrays["lammps_id"], dtype=int)
symbols = np.asarray(atoms.get_chemical_symbols())
central = np.asarray(atoms.arrays["central_pair"], dtype=bool)
ia = int(np.flatnonzero(central & (symbols == "Ag"))[0])
ib = int(np.flatnonzero(central & (symbols == contact))[0])
vec = atoms.positions[ib] - atoms.positions[ia]
distance = float(np.linalg.norm(vec))
started = time.time()
progress_path = OUT / "progress.json"
state_path = OUT / "state.gpw"  # rolling checkpoint; never create per-iteration multi-GB copies

calc = GPAW(
    mode=PW(500), xc="PBE", kpts=(1, 1, 1), occupations=FermiDirac(0.1),
    mixer=Mixer(0.05, 8, 100),
    convergence={"energy": 1e-5, "density": 1e-5, "eigenstates": 1e-5},
    maxiter=160, txt=str(OUT / "gpaw.log"),
)
atoms.calc = calc

def checkpoint(ctx):
    niter = int(ctx.niter)
    if niter % 10 == 0:
        progress_path.write_text(json.dumps({"label": LABEL, "iteration": niter, "elapsed_s": time.time() - started}, indent=2) + "\n")
    if niter % 20 == 0:
        calc.write(str(state_path), mode="all")

calc.hooks["scf_step"] = checkpoint
convergence_status = {"converged": False}
calc.hooks["converged"] = lambda: convergence_status.__setitem__("converged", True)
energy = float(atoms.get_potential_energy())
forces = np.asarray(atoms.get_forces(), dtype=float)
converged = bool(convergence_status["converged"])
scf_iterations = int(calc.get_number_of_iterations())
norms = np.linalg.norm(forces, axis=1)
pair_force = float(np.dot(forces[ib] - forces[ia], vec / distance))
atoms.info["config_type"] = f"{LABEL}_periodic_PW_PBE"
atoms.info["PW_PBE_energy_eV"] = energy
atoms.arrays["PW_PBE_forces"] = forces.copy()
write(OUT / f"{LABEL}_PW_PBE.extxyz", atoms, format="extxyz")
calc.write(str(state_path), mode="all")
summary = {
    "label": LABEL, "source": str(SOURCE), "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    "atoms": len(atoms), "formula": atoms.get_chemical_formula(), "pbc": atoms.pbc.tolist(),
    "central_pair": {"Ag_id": int(ids[ia]), f"{contact}_id": int(ids[ib]), "distance_A": distance,
                     "separating_force_eV_A": pair_force},
    "energy_eV_cell": energy, "scf_converged": converged, "scf_iterations": scf_iterations,
    "fmax_eV_A": float(norms.max()), "rms_force_eV_A": float(np.sqrt(np.mean(norms**2))),
    "method": {"xc": "PBE", "basis": "plane wave", "cutoff_eV": 500, "kpts": [1, 1, 1],
               "smearing_eV": 0.1, "scf_thresholds": "energy/density/eigenstates 1e-5",
               "mixer": "Pulay beta=0.05, nmaxold=8, weight=100", "GPAW_version": __import__("gpaw").__version__},
    "elapsed_s": time.time() - started, "output_extxyz": str(OUT / f"{LABEL}_PW_PBE.extxyz"),
    "state_gpw": str(state_path),
    "note": "Unrelaxed single-point energy/force on a fixed existing interface geometry; external midpoint AgSi/AgC structures remain excluded from training.",
}
(OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
progress_path.write_text(json.dumps({"label": LABEL, "status": "complete", "iteration": scf_iterations, "elapsed_s": summary["elapsed_s"], "scf_converged": converged}, indent=2) + "\n")
if converged:
    archive_root = Path(__file__).resolve().parent / "calculations"
    archive_root.mkdir(exist_ok=True)
    archive_dir = archive_root / LABEL
    if archive_dir.exists():
        raise SystemExit(f"Refusing to overwrite persistent result archive {archive_dir}")
    staging = archive_root / f".{LABEL}.staging-{os.getpid()}"
    staging.mkdir()
    for name in [f"{LABEL}_PW_PBE.extxyz", "summary.json", "gpaw.log", "progress.json"]:
        shutil.copy2(OUT / name, staging / name)
    staging.rename(archive_dir)
print(json.dumps(summary, indent=2), flush=True)
if not converged:
    raise SystemExit(f"SCF did not converge for {LABEL}")
