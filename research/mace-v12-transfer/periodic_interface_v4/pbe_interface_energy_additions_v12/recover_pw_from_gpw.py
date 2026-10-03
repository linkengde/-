#!/usr/bin/env python3
"""Finish a converged GPAW checkpoint and export energy/forces with v26 API."""
from pathlib import Path
import hashlib
import json
import os
import re
import shutil
import sys
import time

import numpy as np
from ase.io import read, write

if len(sys.argv) != 6:
    raise SystemExit("usage: recover_pw_from_gpw.py LABEL INPUT.extxyz CHECKPOINT.gpw OUTPUT_DIR ORIGINAL_LOG")
LABEL = sys.argv[1]
SOURCE = Path(sys.argv[2]).resolve()
CHECKPOINT = Path(sys.argv[3]).resolve()
OUT = Path(sys.argv[4]).resolve()
ORIGINAL_LOG = Path(sys.argv[5]).resolve()
if OUT.exists():
    raise SystemExit(f"Refusing to overwrite {OUT}")
if LABEL not in {"AgSi_d2p2", "AgSi_d2p8", "AgC_d2p0", "AgC_d2p4"}:
    raise SystemExit(f"Unexpected label: {LABEL}")
OUT.mkdir(parents=True)

import gpaw_data
import os
os.environ["GPAW_SETUP_PATH"] = str(gpaw_data.datapath())
from gpaw import GPAW

atoms = read(SOURCE)
atoms.pbc = True
started = time.time()
calc = GPAW(str(CHECKPOINT), txt=str(OUT / "gpaw.log"))
checkpoint_atoms = calc.get_atoms()
if len(atoms) != len(checkpoint_atoms):
    raise SystemExit("Atom count differs between source geometry and checkpoint")
if not np.array_equal(atoms.numbers, checkpoint_atoms.numbers):
    raise SystemExit("Atomic numbers differ between source geometry and checkpoint")
if not np.allclose(atoms.cell.array, checkpoint_atoms.cell.array, rtol=0, atol=1e-12):
    raise SystemExit("Cell differs between source geometry and checkpoint")
if not np.allclose(atoms.positions, checkpoint_atoms.positions, rtol=0, atol=1e-12):
    raise SystemExit("Positions differ between source geometry and checkpoint")
convergence_status = {"converged": False}
calc.hooks["converged"] = lambda: convergence_status.__setitem__("converged", True)
atoms.calc = calc
energy = float(atoms.get_potential_energy())
forces = np.asarray(atoms.get_forces(), dtype=float)
original_log = ORIGINAL_LOG.read_text(errors="replace")
converged_match = re.search(r"Converged in\s+(\d+)\s+steps", original_log)
if not convergence_status["converged"] and not converged_match:
    raise SystemExit("Neither recovery calculation nor original log proves SCF convergence")
iterations = int(converged_match.group(1)) if converged_match else int(calc.get_number_of_iterations())
energy_match = re.findall(r"Extrapolated:\s*([-+0-9.eE]+)", original_log)
if energy_match and abs(energy - float(energy_match[-1])) > 1e-4:
    raise SystemExit(f"Recovered energy differs from converged source log by {energy-float(energy_match[-1]):.6g} eV")
if not np.isfinite(energy) or not np.isfinite(forces).all():
    raise SystemExit("Non-finite DFT energy or forces")

contact = "Si" if LABEL.startswith("AgSi") else "C"
ids = np.asarray(atoms.arrays["lammps_id"], dtype=int)
symbols = np.asarray(atoms.get_chemical_symbols())
central = np.asarray(atoms.arrays["central_pair"], dtype=bool)
ia = int(np.flatnonzero(central & (symbols == "Ag"))[0])
ib = int(np.flatnonzero(central & (symbols == contact))[0])
vec = atoms.positions[ib] - atoms.positions[ia]
distance = float(np.linalg.norm(vec))
norms = np.linalg.norm(forces, axis=1)
pair_force = float(np.dot(forces[ib] - forces[ia], vec / distance))
atoms.info["config_type"] = f"{LABEL}_periodic_PW_PBE"
atoms.info["PW_PBE_energy_eV"] = energy
atoms.arrays["PW_PBE_forces"] = forces.copy()
write(OUT / f"{LABEL}_PW_PBE.extxyz", atoms, format="extxyz")
calc.write(str(OUT / "state.gpw"), mode="all")
summary = {
    "label": LABEL,
    "source": str(SOURCE),
    "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    "checkpoint_source": str(CHECKPOINT),
    "checkpoint_source_sha256": hashlib.sha256(CHECKPOINT.read_bytes()).hexdigest(),
    "original_log": str(ORIGINAL_LOG),
    "original_log_sha256": hashlib.sha256(ORIGINAL_LOG.read_bytes()).hexdigest(),
    "atoms": len(atoms),
    "formula": atoms.get_chemical_formula(),
    "pbc": atoms.pbc.tolist(),
    "central_pair": {"Ag_id": int(ids[ia]), f"{contact}_id": int(ids[ib]), "distance_A": distance,
                     "separating_force_eV_A": pair_force},
    "energy_eV_cell": energy,
    "scf_converged": True,
    "scf_iterations": iterations,
    "fmax_eV_A": float(norms.max()),
    "rms_force_eV_A": float(np.sqrt(np.mean(norms**2))),
    "method": {"xc": "PBE", "basis": "plane wave", "cutoff_eV": 500, "kpts": [1, 1, 1],
               "smearing_eV": 0.1, "scf_thresholds": "energy/density/eigenstates 1e-5",
               "mixer": "Pulay beta=0.05, nmaxold=8, weight=100", "GPAW_version": __import__("gpaw").__version__},
    "elapsed_s": time.time() - started,
    "output_extxyz": str(OUT / f"{LABEL}_PW_PBE.extxyz"),
    "state_gpw": str(OUT / "state.gpw"),
    "note": "Recovered from the same-method rolling GPAW checkpoint after the original driver finished SCF but failed during summary bookkeeping; source geometry and checkpoint atom positions/cell were verified identical, the original log's convergence marker and extrapolated energy were checked, and the recovered forces were evaluated from the saved state.",
}
(OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
(OUT / "progress.json").write_text(json.dumps({"label": LABEL, "status": "complete", "iteration": iterations,
                                                   "elapsed_s": summary["elapsed_s"], "scf_converged": True}, indent=2) + "\n")
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
