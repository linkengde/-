#!/usr/bin/env python3
"""Run one registered fixed-geometry GPAW label for the v21 slab preflight."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time

import numpy as np
from gpaw.mpi import world

ROOT = Path(__file__).resolve().parent
if len(sys.argv) != 3:
    raise SystemExit("usage: run_preflight_label.py LABEL SCRATCH_DIR")
LABEL, SCRATCH_ARG = sys.argv[1:]
manifest = json.loads((ROOT / "input_manifest.json").read_text())
job = next(row for row in manifest["preflight_jobs"] if row["label"] == LABEL)
assert job["launch_enabled"] is True
assert world.size == 4, f"expected four MPI ranks, got {world.size}"
source = ROOT / job["input"]
assert hashlib.sha256(source.read_bytes()).hexdigest() == job["input_sha256"]
scratch = Path(SCRATCH_ARG).resolve()
if world.rank == 0:
    scratch.mkdir(parents=True, exist_ok=False)
world.barrier()

import gpaw_data

os.environ["GPAW_SETUP_PATH"] = str(gpaw_data.datapath())
from ase.io import read, write
from gpaw import GPAW, PW
from gpaw.mixer import Mixer
from gpaw.occupations import FermiDirac

atoms = read(source)
assert len(atoms) == 52 and atoms.pbc.tolist() == [True, True, False]
assert not any(k in atoms.info for k in ("energy", "REF_energy", "PW_PBE_energy_eV"))
assert not any(k in atoms.arrays for k in ("forces", "REF_forces", "PW_PBE_forces"))

started = time.time()
progress_path = scratch / "progress.json"
if world.rank == 0:
    progress_path.write_text(json.dumps({"label": LABEL, "status": "running", "iteration": 0,
        "source_sha256": job["input_sha256"], "mpi_ranks": world.size,
        "started_epoch_s": started}, indent=2) + "\n")
calc = GPAW(
    mode=PW(500), xc="PBE", kpts=tuple(job["kpts"]),
    poissonsolver={"dipolelayer": "xy"},
    occupations=FermiDirac(job["smearing_eV"]),
    mixer=Mixer(0.05, 8, 100), parallel={"sl_auto": True},
    convergence={"energy": 1e-5, "density": 1e-5, "eigenstates": 1e-5},
    maxiter=160, txt=str(scratch / "gpaw.log"),
)
atoms.calc = calc
converged = {"value": False}

def scf_step(context):
    iteration = int(context.niter)
    if world.rank == 0 and iteration % 10 == 0:
        progress_path.write_text(json.dumps({"label": LABEL, "status": "running", "iteration": iteration,
            "elapsed_s": time.time() - started, "source_sha256": job["input_sha256"],
            "mpi_ranks": world.size}, indent=2) + "\n")
    if iteration % 20 == 0:
        calc.write(str(scratch / "state.gpw"), mode="all")

calc.hooks["scf_step"] = scf_step
calc.hooks["converged"] = lambda: converged.__setitem__("value", True)
native = float(atoms.get_potential_energy(force_consistent=False))
free = float(calc.results["free_energy"])
forces = np.asarray(atoms.get_forces(), dtype=float)
iterations = int(calc.get_number_of_iterations())
if not converged["value"]:
    raise SystemExit(f"SCF failed to converge within {iterations} iterations")
if not np.isfinite(native) or not np.isfinite(free) or not np.isfinite(forces).all():
    raise SystemExit("nonfinite DFT result")

if world.rank == 0:
    out = atoms.copy()
    out.calc = None
    out.info.update({"config_type": LABEL, "PW_PBE_energy_eV": native,
                     "PW_PBE_free_energy_eV": free})
    out.arrays["PW_PBE_forces"] = forces.copy()
    write(scratch / f"{LABEL}_PW_PBE.extxyz", out, format="extxyz")
    summary = {
        "label": LABEL, "geometry_label": job["geometry_label"],
        "input": str(source.relative_to(ROOT)), "source_sha256": job["input_sha256"],
        "source_family": "Jeitschko_Nowotny_1967_single_crystal",
        "source_doi": "10.1007/BF00899949", "atom_count": len(atoms),
        "formula": atoms.get_chemical_formula(), "cell_angstrom": atoms.cell.array.tolist(),
        "pbc": atoms.pbc.tolist(), "native_extrapolated_energy_eV_cell": native,
        "free_energy_eV_cell": free,
        "forces_method": "GPAW ASE forces at Fermi-Dirac width 0.10 eV; free energy retained separately",
        "energy_convention": "native extrapolated energy (force_consistent=False); free energy separately recorded",
        "scf_converged": True, "scf_iterations": iterations,
        "max_force_eV_A": float(np.linalg.norm(forces, axis=1).max()),
        "rms_force_vector_eV_A": float(np.sqrt(np.mean(np.sum(forces**2, axis=1)))),
        "method": {"xc": "PBE", "basis": "plane wave", "cutoff_eV": 500,
                   "kpts": list(job["kpts"]), "smearing_eV": job["smearing_eV"],
                   "pbc": atoms.pbc.tolist(), "poissonsolver": {"dipolelayer": "xy"},
                   "mixer": {"type": "Pulay", "beta": 0.05, "nmaxold": 8, "weight": 100},
                   "convergence": {"energy": 1e-5, "density": 1e-5, "eigenstates": 1e-5},
                   "GPAW_version": "26.7.0", "ASE_version": "3.29.0", "mpi_ranks": world.size},
        "elapsed_s": time.time() - started, "dataset_role": job["role"],
        "training_integration": "not performed",
    }
    (scratch / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    progress_path.write_text(json.dumps({"label": LABEL, "status": "complete", "iteration": iterations,
        "elapsed_s": summary["elapsed_s"], "scf_converged": True,
        "source_sha256": job["input_sha256"], "mpi_ranks": world.size}, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
world.barrier()
