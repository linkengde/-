#!/usr/bin/env python3
"""Run a fixed-geometry GPAW PW-PBE label for a wider-span v17 candidate."""
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import sys
import time

import numpy as np
from ase.io import read, write

if len(sys.argv) != 4:
    raise SystemExit("usage: run_pw_reference.py LABEL INPUT.extxyz OUTPUT_DIR")
LABEL = sys.argv[1]
SOURCE = Path(sys.argv[2]).resolve()
OUT = Path(sys.argv[3]).resolve()
ROOT = Path(__file__).resolve().parent
manifest = json.loads((ROOT / "input_manifest.json").read_text())
record = next((r for r in manifest["records"] if r["label"] == LABEL), None)
if record is None or sha256(SOURCE.read_bytes()).hexdigest() != record["input_sha256"]:
    raise SystemExit(f"Input hash mismatch for {LABEL}")

import gpaw_data
os.environ["GPAW_SETUP_PATH"] = str(gpaw_data.datapath())
from gpaw import GPAW, PW
from gpaw.mixer import Mixer
from gpaw.mpi import world
from gpaw.occupations import FermiDirac

if world.size != 4:
    raise SystemExit(f"Expected 4 MPI ranks, got {world.size}")
if world.rank == 0:
    OUT.mkdir(parents=True, exist_ok=False)
world.barrier()

atoms = read(SOURCE)
if not np.all(atoms.pbc):
    raise SystemExit("Input must be periodic in all directions")
if any(k in atoms.info for k in ("REF_energy", "PW_PBE_energy_eV", "energy")):
    raise SystemExit("Unexpected energy label in blind input")
if any(k in atoms.arrays for k in ("REF_forces", "PW_PBE_forces", "forces")):
    raise SystemExit("Unexpected force label in blind input")
started = time.time()
progress_path = OUT / "progress.json"
state_path = OUT / "state.gpw"  # kept in /tmp; do not archive this large checkpoint
if world.rank == 0:
    progress_path.write_text(json.dumps({
        "label": LABEL, "status": "running", "iteration": 0,
        "source_sha256": record["input_sha256"], "mpi_ranks": world.size,
        "started_utc_epoch_s": started,
    }, indent=2) + "\n")
calc = GPAW(
    mode=PW(500), xc="PBE", kpts=(1, 1, 1),
    occupations=FermiDirac(0.1), mixer=Mixer(0.05, 8, 100),
    parallel={"sl_auto": True},
    convergence={"energy": 1e-5, "density": 1e-5, "eigenstates": 1e-5},
    maxiter=160, txt=str(OUT / "gpaw.log"),
)
atoms.calc = calc

def checkpoint(ctx):
    niter = int(ctx.niter)
    if niter % 10 == 0 and world.rank == 0:
        progress_path.write_text(json.dumps({
            "label": LABEL, "status": "running", "iteration": niter,
            "source_sha256": record["input_sha256"], "mpi_ranks": world.size,
            "elapsed_s": time.time() - started,
        }, indent=2) + "\n")
    if niter % 20 == 0:
        calc.write(str(state_path), mode="all")

calc.hooks["scf_step"] = checkpoint
converged_status = {"value": False}
calc.hooks["converged"] = lambda: converged_status.__setitem__("value", True)
energy = float(atoms.get_potential_energy(force_consistent=False))
free_energy = float(calc.results["free_energy"]) if "free_energy" in calc.results else None
forces = np.asarray(atoms.get_forces(), dtype=float)
if not np.isfinite(energy) or not np.isfinite(forces).all() or (free_energy is not None and not np.isfinite(free_energy)):
    raise SystemExit("Nonfinite DFT label")
converged = bool(converged_status["value"])
iterations = int(calc.get_number_of_iterations())
central = np.asarray(atoms.arrays["central_pair"], dtype=bool)
symbols = np.asarray(atoms.get_chemical_symbols())
ia = int(np.flatnonzero(central & (symbols == "Ag"))[0])
contact = record["central_pair"]["elements"].split("-")[1]
ib = int(np.flatnonzero(central & (symbols == contact))[0])
vector = atoms.positions[ib] - atoms.positions[ia]
distance = float(np.linalg.norm(vector))
separating_force = float(np.dot(forces[ib] - forces[ia], vector / distance))
atoms.info["config_type"] = LABEL
atoms.info["PW_PBE_energy_eV"] = energy
atoms.info["source_method"] = "GPAW 26.7.0 PW-PBE 500 eV Gamma single point"
atoms.arrays["PW_PBE_forces"] = forces.copy()
calc.write(str(state_path), mode="all")

if world.rank == 0:
    write(OUT / f"{LABEL}_PW_PBE.extxyz", atoms, format="extxyz")
    summary = {
        "label": LABEL, "source": str(SOURCE), "source_sha256": record["input_sha256"],
        "atoms": len(atoms), "formula": atoms.get_chemical_formula(),
        "central_pair": {
            "elements": record["central_pair"]["elements"],
            "Ag_id": int(atoms.arrays["lammps_id"][ia]),
            "contact_id": int(atoms.arrays["lammps_id"][ib]),
            "distance_A": distance, "separating_force_eV_A": separating_force,
        },
        "energy_eV_cell": energy, "free_energy_eV_cell": free_energy,
        "energy_convention": "GPAW native extrapolated energy, force_consistent=False; free_energy recorded separately", "scf_converged": converged,
        "scf_iterations": iterations,
        "fmax_eV_A": float(np.linalg.norm(forces, axis=1).max()),
        "rms_force_eV_A": float(np.sqrt(np.mean(np.sum(forces**2, axis=1)))),
        "method": {
            "xc": "PBE", "basis": "plane wave", "cutoff_eV": 500,
            "kpts": [1, 1, 1], "smearing_eV": 0.1,
            "parallel": {"sl_auto": True},
            "scf_thresholds": "energy/density/eigenstates 1e-5",
            "mixer": "Pulay beta=0.05, nmaxold=8, weight=100",
            "GPAW_version": __import__("gpaw").__version__,
        },
        "mpi_ranks": world.size, "elapsed_s": time.time() - started,
        "dataset_role": record["role"],
        "note": "V17 paired diagnostic/training acquisition; deliberate near-parent geometry, no independent-validation claim.",
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    progress_path.write_text(json.dumps({
        "label": LABEL, "status": "complete" if converged else "not_converged",
        "iteration": iterations, "elapsed_s": summary["elapsed_s"],
        "scf_converged": converged, "source_sha256": record["input_sha256"],
        "mpi_ranks": world.size,
    }, indent=2) + "\n")
    if converged:
        archive = ROOT / "calculations" / LABEL
        if archive.exists():
            raise SystemExit(f"Refusing to overwrite existing archive: {archive}")
        archive.parent.mkdir(parents=True, exist_ok=True)
        staging = archive.parent / f".{LABEL}.staging-{os.getpid()}"
        staging.mkdir()
        for name in (f"{LABEL}_PW_PBE.extxyz", "summary.json", "gpaw.log", "progress.json"):
            shutil.copy2(OUT / name, staging / name)
        staging.rename(archive)
    print(json.dumps(summary, indent=2), flush=True)
world.barrier()
if not converged:
    raise SystemExit(f"SCF did not converge for {LABEL}")
