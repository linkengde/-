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
from entry import load
# Reject disabled entries before MPI setup; initialize MPI before ASE input reads.
preflight=json.loads((ROOT/'input_manifest.json').read_text())
assert preflight['launch_enabled'] and next(r for r in preflight['records'] if r['label']==LABEL)['launch_enabled'], 'Launch disabled'
from gpaw.mpi import world
manifest, record = load(LABEL, launch=True)
assert SOURCE == (ROOT / record['input']).resolve(), 'Wrong input path'
if OUT.exists(): raise SystemExit('Existing run: preserve it; no duplicate launch')

import gpaw_data
os.environ["GPAW_SETUP_PATH"] = str(gpaw_data.datapath())
import _gpaw
assert hasattr(_gpaw, 'scalapack_diagonalize_dc'), 'ScaLAPACK build required'
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
if atoms.pbc.tolist() != record['pbc']:
    raise SystemExit("Input PBC mismatch")
if any(k in atoms.info for k in ("REF_energy", "PW_PBE_energy_eV", "energy")):
    raise SystemExit("Unexpected energy label in blind input")
if any(k in atoms.arrays for k in ("REF_forces", "PW_PBE_forces", "forces")):
    raise SystemExit("Unexpected force label in blind input")
started = time.time()
progress_path = OUT / "progress.json"
state_path = OUT / "state.gpw"  # retained in local workspace run; never archive the large checkpoint
if world.rank == 0:
    progress_path.write_text(json.dumps({
        "label": LABEL, "status": "running", "iteration": 0,
        "source_sha256": record["input_sha256"], "mpi_ranks": world.size,
        "started_utc_epoch_s": started,
    }, indent=2) + "\n")
calc = GPAW(
    mode=PW(record['method']['cutoff_eV']), xc=record['method']['xc'], kpts=tuple(record['method']['kpts']), poissonsolver=record['poissonsolver'],
    occupations=FermiDirac(record['method']['smearing_eV']), mixer=Mixer(0.05, 8, 100),
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
    if shutil.disk_usage(OUT).free < manifest["disk_budget"]["minimum_before_checkpoint_bytes"]:
        raise SystemExit("Disk reserve depleted during SCF; preserve same run, no new checkpoint write")
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
atoms.info["config_type"] = LABEL
atoms.info["PW_PBE_energy_eV"] = energy
atoms.info["source_method"] = "GPAW 26.7.0 PW-PBE bulk pure-phase numerical control"
atoms.arrays["PW_PBE_forces"] = forces.copy()

if world.rank == 0:
    header = 'Lattice="' + ' '.join(format(v, '.17g') for v in atoms.cell.array.ravel()) + '" Properties=species:S:1:pos:R:3:PW_PBE_forces:R:3:local_new_id:I:1 pbc="' + ' '.join('T' if p else 'F' for p in atoms.pbc) + '" PW_PBE_energy_eV=' + format(energy, '.17g') + ' config_type=' + LABEL
    (OUT / f"{LABEL}_PW_PBE.extxyz").write_text(str(len(atoms)) + '\n' + header + '\n' + ''.join(symbol + ' ' + ' '.join(format(v, '.17g') for v in list(xyz) + list(force)) + f' {int(atom_id)}\n' for symbol, xyz, force, atom_id in zip(atoms.get_chemical_symbols(), atoms.positions, forces, atoms.arrays["local_new_id"])))
    summary = {
        "label": LABEL, "source": str(SOURCE), "source_sha256": record["input_sha256"],
        "atoms": len(atoms), "formula": atoms.get_chemical_formula(),
        "energy_eV_cell": energy, "free_energy_eV_cell": free_energy,
        "energy_convention": "GPAW native extrapolated energy, force_consistent=False; free_energy recorded separately", "scf_converged": converged,
        "scf_iterations": iterations,
        "fmax_eV_A": float(np.linalg.norm(forces, axis=1).max()),
        "rms_force_eV_A": float(np.sqrt(np.mean(np.sum(forces**2, axis=1)))),
        "method": record['method'],
        "mpi_ranks": world.size, "elapsed_s": time.time() - started,
        "dataset_role": record["role"],
        "note": "CIF-derived bulk pure-phase numerical control; deliberate near-parent geometry, no independent-validation claim.",
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
