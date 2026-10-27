#!/usr/bin/env python3
"""Run one registered fixed-geometry four-rank GPAW PBE/PW500 label."""
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
    raise SystemExit("usage: run_dft.py LABEL SCRATCH_DIR")
LABEL, SCRATCH_ARG = sys.argv[1:]
manifest = json.loads((ROOT / "input_manifest.json").read_text())
record = next(item for item in manifest["records"] if item["label"] == LABEL)
assert manifest["launch_enabled"] and record["launch_enabled"]
assert record["status"] == "queued"
assert world.size == 4, f"expected four MPI ranks, got {world.size}"

source = ROOT / record["input"]
if hashlib.sha256(source.read_bytes()).hexdigest() != record["input_sha256"]:
    raise SystemExit("input SHA256 changed after registration")
scratch = Path(SCRATCH_ARG).resolve()
if world.rank == 0:
    scratch.mkdir(parents=True, exist_ok=False)
world.barrier()

import gpaw_data

os.environ["GPAW_SETUP_PATH"] = str(gpaw_data.datapath())
from ase.io import read
from gpaw import GPAW, PW
from gpaw.mixer import Mixer
from gpaw.occupations import FermiDirac

atoms = read(source)
method = record["method"]
assert len(atoms) == record["atom_count"]
assert atoms.get_chemical_formula() == record["formula"]
assert atoms.pbc.tolist() == record["method"]["pbc"]
assert not any(key in atoms.info for key in ("energy", "REF_energy", "PW_PBE_energy_eV"))
assert not any(key in atoms.arrays for key in ("forces", "REF_forces", "PW_PBE_forces"))

started = time.time()
progress_path = scratch / "progress.json"
state_path = scratch / "state.gpw"
if world.rank == 0:
    progress_path.write_text(json.dumps({
        "label": LABEL,
        "status": "running",
        "iteration": 0,
        "source_sha256": record["input_sha256"],
        "mpi_ranks": world.size,
        "started_epoch_s": started,
    }, indent=2) + "\n")

calc = GPAW(
    mode=PW(method["cutoff_eV"]),
    xc=method["xc"],
    kpts=tuple(method["kpts"]),
    occupations=FermiDirac(method["smearing_eV"]),
    mixer=Mixer(method["mixer"]["beta"], method["mixer"]["nmaxold"], method["mixer"]["weight"]),
    parallel={"sl_auto": True},
    convergence=method["scf_thresholds"],
    maxiter=method["maxiter"],
    txt=str(scratch / "gpaw.log"),
)
atoms.calc = calc
converged = {"value": False}


def scf_step(context) -> None:
    iteration = int(context.niter)
    if world.rank == 0 and iteration % 10 == 0:
        progress_path.write_text(json.dumps({
            "label": LABEL,
            "status": "running",
            "iteration": iteration,
            "elapsed_s": time.time() - started,
            "source_sha256": record["input_sha256"],
            "mpi_ranks": world.size,
        }, indent=2) + "\n")
    if iteration % 20 == 0:
        calc.write(str(state_path), mode="all")


calc.hooks["scf_step"] = scf_step
calc.hooks["converged"] = lambda: converged.__setitem__("value", True)
native_energy = float(atoms.get_potential_energy(force_consistent=False))
free_energy = float(calc.results["free_energy"])
forces = np.asarray(atoms.get_forces(), dtype=float)
iterations = int(calc.get_number_of_iterations())
if not converged["value"]:
    raise SystemExit(f"SCF failed to converge within {iterations} iterations")
if not np.isfinite(native_energy) or not np.isfinite(free_energy) or not np.isfinite(forces).all():
    raise SystemExit("nonfinite DFT result")

if world.rank == 0:
    flat_cell = " ".join(format(float(value), ".17g") for value in atoms.cell.array.ravel())
    header = (
        f'Lattice="{flat_cell}" '
        'Properties=species:S:1:pos:R:3:PW_PBE_forces:R:3:local_new_id:I:1 '
        'pbc="T T T" '
        f'PW_PBE_energy_eV={native_energy:.17g} '
        f'PW_PBE_free_energy_eV={free_energy:.17g} '
        f'config_type={LABEL}'
    )
    lines = [str(len(atoms)), header]
    for symbol, position, force, local_id in zip(
        atoms.get_chemical_symbols(), atoms.positions, forces, atoms.arrays["local_new_id"]
    ):
        values = list(position) + list(force)
        lines.append(symbol + " " + " ".join(format(float(value), ".17g") for value in values) + f" {int(local_id)}")
    (scratch / f"{LABEL}_PW_PBE.extxyz").write_text("\n".join(lines) + "\n")
    summary = {
        "label": LABEL,
        "input": str(source.relative_to(ROOT)),
        "input_sha256": record["input_sha256"],
        "source_family": "Jeitschko_Nowotny_1967_single_crystal",
        "source_doi": "10.1007/BF00899949",
        "source_cif_sha256": record["source_cif_sha256"],
        "atom_count": len(atoms),
        "formula": atoms.get_chemical_formula(),
        "cell_angstrom": atoms.cell.array.tolist(),
        "pbc": atoms.pbc.tolist(),
        "native_extrapolated_energy_eV_cell": native_energy,
        "free_energy_eV_cell": free_energy,
        "forces_method": "GPAW ASE forces at Fermi-Dirac width 0.1 eV; preserve free_energy separately",
        "energy_convention": "native extrapolated energy (force_consistent=False); free energy separately recorded",
        "scf_converged": True,
        "scf_iterations": iterations,
        "max_force_eV_A": float(np.linalg.norm(forces, axis=1).max()),
        "rms_force_vector_eV_A": float(np.sqrt(np.mean(np.sum(forces**2, axis=1)))),
        "method": method,
        "mpi_ranks": world.size,
        "elapsed_s": time.time() - started,
        "dataset_role": record["role"],
        "training_integration": "not performed; pending separate audit and versioned dataset review",
    }
    (scratch / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    progress_path.write_text(json.dumps({
        "label": LABEL,
        "status": "complete",
        "iteration": iterations,
        "elapsed_s": summary["elapsed_s"],
        "scf_converged": True,
        "source_sha256": record["input_sha256"],
        "mpi_ranks": world.size,
    }, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
world.barrier()
