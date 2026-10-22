#!/usr/bin/env python3
"""Four-rank fixed-geometry GPAW PW-PBE label for a registered V19 dev row."""
import json
import os
import shutil
import sys
import time
from pathlib import Path

import numpy as np

if len(sys.argv) != 4:
    raise SystemExit("usage: run_pw_reference.py LABEL INPUT.extxyz OUTPUT_DIR")
LABEL = sys.argv[1]
SOURCE = Path(sys.argv[2]).resolve()
OUT = Path(sys.argv[3]).resolve()
ROOT = Path(__file__).resolve().parent
from entry import load

manifest, record = load(LABEL, launch=True)
assert SOURCE == (ROOT.parents[3] / record["input"]).resolve(), "Wrong registered input"
assert not OUT.exists(), "Existing run directory: preserve and inspect it"

from gpaw.mpi import world

if world.size != 4:
    raise SystemExit(f"Expected four MPI ranks; observed {world.size}")

import gpaw_data
os.environ["GPAW_SETUP_PATH"] = str(gpaw_data.datapath())
import _gpaw

assert hasattr(_gpaw, "scalapack_diagonalize_dc"), "GPAW build lacks ScaLAPACK"
from ase.io import read
from gpaw import GPAW, PW
from gpaw.mixer import Mixer
from gpaw.occupations import FermiDirac

if world.rank == 0:
    OUT.mkdir(parents=True, exist_ok=False)
world.barrier()
atoms = read(SOURCE)
assert atoms.pbc.tolist() == record["pbc"]
assert len(atoms) == record["atoms"] and atoms.get_chemical_formula() == record["formula"]
assert "REF_energy" not in atoms.info and "energy" not in atoms.info
assert "REF_forces" not in atoms.arrays and "forces" not in atoms.arrays

started = time.time()
progress_path = OUT / "progress.json"
checkpoint = OUT / "state.gpw"
if world.rank == 0:
    progress_path.write_text(json.dumps({"label": LABEL, "status": "running", "iteration": 0, "input_sha256": record["input_sha256"], "mpi_ranks": 4, "started_epoch_s": started}, indent=2) + "\n")

calc = GPAW(
    mode=PW(record["method"]["cutoff_eV"]),
    xc=record["method"]["xc"],
    kpts=tuple(record["method"]["kpts"]),
    poissonsolver=record["poissonsolver"],
    occupations=FermiDirac(record["method"]["smearing_eV"]),
    mixer=Mixer(0.05, 8, 100),
    parallel={"sl_auto": True},
    convergence={"energy": 1e-5, "density": 1e-5, "eigenstates": 1e-5},
    maxiter=160,
    txt=str(OUT / "gpaw.log"),
)
atoms.calc = calc
converged = {"value": False}


def scf_checkpoint(context):
    iteration = int(context.niter)
    if world.rank == 0 and iteration % 10 == 0:
        progress_path.write_text(json.dumps({"label": LABEL, "status": "running", "iteration": iteration, "input_sha256": record["input_sha256"], "mpi_ranks": 4, "elapsed_s": time.time() - started}, indent=2) + "\n")
    if shutil.disk_usage(OUT).free < manifest["minimum_free_disk_during_scf_bytes"]:
        raise SystemExit("Disk reserve low during SCF; preserve this run and checkpoint")
    if iteration > 0 and iteration % 20 == 0:
        calc.write(str(checkpoint), mode="all")


calc.hooks["scf_step"] = scf_checkpoint
calc.hooks["converged"] = lambda: converged.__setitem__("value", True)
energy = float(atoms.get_potential_energy(force_consistent=False))
free_energy = float(calc.results["free_energy"])
forces = np.asarray(atoms.get_forces(), dtype=float)
iterations = int(calc.get_number_of_iterations())
if not converged["value"]:
    raise SystemExit(f"SCF not converged after {iterations} iterations; preserve outputs")
if not np.isfinite(energy) or not np.isfinite(free_energy) or not np.isfinite(forces).all():
    raise SystemExit("Non-finite DFT result; preserve outputs")

if world.rank == 0:
    atoms.info["config_type"] = LABEL
    atoms.info["PW_PBE_energy_eV"] = energy
    atoms.arrays["PW_PBE_forces"] = forces.copy()
    fields = "Properties=species:S:1:pos:R:3:PW_PBE_forces:R:3:lammps_id:I:1:lammps_type:I:1:central_pair:I:1"
    cell = " ".join(format(float(v), ".17g") for v in atoms.cell.array.ravel())
    pbc = " ".join("T" if value else "F" for value in atoms.pbc)
    header = f'Lattice="{cell}" {fields} pbc="{pbc}" PW_PBE_energy_eV={energy:.17g} config_type={LABEL}'
    rows = [f"{len(atoms)}\n{header}\n"]
    for symbol, position, force, atom_id, atom_type, marker in zip(
        atoms.get_chemical_symbols(), atoms.positions, forces, atoms.arrays["lammps_id"], atoms.arrays["lammps_type"], atoms.arrays["central_pair"]
    ):
        vals = " ".join(format(float(v), ".17g") for v in list(position) + list(force))
        rows.append(f"{symbol} {vals} {int(atom_id)} {int(atom_type)} {int(marker)}\n")
    (OUT / f"{LABEL}_PW_PBE.extxyz").write_text("".join(rows))
    summary = {
        "label": LABEL,
        "input": str(SOURCE.relative_to(ROOT.parents[3])),
        "input_sha256": record["input_sha256"],
        "atoms": len(atoms),
        "formula": atoms.get_chemical_formula(),
        "dataset_role": record["role"],
        "energy_eV_cell": energy,
        "free_energy_eV_cell": free_energy,
        "energy_convention": "GPAW native extrapolated energy, force_consistent=False; free energy separately recorded",
        "scf_converged": True,
        "scf_iterations": iterations,
        "fmax_eV_A": float(np.linalg.norm(forces, axis=1).max()),
        "rms_force_eV_A": float(np.sqrt(np.mean(np.sum(forces ** 2, axis=1)))),
        "method": record["method"],
        "pbc": record["pbc"],
        "poissonsolver": record["poissonsolver"],
        "mpi_ranks": 4,
        "elapsed_s": time.time() - started,
        "frozen_geometry": True,
        "source_CIF_sha256": record["source_CIF_sha256"],
        "geometry_sha256": record["geometry_sha256"],
        "limitations": ["The provisional slab recipe has unresolved k-point/smearing/vacuum convergence; this is a consistent fixed-geometry development label, not a production accuracy claim."],
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    progress_path.write_text(json.dumps({"label": LABEL, "status": "complete", "iteration": iterations, "elapsed_s": summary["elapsed_s"], "scf_converged": True, "input_sha256": record["input_sha256"], "mpi_ranks": 4}, indent=2) + "\n")
    archive = ROOT.parents[3] / record["result_directory"] / "calculations" / LABEL
    archive.parent.mkdir(parents=True, exist_ok=True)
    if archive.exists():
        raise SystemExit(f"Refusing to replace archive {archive}")
    staging = archive.parent / f".{LABEL}.staging-{os.getpid()}"
    staging.mkdir()
    for filename in (f"{LABEL}_PW_PBE.extxyz", "summary.json", "gpaw.log", "progress.json"):
        shutil.copy2(OUT / filename, staging / filename)
    staging.rename(archive)
    print(json.dumps(summary, indent=2), flush=True)
world.barrier()
