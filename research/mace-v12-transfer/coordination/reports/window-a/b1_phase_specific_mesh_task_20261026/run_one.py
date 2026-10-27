#!/usr/bin/env python3
"""Run one registered fixed-geometry GPAW mesh control under 4 MPI ranks."""
import argparse
import hashlib
import json
import shutil
import time
from pathlib import Path

import numpy as np
from ase.io import read, write

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[5]
BASE = ROOT / "research/mace-v12-transfer"
TASK = json.loads((HERE / "task_manifest.json").read_text())
parser = argparse.ArgumentParser()
parser.add_argument("label")
parser.add_argument("--conditional", action="store_true")
args = parser.parse_args()

records = TASK["records"] if not args.conditional else TASK["conditional_records"]
matches = [record for record in records if record["label"] == args.label]
if len(matches) != 1:
    raise SystemExit(f"Unregistered label or wrong conditional mode: {args.label}")
record = matches[0]
if not TASK["global_launch_enabled"] or not TASK["global_dft_authorized"]:
    raise SystemExit("Task launch is disabled")
if args.conditional:
    compare_path = HERE / "comparisons" / f"{record['pair_id']}.json"
    if record["launch_enabled"] is not False or not compare_path.is_file():
        raise SystemExit("Conditional liquid endpoint lacks its trigger comparison")
    prior = json.loads(compare_path.read_text())
    if prior["pass"] is not False or prior["fine_label"] != record["compare_to_label"]:
        raise SystemExit("Conditional endpoint trigger does not match a failed Gamma-to-2x2x2 pair")
elif not record["launch_enabled"] or not record["dft_authorized"]:
    raise SystemExit("Mandatory label is not launch-authorized")

source = (ROOT / record["input"]).resolve()
if not source.is_relative_to(ROOT):
    raise SystemExit("Input path escapes repository root")
source_sha = hashlib.sha256(source.read_bytes()).hexdigest()
if source_sha != record["input_sha256"]:
    raise SystemExit("Input SHA256 changed after registration")
OUT = HERE / TASK["outputs_directory"] / record["label"]
if OUT.exists():
    raise SystemExit(f"Existing run directory; preserve and review: {OUT}")

import gpaw_data
import os
os.environ["GPAW_SETUP_PATH"] = str(gpaw_data.datapath())
from gpaw.mpi import world
if world.size != TASK["resource_cap"]["mpi_ranks_per_endpoint"]:
    raise SystemExit(f"Expected four MPI ranks, got {world.size}")
import _gpaw
if not hasattr(_gpaw, "scalapack_diagonalize_dc"):
    raise SystemExit("This B1 task requires the verified ScaLAPACK GPAW build")
from gpaw import GPAW, PW
from gpaw.mixer import Mixer
from gpaw.occupations import FermiDirac

setup_dir = Path(gpaw_data.datapath())
for filename, expected in TASK["method"]["PAW_setup_sha256"].items():
    actual = hashlib.sha256((setup_dir / filename).read_bytes()).hexdigest()
    if actual != expected:
        raise SystemExit(f"PAW setup hash mismatch: {filename}")

if world.rank == 0:
    OUT.mkdir(parents=True, exist_ok=False)
world.barrier()
atoms = read(source)
if atoms.pbc.tolist() != record["pbc"] or len(atoms) != record["natoms"]:
    raise SystemExit("Registered input metadata changed")
if any(key in atoms.info for key in ("REF_energy", "PW_PBE_energy_eV", "free_energy_eV", "energy")):
    raise SystemExit("Input unexpectedly contains an energy label")
if any(key in atoms.arrays for key in ("REF_forces", "PW_PBE_forces", "forces")):
    raise SystemExit("Input unexpectedly contains force labels")

started = time.time()
progress_path = OUT / "progress.json"
state_path = OUT / "state.gpw"
if world.rank == 0:
    progress_path.write_text(json.dumps({
        "label": record["label"], "status": "running", "iteration": 0,
        "source_sha256": source_sha, "mpi_ranks": world.size,
        "started_utc_epoch_s": started,
    }, indent=2) + "\n")
method = record["method"]
initial_wavefunctions = record.get("initial_wavefunctions", "lcao")
if initial_wavefunctions not in {"lcao", "random"}:
    raise SystemExit(f"Unsupported initial wavefunction mode: {initial_wavefunctions}")
calc = GPAW(
    mode=PW(method["cutoff_eV"]), xc=method["xc"], kpts=tuple(method["kpts"]),
    occupations=FermiDirac(method["smearing_eV"]), mixer=Mixer(0.05, 8, 100),
    parallel={"sl_auto": True}, random=(initial_wavefunctions == "random"),
    convergence={"energy": 1e-5, "density": 1e-5, "eigenstates": 1e-5},
    maxiter=method["maxiter"], txt=str(OUT / "gpaw.log"),
)
atoms.calc = calc

def checkpoint(ctx):
    iteration = int(ctx.niter)
    if world.rank == 0 and iteration % 10 == 0:
        progress_path.write_text(json.dumps({
            "label": record["label"], "status": "running", "iteration": iteration,
            "source_sha256": source_sha, "mpi_ranks": world.size,
            "elapsed_s": time.time() - started,
        }, indent=2) + "\n")
    if shutil.disk_usage(OUT).free < 3_000_000_000:
        raise SystemExit("Less than 3 GB free during SCF; stop and preserve this run")
    if iteration % 20 == 0:
        calc.write(str(state_path), mode="all")

calc.hooks["scf_step"] = checkpoint
converged_flag = {"value": False}
calc.hooks["converged"] = lambda: converged_flag.__setitem__("value", True)
native_energy = float(atoms.get_potential_energy(force_consistent=False))
free_energy = float(calc.results["free_energy"])
forces = np.asarray(atoms.get_forces(), dtype=float)
if not np.isfinite(native_energy) or not np.isfinite(free_energy) or not np.isfinite(forces).all():
    raise SystemExit("Nonfinite result")
iterations = int(calc.get_number_of_iterations())
converged = bool(converged_flag["value"])
if world.rank == 0:
    result = atoms.copy()
    result.calc = None
    result.info["config_type"] = record["label"]
    result.info["source_sha256"] = source_sha
    result.info["PW_PBE_energy_eV"] = native_energy
    result.info["free_energy_eV"] = free_energy
    result.info["dataset_role"] = "NUMERICAL_CONTROL_ONLY"
    result.arrays["PW_PBE_forces"] = forces.copy()
    result.arrays["free_energy_forces"] = forces.copy()
    write(OUT / "result.extxyz", result, format="extxyz")
    summary = {
        "task_id": TASK["task_id"], "label": record["label"], "phase": record["phase"],
        "pair_id": record["pair_id"], "pair_member": record["pair_member"],
        "source": str(source.relative_to(ROOT)), "source_sha256": source_sha,
        "atoms": len(atoms), "formula": atoms.get_chemical_formula(),
        "native_energy_eV_cell": native_energy, "free_energy_eV_cell": free_energy,
        "scf_converged": converged, "scf_iterations": iterations,
        "fmax_eV_A": float(np.linalg.norm(forces, axis=1).max()),
        "rms_force_eV_A": float(np.sqrt(np.mean(np.sum(forces * forces, axis=1)))),
        "method": method, "initial_wavefunctions": initial_wavefunctions,
        "mpi_ranks": world.size, "elapsed_s": time.time() - started,
        "dataset_role": "NUMERICAL_CONTROL_ONLY",
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    progress_path.write_text(json.dumps({
        "label": record["label"], "status": "complete" if converged else "not_converged",
        "iteration": iterations, "elapsed_s": summary["elapsed_s"],
        "scf_converged": converged, "source_sha256": source_sha, "mpi_ranks": world.size,
    }, indent=2) + "\n")
world.barrier()
if not converged:
    raise SystemExit(f"SCF failed to converge for {record['label']}")
