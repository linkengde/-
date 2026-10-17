#!/usr/bin/env python3
"""Run one C-reviewed, fixed-geometry Ti3SiC2 surface single point under MPI."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import time
from pathlib import Path

import ase
import gpaw as gpaw_module
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write
from gpaw import FermiDirac, GPAW, PW
from gpaw.dft import MSR1
from gpaw.mpi import world


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def geometry_sha(atoms) -> str:
    h = hashlib.sha256()
    h.update(" ".join(atoms.get_chemical_symbols()).encode())
    h.update(np.asarray(atoms.cell, dtype="<f8").tobytes(order="C"))
    h.update(np.asarray(np.round(atoms.positions, 8), dtype="<f8").tobytes(order="C"))
    h.update(np.asarray(atoms.pbc, dtype=np.uint8).tobytes())
    h.update(np.asarray(atoms.arrays["source_atom_id"], dtype="<i8").tobytes())
    return h.hexdigest()


parser = argparse.ArgumentParser()
parser.add_argument("--task-id", required=True)
parser.add_argument("--input", required=True)
parser.add_argument("--run-dir", required=True)
parser.add_argument("--config", required=True)
parser.add_argument("--config-sha256", required=True)
parser.add_argument("--execution-lock", required=True)
parser.add_argument("--execution-lock-sha256", required=True)
args = parser.parse_args()

input_path = Path(args.input).resolve()
run_dir = Path(args.run_dir).resolve()
config_path = Path(args.config).resolve()
lock_path = Path(args.execution_lock).resolve()
if not run_dir.is_dir():
    raise RuntimeError(f"Run directory must be prepared by the supervisor: {run_dir}")
if sha256(config_path) != args.config_sha256:
    raise RuntimeError("Scientific configuration SHA mismatch")
if sha256(lock_path) != args.execution_lock_sha256:
    raise RuntimeError("Execution lock SHA mismatch")

config = json.loads(config_path.read_text())
lock = json.loads(lock_path.read_text())
assert lock.get("execution_authorized") is True
assert lock.get("checkpoint_policy", {}).get("checkpoint_enabled") is False
assert lock.get("role") == "TRAIN_FORCE_LABEL_CANDIDATE"
assert ase.__version__ == config["software"]["ASE_version"]
assert gpaw_module.__version__ == config["software"]["GPAW_version"]
assert int(os.environ.get("OMP_NUM_THREADS", "0")) == config["parallel"]["OMP_NUM_THREADS"]
assert os.environ.get("GPAW_MPI_BACKEND") == config["parallel"]["GPAW_MPI_BACKEND"]
assert world.size == config["parallel"]["MPI_ranks"]
assert sha256(Path(__file__).resolve()) == lock["runner_sha256"]

task = next(t for t in config["tasks"] if t["structure_id"] == args.task_id)
assert sha256(input_path) == task["input_sha256"], "Input SHA mismatch"
assert input_path.name == "input.extxyz"
atoms = read(input_path, format="extxyz")
symbols = atoms.get_chemical_symbols()
assert len(atoms) == task["atom_count"] == 100
assert {s: symbols.count(s) for s in sorted(set(symbols))} == task["composition"]
assert atoms.pbc.tolist() == [True, True, True]
assert np.array_equal(atoms.arrays.get("source_atom_id"), np.arange(1, 101))
assert "energy" not in atoms.info and "free_energy" not in atoms.info
assert "forces" not in atoms.arrays and "stress" not in atoms.info
before_geometry_sha = geometry_sha(atoms)

method = config["method"]
mix = method["mixer"]
t0 = time.monotonic()
calc = GPAW(
    mode=PW(method["cutoff_eV"]),
    xc=method["xc"],
    kpts=tuple(method["kpoints_size"]),
    occupations=FermiDirac(method["fermi_dirac_width_eV"]),
    spinpol=method["spin_polarized"],
    convergence=method["scf_convergence"],
    maxiter=method["maxiter"],
    mixer=MSR1(
        nmaxold=mix["nmaxold"], beta=mix["beta"], reg=mix["reg"],
        gb_scale=mix["gb_scale"], max_A=mix["max_A"],
        trust_scale=mix["trust_scale"], soft_lim=mix["soft_lim"],
        hard_lim=mix["hard_lim"], weight=mix["weight"], sigma=mix["sigma"],
    ),
    parallel={"sl_auto": method["sl_auto"]},
    txt=str(run_dir / "gpaw.log"),
)
atoms.calc = calc

energy = float(atoms.get_potential_energy(force_consistent=False))
free_energy = float(atoms.get_potential_energy(force_consistent=True))
forces = np.asarray(atoms.get_forces(), dtype=float)
stress = np.asarray(atoms.get_stress(voigt=False), dtype=float)
log_text = (run_dir / "gpaw.log").read_text(errors="replace")
matches = re.findall(r"Converged in (\d+) steps", log_text)
if not matches:
    raise RuntimeError("No converged-SCF marker; final labels are not written")
iterations = int(matches[-1])
reported_iterations = int(calc.get_number_of_iterations())
if reported_iterations != iterations:
    raise RuntimeError(f"SCF iteration mismatch: log={iterations}, GPAW={reported_iterations}")
if forces.shape != (100, 3) or stress.shape != (3, 3):
    raise RuntimeError(f"Unexpected result dimensions: forces={forces.shape}, stress={stress.shape}")
if not (np.isfinite(energy) and np.isfinite(free_energy)
        and np.all(np.isfinite(forces)) and np.all(np.isfinite(stress))):
    raise RuntimeError("Non-finite DFT result; refusing to write labels")
if geometry_sha(atoms) != before_geometry_sha:
    raise RuntimeError("Coordinates/cell/PBC/atom order changed during fixed-geometry calculation")

if world.rank == 0:
    labelled = atoms.copy()
    labelled.calc = SinglePointCalculator(
        labelled, energy=energy, free_energy=free_energy, forces=forces, stress=stress
    )
    output_path = run_dir / "labelled.extxyz"
    write(output_path, labelled, format="extxyz")
    result = {
        "structure_id": args.task_id,
        "status": "SUCCESS",
        "scf_converged": True,
        "scf_iterations": iterations,
        "energy_eV": energy,
        "force_consistent_free_energy_eV": free_energy,
        "forces_eV_per_angstrom": forces.tolist(),
        "stress_eV_per_angstrom3": stress.tolist(),
        "stress_GPa": (stress * 160.21766208).tolist(),
        "max_force_eV_per_angstrom": float(np.linalg.norm(forces, axis=1).max()),
        "atom_count": len(atoms),
        "composition": task["composition"],
        "input_sha256": task["input_sha256"],
        "geometry_sha256": task["geometry_sha256"],
        "config_sha256": args.config_sha256,
        "execution_lock_sha256": args.execution_lock_sha256,
        "runner_sha256": lock["runner_sha256"],
        "kpoint_type": method["kpoint_type"],
        "kpoints": method["kpoints_size"],
        "GPAW_version": gpaw_module.__version__,
        "ASE_version": ase.__version__,
        "elapsed_seconds": time.monotonic() - t0,
        "role": lock["role"],
        "units": {"energy": "eV", "forces": "eV/Angstrom", "stress": "eV/Angstrom^3 and GPa"},
    }
    (run_dir / "RESULT.json").write_text(json.dumps(result, indent=2) + "\n")
world.barrier()
