#!/usr/bin/env python3
"""Future fixed-geometry GPAW entry point; run only through the guarded shell draft."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import time

import numpy as np
from ase.io import read, write


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--workspace-root", required=True, type=Path)
    p.add_argument("--pack-root", required=True, type=Path)
    p.add_argument("--input-root", required=True, type=Path)
    p.add_argument("--run-root", required=True, type=Path)
    p.add_argument("--authorization", required=True, type=Path)
    p.add_argument("--owner", required=True, choices=("window-a", "window-b"))
    p.add_argument("--label", required=True)
    return p.parse_args()


args = parse_args()
workspace = args.workspace_root.resolve()
pack = args.pack_root.resolve()
input_root = args.input_root.resolve()
run_root = args.run_root.resolve()
auth_path = args.authorization.resolve()
manifest = json.loads((pack / "input_manifest.json").read_text())
authorization = json.loads(auth_path.read_text())
record = next((row for row in manifest["records"] if row["label"] == args.label), None)
if record is None:
    raise SystemExit(f"Unknown label: {args.label}")
if authorization.get("status") != "authorized" or authorization.get("reviewed_by") != "window-a":
    raise SystemExit("A role review/owner authorization is required")
if authorization.get("owner_assignments", {}).get(args.label) != args.owner:
    raise SystemExit("This owner is not assigned to the requested label")
if authorization.get("role_assignments", {}).get(args.label) != record["role"]:
    raise SystemExit("A-approved role does not match the prepared role")
if authorization.get("role_manifest_sha256") != manifest["source_geometry_bundle"]["role_manifest_sha256"]:
    raise SystemExit("A authorization refers to a different role manifest")
if authorization.get("input_sha256", {}).get(args.label) != record["input_sha256"]:
    raise SystemExit("A authorization refers to a different candidate input")
if record["role"] == "withheld_test":
    frozen = authorization.get("frozen_selection", {})
    if frozen.get("frozen_before_test_label_read") is not True or not all(
        frozen.get(key) for key in ("selected_model_sha256", "selection_record_sha256", "frozen_at_utc")
    ):
        raise SystemExit("A must freeze and hash the selected checkpoint before withheld-test label access")

input_path = input_root / Path(record["input"]).name
if digest(input_path) != record["input_sha256"]:
    raise SystemExit("Candidate input SHA256 mismatch")
atoms = read(input_path, format="extxyz")
for key in ("energy", "free_energy", "ref_energy", "pw_pbe_energy_ev"):
    if key in {k.lower() for k in atoms.info}:
        raise SystemExit(f"Candidate input contains unexpected label field: {key}")
for key in ("forces", "ref_forces", "pw_pbe_forces"):
    if key in {k.lower() for k in atoms.arrays}:
        raise SystemExit(f"Candidate input contains unexpected force field: {key}")
if len(atoms) != record["atom_count"] or atoms.get_chemical_symbols() != record["symbols_in_file_order"]:
    raise SystemExit("Candidate atom count/order differs from prepared manifest")
if np.asarray(atoms.arrays["lammps_id"], dtype=int).tolist() != record["atom_ids_in_file_order"]:
    raise SystemExit("Candidate atom IDs/order differ from prepared manifest")
if not np.allclose(atoms.cell.array, record["cell_A"], atol=1e-10, rtol=0) or np.asarray(atoms.pbc, dtype=bool).tolist() != record["pbc"]:
    raise SystemExit("Candidate cell or PBC differs from prepared manifest")

# Deliberately imported only after role authorization and geometry checks.
import gpaw_data
os.environ["GPAW_SETUP_PATH"] = str(gpaw_data.datapath())
from gpaw import GPAW, PW
from gpaw.mixer import Mixer
from gpaw.mpi import world
from gpaw.occupations import FermiDirac

if world.size != 4:
    raise SystemExit(f"Expected exactly four MPI ranks, got {world.size}")
out = run_root / args.label
if world.rank == 0:
    out.mkdir(parents=True, exist_ok=False)
world.barrier()

started = time.time()
progress_path = out / "progress.json"
state_path = out / "state.gpw"  # local scratch only; never copied to the compact archive
if world.rank == 0:
    progress_path.write_text(json.dumps({
        "label": args.label,
        "role": record["role"],
        "owner": args.owner,
        "status": "running",
        "iteration": 0,
        "source_sha256": record["input_sha256"],
        "mpi_ranks": world.size,
        "started_utc_epoch_s": started,
    }, indent=2) + "\n")

calc = GPAW(
    mode=PW(500),
    xc="PBE",
    kpts=(1, 1, 1),
    occupations=FermiDirac(0.1),
    mixer=Mixer(0.05, 8, 100),
    parallel={"sl_auto": True},
    convergence={"energy": 1e-5, "density": 1e-5, "eigenstates": 1e-5},
    maxiter=160,
    txt=str(out / "gpaw.log"),
)
atoms.calc = calc

def checkpoint(ctx):
    niter = int(ctx.niter)
    if niter % 10 == 0 and world.rank == 0:
        progress_path.write_text(json.dumps({
            "label": args.label,
            "role": record["role"],
            "owner": args.owner,
            "status": "running",
            "iteration": niter,
            "source_sha256": record["input_sha256"],
            "mpi_ranks": world.size,
            "elapsed_s": time.time() - started,
        }, indent=2) + "\n")
    if niter % 20 == 0:
        calc.write(str(state_path), mode="all")

calc.hooks["scf_step"] = checkpoint
converged = {"value": False}
calc.hooks["converged"] = lambda: converged.__setitem__("value", True)
energy = float(atoms.get_potential_energy(force_consistent=False))
free_energy = float(calc.results["free_energy"]) if "free_energy" in calc.results else None
forces = np.asarray(atoms.get_forces(), dtype=float)
iterations = int(calc.get_number_of_iterations())
finite = (
    np.isfinite(energy)
    and np.isfinite(forces).all()
    and (free_energy is None or np.isfinite(free_energy))
)
if world.rank == 0:
    summary = {
        "label": args.label,
        "role": record["role"],
        "owner": args.owner,
        "source": str(input_path),
        "source_sha256": record["input_sha256"],
        "atom_count": len(atoms),
        "formula": atoms.get_chemical_formula(),
        "symbols_in_file_order": atoms.get_chemical_symbols(),
        "atom_ids_in_file_order": np.asarray(atoms.arrays["lammps_id"], dtype=int).tolist(),
        "energy_eV_cell": energy,
        "energy_convention": "GPAW native extrapolated energy, force_consistent=False",
        "free_energy_eV_cell": free_energy,
        "free_energy_supported": free_energy is not None,
        "scf_converged": bool(converged["value"]),
        "finite": bool(finite),
        "scf_iterations": iterations,
        "mpi_ranks": world.size,
        "elapsed_s": time.time() - started,
        "method": {
            "GPAW_version": __import__("gpaw").__version__,
            "xc": "PBE",
            "basis": "plane wave",
            "cutoff_eV": 500,
            "kpts": [1, 1, 1],
            "smearing": "FermiDirac",
            "smearing_eV": 0.1,
            "parallel": {"sl_auto": True},
            "scf_thresholds": {"energy": 1e-5, "density": 1e-5, "eigenstates": 1e-5},
            "mixer": "Pulay beta=0.05, nmaxold=8, weight=100",
        },
        "withheld_test_policy": (
            "Mechanical archive integrity only; B must not model-score or analyze this label."
            if record["role"] == "withheld_test" else
            "Development-validation role; A review and assignment required."
        ),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    progress_path.write_text(json.dumps({
        "label": args.label,
        "role": record["role"],
        "owner": args.owner,
        "status": "complete" if converged["value"] and finite else "not_converged_or_nonfinite",
        "iteration": iterations,
        "elapsed_s": summary["elapsed_s"],
        "scf_converged": bool(converged["value"]),
        "finite": bool(finite),
        "source_sha256": record["input_sha256"],
        "mpi_ranks": world.size,
    }, indent=2) + "\n")
    if finite:
        atoms.info["config_type"] = args.label
        atoms.info["PW_PBE_energy_eV"] = energy
        atoms.info["source_method"] = "GPAW 26.7.0 PW-PBE 500 eV Gamma fixed-geometry single point"
        atoms.arrays["PW_PBE_forces"] = forces.copy()
        write(out / f"{args.label}_PW_PBE.extxyz", atoms, format="extxyz")
        calc.write(str(state_path), mode="all")
world.barrier()
if not converged["value"] or not finite:
    raise SystemExit(f"DFT did not produce a converged finite label: {args.label}")
