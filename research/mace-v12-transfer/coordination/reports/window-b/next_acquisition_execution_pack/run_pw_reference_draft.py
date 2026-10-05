#!/usr/bin/env python3
"""REVIEW ONLY: disabled future four-rank DFT runner, not a launcher."""
from hashlib import sha256
import argparse
import json
import os
from pathlib import Path
import sys
import time

# Hard stop before GPAW imports, filesystem output or any calculation.
# A must publish production integration and owner assignments; no flag/env can enable this draft.
raise SystemExit("blocked_pending_A_review_and_V15_screen: review-only draft; A must integrate and assign jobs before execution")

import numpy as np
from ase.io import read, write
parser = argparse.ArgumentParser()
parser.add_argument("--repo-root", required=True, type=Path)
parser.add_argument("--label", required=True)
parser.add_argument("--output-dir", required=True, type=Path)
args = parser.parse_args()
REPO = args.repo_root.resolve()
ROOT = REPO / "research/mace-v12-transfer/coordination/reports/window-b/next_acquisition_execution_pack"
OUT = args.output_dir.resolve()
LABEL = args.label
manifest = json.loads((ROOT / "input_manifest.json").read_text())
record = next((r for r in manifest["records"] if r["label"] == LABEL), None)
if record is None:
    raise SystemExit("Unknown planned label")
SOURCE = (ROOT / record["input"]).resolve()
if SOURCE.parent != ROOT / "inputs" or sha256(SOURCE.read_bytes()).hexdigest() != record["input_sha256"]:
    raise SystemExit("Input path/hash mismatch")
if not OUT.is_relative_to(Path("/tmp")) or OUT == Path("/tmp"):
    raise SystemExit("Working output must be an explicit new /tmp child; never write checkpoints to production")
try:
    import gpaw_data
    os.environ["GPAW_SETUP_PATH"] = str(gpaw_data.datapath())
    from gpaw import GPAW, PW
    from gpaw.mixer import Mixer
    from gpaw.mpi import world
    from gpaw.occupations import FermiDirac

    if world.size != 4:
        raise SystemExit(f"Expected 4 MPI ranks, got {world.size}")
    if world.rank == 0:
        if OUT.exists() or (OUT / "active_job.json").exists():
            raise ValueError("Existing output/checkpoint/active job: preserve and investigate")
        # Atomic claim of a new output directory rejects nonempty/empty old runs and active markers.
        OUT.mkdir(parents=True, exist_ok=False)
        (OUT / "active_job.json").write_text(json.dumps({"pid": os.getpid(), "label": LABEL, "input_sha256": record["input_sha256"]}) + "\n")
    world.barrier()

    atoms = read(SOURCE)
    if not np.all(atoms.pbc):
        raise SystemExit("Input must be periodic in all directions")
    if any(k in atoms.info for k in ("REF_energy", "PW_PBE_energy_eV", "energy", "free_energy")):
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
    free_energy_value = calc.results.get("free_energy")
    free_energy = None if free_energy_value is None else float(free_energy_value)
    if not np.isfinite(energy) or (free_energy is not None and not np.isfinite(free_energy)):
        raise ValueError("Nonfinite energy/free energy")
    forces = np.asarray(atoms.get_forces(), dtype=float)
    if forces.shape != (len(atoms), 3) or not np.isfinite(forces).all():
        raise ValueError("Invalid forces")
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
    if free_energy is not None:
        atoms.info["PW_PBE_free_energy_eV"] = free_energy
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
            "energy_convention": "native/extrapolated; get_potential_energy(force_consistent=False); free energy recorded separately",
            "scf_converged": converged,
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
            "note": "Future training acquisition, pending explicit A integration; not a blind test.",
        }
        (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        progress_path.write_text(json.dumps({
            "label": LABEL, "status": "complete" if converged else "not_converged",
            "iteration": iterations, "elapsed_s": summary["elapsed_s"],
            "scf_converged": converged, "source_sha256": record["input_sha256"],
            "mpi_ranks": world.size,
        }, indent=2) + "\n")
        print(json.dumps(summary, indent=2), flush=True)
    world.barrier()
    if world.rank == 0 and converged:
        (OUT / "active_job.json").unlink()
    if not converged:
        raise SystemExit(f"SCF did not converge for {LABEL}")
except Exception as error:
    if OUT.is_dir() and globals().get("world") is not None and world.rank == 0:
        (OUT / "diagnostics_failure.json").write_text(json.dumps({"label": LABEL, "error": type(error).__name__, "message": str(error), "input_sha256": record["input_sha256"]}, indent=2) + "\n")
    raise
