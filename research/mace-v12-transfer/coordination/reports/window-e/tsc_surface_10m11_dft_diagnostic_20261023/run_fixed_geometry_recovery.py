#!/usr/bin/env python3
"""Run one frozen-geometry Ti3SiC2 (10-11) surface diagnostic single point."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import time
from pathlib import Path

import ase
import gpaw as gpaw_module
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write
from gpaw import FermiDirac, GPAW, PW, restart
from gpaw.dft import MSR1
from gpaw.eigensolvers import Davidson
from gpaw.mpi import world
from gpaw_data import datapath


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def json_atomic(path: Path, data: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, path)


def geometry_sha(atoms) -> str:
    h = hashlib.sha256()
    h.update(np.asarray(atoms.numbers, dtype="<i4").tobytes())
    h.update(np.asarray(atoms.cell.array, dtype="<f8").tobytes())
    h.update(np.asarray(atoms.positions, dtype="<f8").tobytes())
    h.update(np.asarray(atoms.pbc, dtype=np.uint8).tobytes())
    h.update(np.asarray(atoms.arrays["source_atom_id"], dtype="<i8").tobytes())
    return h.hexdigest()


def free_bytes(path: Path) -> int:
    return shutil.disk_usage(path).free


parser = argparse.ArgumentParser()
parser.add_argument("--structure-id", required=True)
parser.add_argument("--input", required=True)
parser.add_argument("--run-dir", required=True)
parser.add_argument("--config", required=True)
parser.add_argument("--config-sha256", required=True)
parser.add_argument("--lock", required=True)
parser.add_argument("--resume", action="store_true")
args = parser.parse_args()

input_path = Path(args.input).resolve()
run_dir = Path(args.run_dir).resolve()
config_path = Path(args.config).resolve()
lock_path = Path(args.lock).resolve()
run_dir.mkdir(parents=True, exist_ok=True)
config = json.loads(config_path.read_text())
lock = json.loads(lock_path.read_text())
assert sha256(config_path) == args.config_sha256 == lock["config_sha256"]
assert sha256(Path(__file__).resolve()) == lock["runner_sha256"]
assert ase.__version__ == config["software"]["ASE_version"]
assert gpaw_module.__version__ == config["software"]["GPAW_version"]
assert os.environ.get("OMP_NUM_THREADS") == "1"
assert os.environ.get("GPAW_MPI_BACKEND") == "cgpaw"
assert world.size == 4
assert Path(os.environ["GPAW_SETUP_PATH"]).resolve() == Path(datapath()).resolve()
for paw_name, expected_sha in config["paw_sha256"].items():
    assert sha256(Path(datapath()) / paw_name) == expected_sha, f"PAW SHA mismatch: {paw_name}"

task = next(t for t in config["tasks"] if t["structure_id"] == args.structure_id)
assert sha256(input_path) == task["input_sha256"]
atoms = read(input_path, format="extxyz")
symbols = atoms.get_chemical_symbols()
assert len(atoms) == 69
assert {s: symbols.count(s) for s in sorted(set(symbols))} == {"C": 24, "Si": 12, "Ti": 33}
assert atoms.pbc.tolist() == [True, True, True]
ids = np.asarray(atoms.arrays.get("source_atom_id"), dtype=int)
assert ids.shape == (69,) and len(np.unique(ids)) == 69
assert "energy" not in atoms.info and "free_energy" not in atoms.info
assert "forces" not in atoms.arrays and "stress" not in atoms.info
before_geometry_sha = geometry_sha(atoms)
assert before_geometry_sha == task["geometry_sha256"]

status_path = run_dir / "RUN_STATUS.json"
checkpoint_latest = run_dir / "checkpoint_latest.json"
if (run_dir / "RESULT.json").exists():
    raise SystemExit("Existing result found; refusing duplicate calculation")

method = config["method"]
mix = dict(method["mixer"])
mix.pop("type", None)
started = time.monotonic()
prior_iterations = 0
checkpoint = None
checkpoint_meta = None
if args.resume:
    assert checkpoint_latest.is_file(), "resume requested without a committed checkpoint pointer"
    pointer = json.loads(checkpoint_latest.read_text())
    checkpoint = run_dir / pointer["checkpoint_file"]
    checkpoint_meta = run_dir / pointer["metadata_file"]
    assert checkpoint.is_file() and checkpoint_meta.is_file(), "committed checkpoint pair is incomplete"
    meta = json.loads(checkpoint_meta.read_text())
    assert meta["input_sha256"] == task["input_sha256"]
    assert meta["geometry_sha256"] == before_geometry_sha
    assert meta["config_sha256"] == args.config_sha256
    assert meta["checkpoint_sha256"] == sha256(checkpoint)
    assert pointer["checkpoint_sha256"] == meta["checkpoint_sha256"]
    prior_iterations = int(meta["completed_scf_iterations"])
    assert prior_iterations < int(method["maxiter"]), "checkpoint already reached maxiter"
    log_path = run_dir / f"gpaw_attempt_after_{prior_iterations:04d}.log"
    restart_atoms, calc = restart(str(checkpoint), txt=str(log_path))
    assert np.array_equal(restart_atoms.numbers, atoms.numbers)
    assert np.allclose(restart_atoms.positions, atoms.positions, atol=1e-12, rtol=0)
    assert np.allclose(restart_atoms.cell.array, atoms.cell.array, atol=1e-12, rtol=0)
    assert np.array_equal(restart_atoms.pbc, atoms.pbc)
    # GPAW 26.7's new ASE frontend does not expose maxiter through
    # ASECalculator.set().  Resume from the same serialized SCF state and set
    # only the remaining iteration budget on its SCF loop.  Across the initial
    # run and all resumes, the total hard ceiling remains method["maxiter"].
    remaining_iterations = int(method["maxiter"]) - prior_iterations
    assert remaining_iterations > 0, "checkpoint already reached maxiter"
    assert calc.dft.scf_loop.maxiter == int(method["maxiter"]), (
        "checkpoint SCF maxiter differs from the locked configuration")
    calc.dft.scf_loop.maxiter = remaining_iterations
    assert calc.dft.scf_loop.maxiter == remaining_iterations
    atoms.calc = calc
else:
    assert not checkpoint_latest.exists(), "checkpoint exists; use --resume after verification"
    assert not list(run_dir.glob("checkpoint_iter_*.gpw")), "orphan checkpoint requires manual audit"
    log_path = run_dir / "gpaw_attempt_0000.log"
    calc = GPAW(
        mode=PW(method["cutoff_eV"]),
        xc=method["xc"],
        kpts=tuple(method["kpoints_size"]),
        occupations=FermiDirac(method["fermi_dirac_width_eV"]),
        spinpol=method["spin_polarized"],
        eigensolver=Davidson(niter=method["eigensolver_niter"]),
        mixer=MSR1(**mix),
        convergence=method["scf_convergence"],
        maxiter=method["maxiter"],
        parallel={"sl_auto": method["sl_auto"]},
        txt=str(log_path),
    )
    atoms.calc = calc

status = {
    "structure_id": args.structure_id,
    "status": "RUNNING",
    "scientific_role": "UNASSIGNED_SAME_SOURCE_FACET_DIAGNOSTIC",
    "input_sha256": task["input_sha256"],
    "geometry_sha256": before_geometry_sha,
    "config_sha256": args.config_sha256,
    "runner_sha256": lock["runner_sha256"],
    "resume_from_completed_iterations": prior_iterations,
    "mpi_ranks": world.size,
    "omp_threads": 1,
    "checkpoint_interval_scf_steps": config["checkpoint"]["interval_scf_steps"],
}
if world.rank == 0:
    json_atomic(status_path, status)

last_iteration = prior_iterations
converged = {"value": False}
stop_reason = {"value": None}
checkpoint_size = checkpoint.stat().st_size if checkpoint is not None else 0


def save_checkpoint(iteration: int) -> None:
    global checkpoint_size, checkpoint, checkpoint_meta
    reserve = int(config["checkpoint"]["minimum_free_reserve_bytes"])
    transient_reserve = int(config["checkpoint"]["minimum_transient_free_bytes"])
    available = free_bytes(run_dir)
    # A replacement checkpoint is written beside the previous one. Keep a reserve
    # for logs/results; never remove the last valid checkpoint to make room.
    required = checkpoint_size + transient_reserve if checkpoint_size else reserve
    if available < required:
        stop_reason["value"] = "CHECKPOINT_SPACE_GUARD"
        raise RuntimeError("Insufficient disk headroom for safe checkpoint replacement")
    stem = f"checkpoint_iter_{int(iteration):04d}"
    checkpoint_next = run_dir / f"{stem}.next.gpw"
    checkpoint_new = run_dir / f"{stem}.gpw"
    meta_new = run_dir / f"{stem}.json"
    if checkpoint_next.exists():
        checkpoint_next.unlink()
    calc.write(str(checkpoint_next), mode="all")
    world.barrier()
    outcome = np.array([1, 0, 0], dtype=np.int64)
    if world.rank == 0:
        try:
            with checkpoint_next.open("rb") as f:
                os.fsync(f.fileno())
            size = checkpoint_next.stat().st_size
            digest = sha256(checkpoint_next)
            if hasattr(os, "posix_fadvise"):
                fd = os.open(checkpoint_next, os.O_RDONLY)
                try:
                    os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
                finally:
                    os.close(fd)
            free_after_write = free_bytes(run_dir)
            if free_after_write < transient_reserve:
                raise OSError("checkpoint would consume the transient disk headroom")
            if free_after_write + checkpoint_size < reserve:
                raise OSError("checkpoint replacement would consume persistent disk reserve")
            os.replace(checkpoint_next, checkpoint_new)
            meta = {
                "structure_id": args.structure_id,
                "completed_scf_iterations": int(iteration),
                "input_sha256": task["input_sha256"],
                "geometry_sha256": before_geometry_sha,
                "config_sha256": args.config_sha256,
                "checkpoint_sha256": digest,
                "checkpoint_bytes": size,
                "status": "WRITTEN_SHA_VERIFIED_RELOAD_PENDING",
            }
            json_atomic(meta_new, meta)
            json_atomic(checkpoint_latest, {
                "checkpoint_file": checkpoint_new.name,
                "metadata_file": meta_new.name,
                "checkpoint_sha256": digest,
                "completed_scf_iterations": int(iteration),
            })
            old_checkpoint, old_meta = checkpoint, checkpoint_meta
            checkpoint_size = size
            checkpoint, checkpoint_meta = checkpoint_new, meta_new
            outcome[1] = size
            # Only delete the previous checkpoint after the new data, metadata,
            # and atomic latest-pointer commit all exist and match.
            if old_checkpoint is not None and old_checkpoint != checkpoint and old_checkpoint.exists():
                old_checkpoint.unlink()
            if old_meta is not None and old_meta != checkpoint_meta and old_meta.exists():
                old_meta.unlink()
            status.update({
                "status": "RUNNING_CHECKPOINT_SAVED",
                "scf_iteration": int(iteration),
                "checkpoint_sha256": digest,
                "checkpoint_bytes": size,
                "free_disk_bytes_after_checkpoint": free_bytes(run_dir),
            })
            json_atomic(status_path, status)
        except Exception:
            outcome[0] = 0
            outcome[1] = 1
            if checkpoint_next.exists():
                checkpoint_next.unlink()
    world.broadcast(outcome, 0)
    world.barrier()
    if outcome[0] != 1:
        stop_reason["value"] = "CHECKPOINT_WRITE_OR_SPACE_FAILURE"
        raise RuntimeError("Checkpoint could not be safely committed; previous checkpoint retained")
    checkpoint_size = int(outcome[1])
    checkpoint, checkpoint_meta = checkpoint_new, meta_new


def on_scf_step(ctx) -> None:
    global last_iteration
    last_iteration = prior_iterations + int(ctx.niter)
    elapsed = time.monotonic() - started
    try:
        current_memory = int(Path("/sys/fs/cgroup/memory.current").read_text())
    except Exception:
        current_memory = 0
    due = last_iteration % int(config["checkpoint"]["interval_scf_steps"]) == 0
    wall_stop = elapsed >= int(config["resources"]["walltime_seconds"]) - 600
    memory_stop = current_memory >= int(config["resources"]["safe_memory_stop_bytes"])
    if due or wall_stop or memory_stop:
        already_saved = (checkpoint_meta is not None and checkpoint_meta.exists()
                         and last_iteration == int(json.loads(checkpoint_meta.read_text())["completed_scf_iterations"]))
        if not already_saved:
            save_checkpoint(last_iteration)
    if wall_stop:
        stop_reason["value"] = "WALLTIME_SAFE_STOP"
        raise RuntimeError("Safe stop requested before walltime limit")
    if memory_stop:
        stop_reason["value"] = "MEMORY_SAFE_STOP"
        raise RuntimeError("Safe stop requested before cgroup memory limit")


calc.hooks["scf_step"] = on_scf_step
calc.hooks["converged"] = lambda: converged.__setitem__("value", True)

try:
    energy = float(atoms.get_potential_energy(force_consistent=False))
    if hasattr(calc, "log") and getattr(calc.log, "fd", None):
        calc.log.fd.flush()
    log_text = log_path.read_text(errors="replace")
    matches = re.findall(r"Converged in (\d+) steps", log_text)
    # GPAW 26.7 restart() logs an authoritative "Converged in N steps" line;
    # the ASE converged hook can remain false on a resumed calculator.
    if not matches:
        raise RuntimeError("GPAW log has no converged SCF marker")
    iterations = prior_iterations + int(matches[-1])
    if not np.isfinite(energy):
        raise RuntimeError("Non-finite total energy")
    free_energy = float(atoms.get_potential_energy(force_consistent=True))
    forces = np.asarray(atoms.get_forces(), dtype=float)
    stress = np.asarray(atoms.get_stress(voigt=False), dtype=float)
    if forces.shape != (69, 3) or stress.shape != (3, 3):
        raise RuntimeError(f"Unexpected force/stress shape: {forces.shape}, {stress.shape}")
    if not (np.isfinite(free_energy) and np.isfinite(forces).all() and np.isfinite(stress).all()):
        raise RuntimeError("Non-finite free energy, force, or stress")
    if geometry_sha(atoms) != before_geometry_sha:
        raise RuntimeError("Geometry changed during fixed-geometry calculation")

    if world.rank == 0:
        labelled = atoms.copy()
        labelled.calc = SinglePointCalculator(labelled, energy=energy, free_energy=free_energy,
                                              forces=forces, stress=stress)
        output_path = run_dir / "labelled.extxyz"
        write(output_path.with_suffix(".extxyz.tmp"), labelled, format="extxyz")
        os.replace(output_path.with_suffix(".extxyz.tmp"), output_path)
        result = {
            "structure_id": args.structure_id,
            "status": "SUCCESS",
            "scf_converged": True,
            "scf_iterations": iterations,
            "convergence_detected_by_gpaw_log": True,
            "convergence_callback_fired": bool(converged["value"]),
            "energy_eV": energy,
            "force_consistent_free_energy_eV": free_energy,
            "forces_eV_per_angstrom": forces.tolist(),
            "stress_eV_per_angstrom3": stress.tolist(),
            "max_force_eV_per_angstrom": float(np.linalg.norm(forces, axis=1).max()),
            "atom_count": len(atoms),
            "composition": {s: symbols.count(s) for s in sorted(set(symbols))},
            "input_sha256": task["input_sha256"],
            "geometry_sha256": before_geometry_sha,
            "config_sha256": args.config_sha256,
            "checkpoint_sha256_at_last_save": (json.loads(checkpoint_meta.read_text()).get("checkpoint_sha256") if checkpoint_meta is not None and checkpoint_meta.exists() else None),
            "checkpoint_retained_after_success": False,
            "gpaw_log_sha256": sha256(log_path),
            "role": "UNASSIGNED_SAME_SOURCE_FACET_DIAGNOSTIC",
            "GPAW_version": gpaw_module.__version__,
            "ASE_version": ase.__version__,
            "kpoints": method["kpoints_size"],
            "cutoff_eV": method["cutoff_eV"],
            "xc": method["xc"],
            "elapsed_seconds": time.monotonic() - started,
        }
        json_atomic(run_dir / "RESULT.json", result)
        # The converged labels are durable now; remove only this run's temporary
        # recovery state, as requested for successful jobs.
        for p in run_dir.glob("checkpoint_iter_*"):
            if p.is_file():
                p.unlink()
        if checkpoint_latest.exists():
            checkpoint_latest.unlink()
        status.update({"status": "SUCCESS", "scf_iterations": iterations,
                       "checkpoint_cleanup_after_success": True})
        json_atomic(status_path, status)
        files = sorted(p for p in run_dir.iterdir() if p.is_file() and p.name != "SHA256SUMS.txt")
        (run_dir / "SHA256SUMS.txt").write_text("".join(f"{sha256(p)}  {p.name}\n" for p in files))
    world.barrier()
except Exception as exc:
    if world.rank == 0:
        status.update({
            "status": "SAFE_STOP_OR_FAILURE_NO_LABELS",
            "scf_converged": False,
            "scf_iteration": last_iteration,
            "stop_reason": stop_reason["value"] or type(exc).__name__,
            "error": f"{type(exc).__name__}: {exc}",
            "elapsed_seconds": time.monotonic() - started,
            "checkpoint_present": checkpoint is not None and checkpoint.exists(),
            "checkpoint_sha256": sha256(checkpoint) if checkpoint is not None and checkpoint.exists() else None,
            "gpaw_log": log_path.name,
        })
        json_atomic(status_path, status)
        (run_dir / "ERROR.txt").write_text(status["error"] + "\n")
        files = sorted(p for p in run_dir.iterdir() if p.is_file() and p.name != "SHA256SUMS.txt")
        (run_dir / "SHA256SUMS.txt").write_text("".join(f"{sha256(p)}  {p.name}\n" for p in files))
    raise
