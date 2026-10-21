#!/workspace/.venvs/gpaw-mpi/bin/python
"""Run one frozen GPAW single point for the new prismatic DEV candidate."""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import re
import time
from pathlib import Path

import ase
import gpaw
import gpaw_data
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write
from gpaw import FermiDirac, GPAW, PW
from gpaw.dft import MSR1
from gpaw.mpi import world


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[5]


class SafeStop(Exception):
    pass


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


def atomic_json(path: Path, data: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, path)


parser = argparse.ArgumentParser()
parser.add_argument("--task-id", required=True)
parser.add_argument("--attempt-dir", required=True)
parser.add_argument("--resume", action="store_true")
args = parser.parse_args()

attempt = Path(args.attempt_dir).resolve()
if not attempt.is_dir():
    raise RuntimeError(f"Attempt directory must be prepared by supervisor: {attempt}")
cfg_path = ROOT / "SCIENTIFIC_CONFIG.json"
cfg_sha = sha256(cfg_path)
config = json.loads(cfg_path.read_text())
if int(config["resources"].get("disk_stop_guard_bytes", -1)) <= 0:
    raise RuntimeError("SCIENTIFIC_CONFIG is missing a positive disk_stop_guard_bytes guard")
manifest = json.loads((ROOT / "TASKS.json").read_text())
task = next(t for t in manifest["tasks"] if t["structure_id"] == args.task_id)
input_path = ROOT / task["input_path"]
if sha256(input_path) != task["input_sha256"]:
    raise RuntimeError("Input SHA mismatch")
if ase.__version__ != config["software"]["ASE_version"]:
    raise RuntimeError(f"ASE version mismatch: {ase.__version__}")
if gpaw.__version__ != config["software"]["GPAW_version"]:
    raise RuntimeError(f"GPAW version mismatch: {gpaw.__version__}")
if os.environ.get("GPAW_MPI_BACKEND") != config["parallel"]["GPAW_MPI_BACKEND"]:
    raise RuntimeError("GPAW MPI backend mismatch")
if os.environ.get("OMP_NUM_THREADS") != "1" or world.size != config["parallel"]["MPI_ranks"]:
    raise RuntimeError("MPI/OMP launch configuration mismatch")
setup_path = Path(gpaw_data.datapath()).resolve()
os.environ["GPAW_SETUP_PATH"] = str(setup_path)
for filename, digest in config["paw_sha256"].items():
    if sha256(setup_path / filename) != digest:
        raise RuntimeError(f"PAW SHA mismatch: {filename}")

expected = read(input_path, format="extxyz")
if len(expected) != task["atom_count"] or expected.pbc.tolist() != [True, True, True]:
    raise RuntimeError("Input atom-count/PBC check failed")
if geometry_sha(expected) != task["geometry_sha256"]:
    raise RuntimeError("Input geometry SHA mismatch")
if not np.array_equal(expected.arrays.get("source_atom_id"), np.arange(1, len(expected) + 1)):
    raise RuntimeError("Source atom IDs are missing or not contiguous")
if any(k in expected.info for k in ("energy", "free_energy")) or "forces" in expected.arrays:
    raise RuntimeError("Input must be unlabeled before the DFT calculation")

run_tag = "initial"
if args.resume:
    prior_runs = sorted(attempt.glob("RUN_MANIFEST_resume-*.json"))
    run_tag = f"resume-{len(prior_runs) + 1:02d}"
if world.rank == 0:
    with (attempt / f"RUN_MANIFEST_{run_tag}.json").open("w") as f:
        json.dump({
        "task_id": args.task_id,
        "status": "SCF_RUNNING",
        "input_path": task["input_path"],
        "input_sha256": task["input_sha256"],
        "geometry_sha256": task["geometry_sha256"],
        "config_sha256": cfg_sha,
        "runner_sha256": sha256(Path(__file__).resolve()),
        "role_candidate": task["role_candidate"],
        "source_family": task["source_family"],
        "source_cif_sha256": manifest["source_cif_sha256"],
        "resume": args.resume,
        "GPAW_version": gpaw.__version__,
        "ASE_version": ase.__version__,
        "PAW_sha256": config["paw_sha256"],
        "method": config["method"],
        "started_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        }, f, indent=2)
world.barrier()

log_path = attempt / ("gpaw.log" if not args.resume else f"gpaw_{run_tag}.log")
start_step = 0
if args.resume:
    cp_path = attempt / "checkpoint.gpw"
    cp_meta = json.loads((attempt / "CHECKPOINT_STATUS.json").read_text())
    if not cp_path.is_file() or sha256(cp_path) != cp_meta["checkpoint_sha256"]:
        raise RuntimeError("Resume checkpoint missing or SHA mismatch")
    calc = GPAW(str(cp_path), txt=str(log_path))
    atoms = calc.atoms.copy()
    # GPAW GPW preserves chemical ordering and exact geometry but may omit
    # ASE-only metadata arrays. Reattach source IDs only after verifying every
    # stored atom matches the SHA-bound frozen extxyz in order and coordinates.
    if (len(atoms) != len(expected)
            or atoms.get_chemical_symbols() != expected.get_chemical_symbols()
            or not np.array_equal(np.asarray(atoms.positions), np.asarray(expected.positions))
            or not np.array_equal(np.asarray(atoms.cell), np.asarray(expected.cell))
            or not np.array_equal(np.asarray(atoms.pbc), np.asarray(expected.pbc))):
        raise RuntimeError("Checkpoint geometry does not match requested input")
    expected_ids = expected.arrays["source_atom_id"]
    checkpoint_ids = atoms.arrays.get("source_atom_id")
    if checkpoint_ids is None:
        atoms.set_array("source_atom_id", expected_ids.copy())
    elif not np.array_equal(checkpoint_ids, expected_ids):
        raise RuntimeError("Checkpoint source atom IDs do not match frozen input")
    if geometry_sha(atoms) != task["geometry_sha256"]:
        raise RuntimeError("Checkpoint geometry/source-ID SHA does not match requested input")
    start_step = int(cp_meta["logical_scf_iteration"])
    remaining = int(config["method"]["maxiter"]) - start_step
    if remaining <= 0:
        raise RuntimeError("Checkpoint already reached configured SCF iteration cap")
    calc.dft.ibzwfs.make_sure_wfs_are_read_from_gpw_file()
    fd = os.open(cp_path, os.O_RDONLY)
    try:
        os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
    finally:
        os.close(fd)
    calc.dft.scf_loop.maxiter = remaining
else:
    atoms = expected.copy()
    method = config["method"]
    mix = method["mixer"]
    calc = GPAW(
        mode=PW(method["cutoff_eV"]),
        xc=method["xc"],
        kpts=tuple(method["kpoints_size"]),
        occupations=FermiDirac(method["fermi_dirac_width_eV"]),
        spinpol=method["spin_polarized"],
        convergence=method["scf_convergence"],
        maxiter=method["maxiter"],
        eigensolver={"name": "dav", "niter": 2},
        mixer=MSR1(
            nmaxold=mix["nmaxold"], beta=mix["beta"], reg=mix["reg"],
            gb_scale=mix["gb_scale"], max_A=mix["max_A"],
            trust_scale=mix["trust_scale"], soft_lim=mix["soft_lim"],
            hard_lim=mix["hard_lim"], weight=mix["weight"], sigma=mix["sigma"],
        ),
        parallel={"sl_auto": method["sl_auto"]},
        txt=str(log_path),
    )
    atoms.calc = calc

initial_geometry_sha = geometry_sha(expected)
last_checkpoint_step = start_step
safe_stop_requested = False


def scf_hook(ctx):
    global last_checkpoint_step, safe_stop_requested
    step = int(ctx.niter)
    logical_step = start_step + step
    if world.rank == 0:
        atomic_json(attempt / "SCF_PROGRESS.json", {
            "status": "SCF_RUNNING", "logical_scf_iteration": logical_step,
            "iteration_limit": config["method"]["maxiter"],
            "updated_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        })
    world.barrier()
    stop_path = attempt / "STOP_REQUEST.json"
    stop_reason = None
    if stop_path.exists():
        try:
            stop_reason = json.loads(stop_path.read_text()).get("reason", "SAFE_STOP_REQUESTED")
        except Exception:
            stop_reason = "SAFE_STOP_REQUESTED"
        safe_stop_requested = True
    due = (logical_step % 20 == 0) or safe_stop_requested
    if not due:
        return

    try:
        free = os.statvfs(attempt).f_bavail * os.statvfs(attempt).f_frsize
        old_bytes = (attempt / "checkpoint.gpw").stat().st_size if (attempt / "checkpoint.gpw").exists() else 0
        min_free = int(config["resources"]["disk_stop_guard_bytes"])
        # Keep enough room for a rotating checkpoint plus the configured safety reserve.
        estimated_new = int(config["resources"]["estimated_checkpoint_bytes"])
        if free < min_free + old_bytes + (0 if old_bytes else estimated_new):
            if world.rank == 0:
                atomic_json(attempt / "CHECKPOINT_WRITE_REFUSED.json", {
                    "status": "HOLD_INSUFFICIENT_DISK_HEADROOM",
                    "free_bytes": free, "old_checkpoint_bytes": old_bytes,
                    "estimated_new_checkpoint_bytes": estimated_new,
                    "required_reserve_bytes": min_free,
                    "logical_scf_iteration": logical_step,
                })
            world.barrier()
            raise SafeStop("CHECKPOINT_HEADROOM_GUARD")

        next_path = attempt / "checkpoint.next.gpw"
        if world.rank == 0 and next_path.exists():
            next_path.unlink()
        world.barrier()
        calc.write(str(next_path), mode="all", precision="single")
        world.barrier()
        write_status = np.array([0], dtype=np.int64)
        if world.rank == 0:
            if not next_path.is_file() or next_path.stat().st_size < 100_000_000:
                write_status[0] = 1
            else:
                digest = sha256(next_path)
                os.replace(next_path, attempt / "checkpoint.gpw")
                atomic_json(attempt / "CHECKPOINT_STATUS.json", {
                    "status": "WRITTEN_HASHED",
                    "checkpoint_path": "checkpoint.gpw",
                    "checkpoint_sha256": digest,
                    "checkpoint_bytes": (attempt / "checkpoint.gpw").stat().st_size,
                    "logical_scf_iteration": logical_step,
                    "config_sha256": cfg_sha,
                    "input_sha256": task["input_sha256"],
                    "geometry_sha256": task["geometry_sha256"],
                    "written_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                    "readback_validation": "WILL_BE_VERIFIED_BEFORE_RESUME",
                })
                last_checkpoint_step = logical_step
                fd = os.open(attempt / "checkpoint.gpw", os.O_RDONLY)
                try:
                    os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
                finally:
                    os.close(fd)
        world.broadcast(write_status, 0)
        world.barrier()
        if int(write_status[0]):
            raise SafeStop("CHECKPOINT_WRITE_FAILED")
    except SafeStop:
        raise
    except Exception as exc:
        if world.rank == 0:
            atomic_json(attempt / "CHECKPOINT_FAILURE.json", {
                "status": "CHECKPOINT_WRITE_FAILED", "exception_type": type(exc).__name__,
                "message": str(exc), "logical_scf_iteration": logical_step,
                "old_checkpoint_preserved": (attempt / "checkpoint.gpw").is_file(),
            })
        world.barrier()
        raise SafeStop("CHECKPOINT_WRITE_FAILED") from exc

    if safe_stop_requested:
        raise SafeStop(stop_reason or "SAFE_STOP_REQUESTED")


calc.hooks["scf_step"] = scf_hook
if args.resume:
    atoms.calc = calc
try:
    energy = float(atoms.get_potential_energy(force_consistent=False))
except SafeStop as exc:
    if world.rank == 0:
        progress = json.loads((attempt / "SCF_PROGRESS.json").read_text()) if (attempt / "SCF_PROGRESS.json").is_file() else {}
        atomic_json(attempt / "RUN_STATUS.json", {
            "status": "STOPPED_WITH_CHECKPOINT" if (attempt / "checkpoint.gpw").is_file() else "STOPPED_WITHOUT_CHECKPOINT",
            "scf_converged": False,
            "stop_reason": str(exc),
            "logical_scf_iteration": int(progress.get("logical_scf_iteration", last_checkpoint_step)),
            "last_valid_checkpoint_iteration": last_checkpoint_step,
            "checkpoint_sha256": json.loads((attempt / "CHECKPOINT_STATUS.json").read_text()).get("checkpoint_sha256")
                                  if (attempt / "CHECKPOINT_STATUS.json").is_file() else None,
            "labels_written": False,
            "ended_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        })
    world.barrier()
    raise SystemExit(75)
except Exception as exc:
    if world.rank == 0:
        atomic_json(attempt / "RUN_STATUS.json", {
            "status": "FAILED_SCF_OR_RUNTIME", "scf_converged": False,
            "exception_type": type(exc).__name__, "message": str(exc),
            "logical_scf_iteration": start_step + int(calc.get_number_of_iterations()),
            "checkpoint_sha256": json.loads((attempt / "CHECKPOINT_STATUS.json").read_text()).get("checkpoint_sha256")
                                  if (attempt / "CHECKPOINT_STATUS.json").is_file() else None,
            "labels_written": False,
            "ended_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        })
    world.barrier()
    raise

additional_steps = int(calc.get_number_of_iterations())
logical_steps = start_step + additional_steps
calc.log.fd.flush()
world.barrier()
marker = np.array([0, -1], dtype=np.int64)
if world.rank == 0:
    matches = re.findall(r"Converged in (\d+) steps", log_path.read_text(errors="replace"))
    marker[0] = int(bool(matches))
    marker[1] = int(matches[-1]) if matches else -1
world.broadcast(marker, 0)
if not marker[0] or int(marker[1]) != additional_steps:
    raise RuntimeError(f"Convergence marker mismatch: log={int(marker[1])}, GPAW={additional_steps}")
if logical_steps > int(config["method"]["maxiter"]):
    raise RuntimeError("Logical SCF iteration cap exceeded")
if not np.isfinite(energy):
    raise RuntimeError("Non-finite energy")

free_energy = float(atoms.get_potential_energy(force_consistent=True))
forces = np.asarray(atoms.get_forces(), dtype=float)
stress = np.asarray(atoms.get_stress(voigt=False), dtype=float)
if forces.shape != (len(expected), 3) or not np.isfinite(forces).all():
    raise RuntimeError("Invalid forces")
if stress.shape != (3, 3) or not np.isfinite(stress).all():
    raise RuntimeError("Invalid stress")
if geometry_sha(atoms) != initial_geometry_sha or geometry_sha(expected) != task["geometry_sha256"]:
    raise RuntimeError("Coordinates/cell/PBC/atom order changed during fixed-geometry SCF")

if world.rank == 0:
    labelled = atoms.copy()
    labelled.calc = SinglePointCalculator(labelled, energy=energy, free_energy=free_energy,
                                          forces=forces, stress=stress)
    write(attempt / "labelled.extxyz", labelled, format="extxyz")
    atomic_json(attempt / "RESULT.json", {
        "status": "SUCCESS", "scf_converged": True,
        "structure_id": args.task_id, "source_family": task["source_family"],
        "role_candidate": task["role_candidate"],
        "warning": "Different surface orientation from TRAIN but shares COD 9009647 bulk source; not an independent crystallographic-source blind test.",
        "scf_steps": logical_steps,
        "restart_start_step": start_step,
        "total_energy_eV": energy,
        "force_consistent_free_energy_eV": free_energy,
        "forces_eV_per_angstrom": forces.tolist(),
        "stress_eV_per_angstrom3": stress.tolist(),
        "stress_GPa": (stress * 160.21766208).tolist(),
        "max_force_eV_per_angstrom": float(np.linalg.norm(forces, axis=1).max()),
        "atom_count": len(atoms), "composition": task["composition"],
        "input_sha256": task["input_sha256"], "geometry_sha256": task["geometry_sha256"],
        "config_sha256": cfg_sha,
        "runner_sha256": sha256(Path(__file__).resolve()),
        "kpoints_type": config["method"]["kpoint_type"],
        "kpoints": config["method"]["kpoints_size"],
        "GPAW_version": gpaw.__version__, "ASE_version": ase.__version__,
        "units": {"energy": "eV", "forces": "eV/Angstrom", "stress": "eV/Angstrom^3"},
        "role_assignment_changed": False,
        "ended_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
    })
    atomic_json(attempt / "RUN_STATUS.json", {
        "status": "SUCCESS", "scf_converged": True,
        "scf_steps": logical_steps,
        "labels_written": True,
        "checkpoint_status": "retained temporarily until successful result archive is verified",
        "ended_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
    })
world.barrier()
