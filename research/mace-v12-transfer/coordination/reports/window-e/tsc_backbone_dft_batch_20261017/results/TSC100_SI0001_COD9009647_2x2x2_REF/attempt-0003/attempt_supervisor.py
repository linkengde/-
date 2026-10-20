#!/usr/bin/env python3
"""One bounded, independent SCF attempt with durable latest-checkpoint retention."""
from __future__ import annotations

import datetime as dt
import argparse
import hashlib
import importlib.metadata
import json
import os
import shutil
import signal
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

import ase
import gpaw
import gpaw_data
import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parent
TASK = "TSC100_SI0001_COD9009647_2x2x2_REF"
LIMIT = 18000
STOP_RESERVE = 600
MEMORY_GUARD = 28 * 1024**3
MIN_RUNTIME_FREE = 6 * 1024**3
ENV_ROOT = Path("/workspace/.venvs/gpaw-mpi")
PYTHON = ENV_ROOT / "bin/python"
INPUT = ROOT / "input.extxyz"
CONFIG = ROOT / "diagnostic_config.json"
LOCK = ROOT / "EXECUTION_LOCK.json"
WORKER = ROOT / "attempt_worker.py"
PROFILE = ROOT / "RESOURCE_PROFILE.jsonl"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def utc() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def read_int(name: str) -> int:
    value = Path("/sys/fs/cgroup", name).read_text().strip()
    return 2**63 - 1 if value == "max" else int(value)


def geometry_sha(atoms) -> str:
    h = hashlib.sha256()
    h.update(" ".join(atoms.get_chemical_symbols()).encode())
    h.update(np.asarray(atoms.cell, dtype="<f8").tobytes(order="C"))
    h.update(np.asarray(np.round(atoms.positions, 8), dtype="<f8").tobytes(order="C"))
    h.update(np.asarray(atoms.pbc, dtype=np.uint8).tobytes())
    h.update(np.asarray(atoms.arrays["source_atom_id"], dtype="<i8").tobytes())
    return h.hexdigest()


def active_compute_processes():
    found = []
    for proc in Path("/proc").iterdir():
        if not proc.name.isdigit():
            continue
        try:
            comm = (proc / "comm").read_text().strip().lower()
            state = next(x.split()[1] for x in (proc / "status").read_text().splitlines()
                         if x.startswith("State:"))
            if state != "Z" and comm in {"gpaw", "mpiexec", "mpirun", "orted", "prte"}:
                found.append({"pid": int(proc.name), "comm": comm, "state": state})
        except (OSError, StopIteration):
            pass
    return found


def atomic_json(path: Path, data):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, path)


def preflight():
    cfg = json.loads(CONFIG.read_text())
    lock = json.loads(LOCK.read_text())
    expected_config_sha = lock["config_sha256"]
    expected_lock_sha = (ROOT / "EXECUTION_LOCK.sha256").read_text().split()[0]
    setup_path = Path(gpaw_data.datapath()).resolve()
    cpu = Path("/sys/fs/cgroup/cpu.max").read_text().split()
    free = shutil.disk_usage(ROOT).free
    memory_current = read_int("memory.current")
    checks = {
        "execution_authorized": lock.get("execution_authorized") is True,
        "one_attempt_only": lock.get("automatic_retry") is False,
        "input_sha256": sha(INPUT) == cfg["input_sha256"],
        "config_sha256": sha(CONFIG) == expected_config_sha,
        "lock_sha256": sha(LOCK) == expected_lock_sha,
        "worker_sha256": sha(WORKER) == lock["worker_sha256"],
        "supervisor_sha256": sha(Path(__file__)) == lock["supervisor_sha256"],
        "GPAW_version": gpaw.__version__ == "26.7.0",
        "ASE_version": ase.__version__ == "3.29.0",
        "gpaw_data_version": importlib.metadata.version("gpaw-data") == "1.2.1",
        "PAW_hashes": all(sha(setup_path / name) == digest
                          for name, digest in cfg["paw_sha256"].items()),
        "cpu_quota_4": len(cpu) == 2 and int(cpu[0]) / int(cpu[1]) >= 4,
        "memory_limit_32GiB": read_int("memory.max") >= 32 * 1024**3,
        "memory_below_guard": memory_current < MEMORY_GUARD,
        "checkpoint_space_budget": free >= cfg["checkpoint_policy"]["checkpoint_and_headroom_required_bytes"],
        "runtime_headroom_available": free >= MIN_RUNTIME_FREE,
        "no_other_DFT_process": not active_compute_processes(),
    }
    atoms = read(INPUT, format="extxyz")
    composition = dict(Counter(atoms.get_chemical_symbols()))
    checks["structure_identity"] = (
        len(atoms) == 100 and composition == {"Ti": 48, "Si": 20, "C": 32}
        and atoms.pbc.tolist() == [True, True, True]
        and np.array_equal(atoms.arrays.get("source_atom_id"), np.arange(1, 101))
        and geometry_sha(atoms) == cfg["geometry_sha256"]
        and atoms.info.get("surface_termination") == "Si"
        and not any(k in atoms.info for k in ("energy", "free_energy", "stress"))
        and not any(k in atoms.arrays for k in ("forces", "stress")))
    env = os.environ.copy()
    env.update({"GPAW_SETUP_PATH": str(setup_path), "GPAW_MPI_BACKEND": "cgpaw",
                "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
                "MKL_NUM_THREADS": "1"})
    mpi_code = ("from gpaw.mpi import world; import gpaw,ase,socket; "
                "print(world.rank,world.size,socket.gethostname(),gpaw.__version__,ase.__version__)")
    probe = subprocess.run([shutil.which("mpiexec") or "mpiexec", "-n", "4",
                            str(PYTHON), "-c", mpi_code], env=env, text=True,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True)
    rows = [line.split() for line in probe.stdout.splitlines() if line.strip()]
    checks["MPI_4_rank_consistency"] = (
        len(rows) == 4 and {int(row[0]) for row in rows} == {0, 1, 2, 3}
        and all(int(row[1]) == 4 and row[3] == "26.7.0" and row[4] == "3.29.0"
                for row in rows) and len({row[2] for row in rows}) == 1)
    details = {
        "utc": utc(), "attempt_id": "attempt-0003", "checks": checks,
        "workspace_free_bytes": free, "required_checkpoint_space_bytes": cfg["checkpoint_policy"]["checkpoint_and_headroom_required_bytes"],
        "memory_current_bytes": memory_current, "memory_limit_bytes": read_int("memory.max"),
        "cpu_max": " ".join(cpu), "active_compute_processes": active_compute_processes(),
        "GPAW_version": gpaw.__version__, "ASE_version": ase.__version__,
        "gpaw_data_version": importlib.metadata.version("gpaw-data"),
        "GPAW_SETUP_PATH": str(setup_path), "MPI_probe": rows,
        "kpoints": cfg["method"]["kpoints_size"], "mixer": cfg["method"]["mixer"]}
    atomic_json(ROOT / "TARGET_PREFLIGHT.json", {
        "status": "PREFLIGHT_PASS" if all(checks.values()) else "PREFLIGHT_HOLD",
        "details": details, "walltime_limit_seconds": LIMIT,
        "checkpoint_policy": cfg["checkpoint_policy"]})
    return checks, details, env


def terminate(proc):
    if proc.poll() is not None:
        return proc.wait()
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        return proc.wait(timeout=35)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        return proc.wait()


def sample(elapsed):
    row = {"utc": utc(), "elapsed_seconds": round(elapsed, 2),
           "memory_current_bytes": read_int("memory.current"),
           "workspace_free_bytes": shutil.disk_usage(ROOT).free}
    with PROFILE.open("a") as f:
        f.write(json.dumps(row, sort_keys=True) + "\n")
        f.flush()
    return row


def write_checksums():
    lines = []
    for path in sorted(ROOT.rglob("*")):
        if (not path.is_file() or path == ROOT / "SHA256SUMS.txt"
                or path.suffix == ".gpw" or path.name.endswith(".tmp")):
            continue
        lines.append(f"{sha(path)}  ./{path.relative_to(ROOT).as_posix()}")
    (ROOT / "SHA256SUMS.txt").write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()
    checks, details, env = preflight()
    print(json.dumps({"status": "PREFLIGHT_PASS" if all(checks.values()) else "PREFLIGHT_HOLD",
                      "checks": checks, "details": details}, indent=2), flush=True)
    if not all(checks.values()):
        return 2
    if args.preflight_only:
        return 0
    cfg = json.loads(CONFIG.read_text())
    lock = json.loads(LOCK.read_text())
    env.update({"ATTEMPT_ID": "attempt-0003", "TMPDIR": str(ROOT / "checkpoints"),
                "TMP": str(ROOT / "checkpoints"),
                "TASK_DEADLINE_MONOTONIC": str(time.monotonic() + LIMIT)})
    lock_sha = sha(LOCK)
    cmd = [shutil.which("mpiexec") or "mpiexec", "--bind-to", "core", "--map-by", "core",
           "-n", "4", str(PYTHON), str(WORKER), "--config", str(CONFIG),
           "--config-sha256", lock["config_sha256"], "--execution-lock", str(LOCK),
           "--execution-lock-sha256", lock_sha]
    atomic_json(ROOT / "RUN_MANIFEST.json", {
        "status": "RUNNING", "attempt_id": "attempt-0003", "structure_id": TASK,
        "input_sha256": cfg["input_sha256"], "geometry_sha256": cfg["geometry_sha256"],
        "config_sha256": lock["config_sha256"], "execution_lock_sha256": lock_sha,
        "worker_sha256": lock["worker_sha256"], "supervisor_sha256": lock["supervisor_sha256"],
        "C_review_commit": lock["C_review_commit"], "base_config_sha256": lock["base_config_sha256"],
        "method": cfg["method"], "checkpoint_policy": cfg["checkpoint_policy"],
        "environment": details, "command": cmd, "created_utc": utc()})
    atomic_json(ROOT / "RUN_STATUS.json", {"status": "RUNNING", "attempt_id": "attempt-0003",
        "scf_converged": False, "labels_written": False, "automatic_retry": False,
        "started_utc": utc(), "walltime_limit_seconds": LIMIT})
    with (ROOT / "stdout.log").open("wb") as out, (ROOT / "stderr.log").open("wb") as err:
        proc = subprocess.Popen(cmd, cwd=ROOT, env=env, stdin=subprocess.DEVNULL,
                                stdout=out, stderr=err, start_new_session=True)
        atomic_json(ROOT / "supervisor_started.json", {"pid": proc.pid, "utc": utc(),
            "walltime_limit_seconds": LIMIT, "stop_reserve_seconds": STOP_RESERVE})
        started = time.monotonic()
        max_memory = 0
        min_free = shutil.disk_usage(ROOT).free
        stop_reason = None
        stop_requested = False
        while proc.poll() is None:
            elapsed = time.monotonic() - started
            row = sample(elapsed)
            max_memory = max(max_memory, row["memory_current_bytes"])
            min_free = min(min_free, row["workspace_free_bytes"])
            if not stop_requested and row["memory_current_bytes"] >= MEMORY_GUARD:
                (ROOT / "STOP_REQUEST").write_text("memory guard approached; checkpoint and safe stop\n")
                stop_requested = True
                stop_reason = "SAFE_MEMORY_GUARD_STOP"
            elif not stop_requested and row["workspace_free_bytes"] < MIN_RUNTIME_FREE:
                (ROOT / "STOP_REQUEST").write_text("disk headroom low; checkpoint and safe stop\n")
                stop_requested = True
                stop_reason = "SAFE_DISK_HEADROOM_STOP"
            elif not stop_requested and elapsed >= LIMIT - STOP_RESERVE:
                (ROOT / "STOP_REQUEST").write_text("walltime reserve reached; checkpoint and safe stop\n")
                stop_requested = True
                stop_reason = "SAFE_WALLTIME_STOP"
            if elapsed >= LIMIT:
                stop_reason = stop_reason or "WALLTIME_HARD_STOP"
                break
            time.sleep(5)
        if proc.poll() is None:
            exit_code = terminate(proc)
        else:
            exit_code = proc.wait()
    status = json.loads((ROOT / "SCF_CONVERGENCE.json").read_text()) if (ROOT / "SCF_CONVERGENCE.json").exists() else None
    success = exit_code == 0 and status is not None and status.get("scf_converged") is True and (ROOT / "RESULT.json").exists()
    cp = json.loads((ROOT / "CHECKPOINT_STATUS.json").read_text()) if (ROOT / "CHECKPOINT_STATUS.json").exists() else None
    final = {"status": "SUCCESS" if success else "FAILED_SCF",
        "attempt_id": "attempt-0003", "structure_id": TASK, "exit_code": exit_code,
        "stop_reason": stop_reason, "scf_converged": bool(success),
        "scf_steps": status.get("scf_iterations") if status else None,
        "labels_written": bool(success), "latest_checkpoint": cp,
        "minimum_workspace_free_bytes_observed": min_free,
        "maximum_cgroup_memory_current_bytes_observed": max_memory,
        "ended_utc": utc(), "automatic_retry": False}
    atomic_json(ROOT / "RUN_STATUS.json", final)
    atomic_json(ROOT / "RUN_SUMMARY.json", {"status": final["status"],
        "stop_reason": stop_reason, "result_sha256": sha(ROOT / "RESULT.json") if success else None,
        "checkpoint_sha256": cp.get("sha256") if cp else None,
        "checkpoint_step": cp.get("step") if cp else None})
    if success:
        # After complete result files are written and hashed, no restart checkpoint is needed.
        checkpoint = ROOT / "checkpoints/latest.gpw"
        if checkpoint.exists():
            checkpoint.unlink()
        atomic_json(ROOT / "CHECKPOINT_CLEANUP.json", {"status": "REMOVED_AFTER_VERIFIED_SUCCESS",
            "utc": utc(), "reason": "converged labels and result manifest exist"})
    write_checksums()
    print(json.dumps(final, indent=2))
    return 0 if success else 3


if __name__ == "__main__":
    raise SystemExit(main())
