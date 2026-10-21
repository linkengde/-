#!/usr/bin/env python3
"""Single-point supervisor with bounded walltime, memory, disk and checkpoint guards."""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[5]
MAX_WALLTIME = 18000
SOFT_WALLTIME = 17400
MEMORY_GUARD = 30064771072
MEMORY_REQUEST_STOP = 27917287424
DISK_GUARD = 4294967296
DISK_REQUEST_STOP = 6442450944


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_int(path: Path) -> int:
    return int(path.read_text().strip())


def atomic_json(path: Path, data: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, path)


def free_bytes(path: Path) -> int:
    st = os.statvfs(path)
    return st.f_bavail * st.f_frsize


def active_compute_processes() -> list[dict]:
    rows = []
    for p in Path("/proc").iterdir():
        if not p.name.isdigit():
            continue
        try:
            lines = (p / "status").read_text(errors="replace").splitlines()
            name = next(x.split("\t", 1)[1].strip() for x in lines if x.startswith("Name:"))
            state = next(x.split()[1] for x in lines if x.startswith("State:"))
            if state != "Z" and name in {"gpaw", "mpiexec", "mpirun", "prterun"}:
                rows.append({"pid": int(p.name), "name": name, "state": state})
        except (OSError, StopIteration, ValueError):
            continue
    return rows


parser = argparse.ArgumentParser()
parser.add_argument("--task-id", required=True)
parser.add_argument("--resume", action="store_true")
parser.add_argument("--preflight-only", action="store_true")
parser.add_argument("--attempt-number", type=int)
args = parser.parse_args()
manifest = json.loads((ROOT / "TASKS.json").read_text())
task = next(t for t in manifest["tasks"] if t["structure_id"] == args.task_id)
task_results = ROOT / "results" / args.task_id
existing_attempts = sorted(task_results.glob("attempt-[0-9][0-9][0-9][0-9]")) if task_results.exists() else []
if args.attempt_number is not None:
    attempt_number = args.attempt_number
elif args.resume:
    resumable = [p for p in existing_attempts if (p / "checkpoint.gpw").is_file()]
    if not resumable:
        raise RuntimeError("--resume requested but no existing attempt has a checkpoint")
    attempt_number = max(int(p.name.split("-")[1]) for p in resumable)
else:
    attempt_number = max((int(p.name.split("-")[1]) for p in existing_attempts), default=0) + 1
attempt = task_results / f"attempt-{attempt_number:04d}"
if attempt_number < 1:
    raise RuntimeError("attempt number must be >= 1")
attempt.mkdir(parents=True, exist_ok=True)
if any(attempt.iterdir()) and not args.resume:
    raise RuntimeError(f"Attempt path is not empty: {attempt}")
if args.resume and not (attempt / "checkpoint.gpw").is_file():
    raise RuntimeError("--resume requested without a current checkpoint")
if args.resume and (attempt / "RESULT.json").is_file():
    raise RuntimeError("Refusing to resume an attempt that already has a final result")

resume_number = len(list(attempt.glob("SUPERVISOR_STARTED_resume-*.json"))) + 1 if args.resume else 0
run_tag = f"resume-{resume_number:02d}" if args.resume else "initial"

cfg = json.loads((ROOT / "SCIENTIFIC_CONFIG.json").read_text())
if int(cfg["resources"].get("disk_stop_guard_bytes", -1)) != DISK_GUARD:
    raise RuntimeError("SCIENTIFIC_CONFIG disk_stop_guard_bytes missing or inconsistent with supervisor")
setup = Path("/workspace/.venvs/gpaw-mpi/lib/python3.12/site-packages/gpaw_data/setups")
start_free = free_bytes(Path("/workspace"))
preflight = {
    "status": "PREFLIGHT_PASS",
    "task_id": args.task_id,
    "input_sha256": sha256(ROOT / task["input_path"]),
    "expected_input_sha256": task["input_sha256"],
    "config_sha256": sha256(ROOT / "SCIENTIFIC_CONFIG.json"),
    "runner_sha256": sha256(ROOT / "run_point.py"),
    "supervisor_sha256": sha256(Path(__file__).resolve()),
    "source_cif_sha256": manifest["source_cif_sha256"],
    "python": sys.version.split()[0],
    "GPAW_version": subprocess.check_output([sys.executable, "-c", "import gpaw; print(gpaw.__version__)"], text=True).strip(),
    "ASE_version": subprocess.check_output([sys.executable, "-c", "import ase; print(ase.__version__)"], text=True).strip(),
    "GPAW_SETUP_PATH": str(setup),
    "PAW_hashes_match": all(sha256(setup / f) == h for f, h in cfg["paw_sha256"].items()),
    "cpu_max": Path("/sys/fs/cgroup/cpu.max").read_text().strip(),
    "cpu_quota_cores": (lambda q: "max" if q[0] == "max" else round(int(q[0]) / int(q[1]), 3))(
        Path("/sys/fs/cgroup/cpu.max").read_text().split()),
    "memory_limit_bytes": read_int(Path("/sys/fs/cgroup/memory.max")),
    "memory_current_bytes": read_int(Path("/sys/fs/cgroup/memory.current")),
    "workspace_free_bytes": start_free,
    "workspace_free_gib": round(start_free / 1024**3, 3),
    "disk_guard_bytes": DISK_GUARD,
    "disk_stop_guard_bytes": int(cfg["resources"]["disk_stop_guard_bytes"]),
    "minimum_start_free_bytes": 13_000_000_000,
    "estimated_checkpoint_bytes": cfg["resources"]["estimated_checkpoint_bytes"],
    "active_compute_processes": active_compute_processes(),
    "mpi_version": subprocess.check_output(["mpiexec", "--version"], text=True, stderr=subprocess.STDOUT).splitlines()[0],
    "mpi_backend": "cgpaw",
    "omp_threads": 1,
    "resuming": args.resume,
    "git_branch": subprocess.check_output(["git", "-C", str(REPO), "branch", "--show-current"], text=True).strip(),
    "git_head": subprocess.check_output(["git", "-C", str(REPO), "rev-parse", "HEAD"], text=True).strip(),
    "checked_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
}
if preflight["GPAW_version"] != cfg["software"]["GPAW_version"] or preflight["ASE_version"] != cfg["software"]["ASE_version"]:
    preflight["status"] = "PREFLIGHT_FAIL_SOFTWARE_VERSION"
if start_free < preflight["minimum_start_free_bytes"]:
    preflight["status"] = "PREFLIGHT_FAIL_DISK_HEADROOM"
if preflight["memory_limit_bytes"] < cfg["resources"]["memory_limit_bytes"]:
    preflight["status"] = "PREFLIGHT_FAIL_MEMORY_LIMIT"
if preflight["memory_current_bytes"] >= MEMORY_REQUEST_STOP:
    preflight["status"] = "PREFLIGHT_FAIL_MEMORY_CURRENT"
if preflight["disk_stop_guard_bytes"] != DISK_GUARD:
    preflight["status"] = "PREFLIGHT_FAIL_DISK_GUARD_CONFIG"
if preflight["input_sha256"] != task["input_sha256"] or not preflight["PAW_hashes_match"]:
    preflight["status"] = "PREFLIGHT_FAIL_SHA"
if preflight["active_compute_processes"]:
    preflight["status"] = "PREFLIGHT_FAIL_CONCURRENT_COMPUTE"
if preflight["cpu_quota_cores"] != 4:
    preflight["status"] = "PREFLIGHT_FAIL_CPU_QUOTA"
preflight_path = (ROOT / f"PREFLIGHT_{args.task_id}.json" if args.preflight_only else
                  attempt / ("PREFLIGHT.json" if not args.resume else f"PREFLIGHT_{run_tag}.json"))
atomic_json(preflight_path, preflight)
if args.preflight_only:
    print(json.dumps(preflight, indent=2), flush=True)
    raise SystemExit(0 if preflight["status"] == "PREFLIGHT_PASS" else 1)
if preflight["status"] != "PREFLIGHT_PASS":
    raise SystemExit(f"Preflight failed: {preflight['status']}")

env = os.environ.copy()
env.update({
    "GPAW_MPI_BACKEND": "cgpaw", "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
    "GPAW_SETUP_PATH": str(setup), "TMPDIR": str(attempt), "TMP": str(attempt),
})
cmd = ["mpiexec", "--bind-to", "core", "--map-by", "core", "-n", "4",
       sys.executable, str(ROOT / "run_point.py"),
       "--task-id", args.task_id, "--attempt-dir", str(attempt)]
if args.resume:
    cmd.append("--resume")

atomic_json(attempt / f"SUPERVISOR_STARTED_{run_tag}.json", {
    "command": cmd, "started_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
    "max_walltime_seconds": MAX_WALLTIME, "soft_stop_seconds": SOFT_WALLTIME,
    "memory_request_stop_bytes": MEMORY_REQUEST_STOP, "memory_hard_guard_bytes": MEMORY_GUARD,
    "disk_request_stop_bytes": DISK_REQUEST_STOP, "disk_hard_guard_bytes": DISK_GUARD,
})
started = time.monotonic()
stop_reason = None
peak_memory = 0
peak_memory_sample = None
samples = 0
with (attempt / f"launcher_stdout_{run_tag}.log").open("wb") as out, (attempt / f"launcher_stderr_{run_tag}.log").open("wb") as err, (attempt / f"RESOURCE_PROFILE_{run_tag}.jsonl").open("w") as prof:
    proc = subprocess.Popen(cmd, cwd=REPO, env=env, stdin=subprocess.DEVNULL,
                            stdout=out, stderr=err, start_new_session=True)
    while proc.poll() is None:
        now = time.monotonic()
        current = read_int(Path("/sys/fs/cgroup/memory.current"))
        free = free_bytes(Path("/workspace"))
        samples += 1
        row = {"utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
               "elapsed_seconds": round(now - started, 1), "memory_current_bytes": current,
               "memory_stat": {}, "workspace_free_bytes": free}
        try:
            for line in Path("/sys/fs/cgroup/memory.stat").read_text().splitlines():
                k, v = line.split()
                if k in {"anon", "file", "inactive_file", "active_file", "kernel", "slab"}:
                    row["memory_stat"][k] = int(v)
        except OSError:
            pass
        prof.write(json.dumps(row, separators=(",", ":")) + "\n")
        prof.flush()
        if current > peak_memory:
            peak_memory = current
            peak_memory_sample = row
        if current >= MEMORY_REQUEST_STOP:
            stop_reason = "MEMORY_SOFT_GUARD"
        elif free <= DISK_REQUEST_STOP:
            stop_reason = "DISK_SOFT_GUARD"
        elif now - started >= SOFT_WALLTIME:
            stop_reason = "WALLTIME_4H50M_SAFE_STOP"
        if stop_reason:
            atomic_json(attempt / "STOP_REQUEST.json", {
                "reason": stop_reason, "requested_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                "memory_current_bytes": current, "workspace_free_bytes": free,
            })
            deadline = time.monotonic() + 180
            while proc.poll() is None and time.monotonic() < deadline:
                if current >= MEMORY_GUARD or free_bytes(Path("/workspace")) <= DISK_GUARD:
                    break
                time.sleep(2)
                current = read_int(Path("/sys/fs/cgroup/memory.current"))
            if proc.poll() is None:
                try:
                    os.killpg(proc.pid, signal.SIGTERM)
                    proc.wait(timeout=20)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()
            break
        if current >= MEMORY_GUARD or free_bytes(Path("/workspace")) <= DISK_GUARD:
            stop_reason = "HARD_RESOURCE_GUARD"
            atomic_json(attempt / "STOP_REQUEST.json", {
                "reason": stop_reason, "requested_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                "memory_current_bytes": current, "workspace_free_bytes": free,
            })
            try:
                os.killpg(proc.pid, signal.SIGTERM)
                proc.wait(timeout=20)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()
            break
        time.sleep(5)
    rc = proc.wait()

status_path = attempt / "RUN_STATUS.json"
if not status_path.exists():
    atomic_json(status_path, {
        "status": "WORKER_EXITED_WITHOUT_FINAL_STATUS",
        "exit_code": rc,
        "stop_reason": stop_reason or "PROCESS_EXIT",
        "scf_converged": False,
        "labels_written": (attempt / "RESULT.json").is_file(),
    })
status = json.loads(status_path.read_text())
summary = {
    "task_id": args.task_id,
    "status": status.get("status"),
    "scf_converged": status.get("scf_converged", False),
    "exit_code": rc,
    "stop_reason": status.get("stop_reason", stop_reason or "PROCESS_EXIT"),
    "elapsed_seconds": round(time.monotonic() - started, 3),
    "max_walltime_seconds": MAX_WALLTIME,
    "memory_current_peak_bytes": peak_memory,
    "memory_peak_sample": peak_memory_sample,
    "samples": samples,
    "workspace_free_bytes_at_end": free_bytes(Path("/workspace")),
    "checkpoint_status": json.loads((attempt / "CHECKPOINT_STATUS.json").read_text())
                          if (attempt / "CHECKPOINT_STATUS.json").is_file() else None,
    "labels_written": (attempt / "RESULT.json").is_file(),
    "ended_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
}
atomic_json(attempt / f"SUPERVISOR_SUMMARY_{run_tag}.json", summary)
print(json.dumps(summary, indent=2), flush=True)
raise SystemExit(0 if status.get("scf_converged") is True else 1)
