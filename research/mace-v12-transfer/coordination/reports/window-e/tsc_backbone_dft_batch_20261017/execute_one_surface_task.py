#!/usr/bin/env python3
"""Preflight, run one fixed surface DFT point, archive it, and verify fresh-fetch."""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import io
import json
import os
import shutil
import signal
import subprocess
import sys
import tarfile
import time
from collections import Counter
from pathlib import Path

import ase
import gpaw
import gpaw_data
import numpy as np
from ase.io import read


BASE = Path(__file__).resolve().parent
REPO = Path(subprocess.check_output(["git", "rev-parse", "--show-toplevel"], cwd=BASE, text=True).strip())
BRANCH = "cloud-e-b1-dft"
CONFIG_PATH = BASE / "SURFACE_DFT_CONFIG.json"
LOCK_PATH = BASE / "EXECUTION_LOCK_SURFACE_TRAIN_20261017.json"
WORKER_PATH = BASE / "run_surface_singlepoint.py"
CONFIG_SHA = "e9a4e27d23c4646b873e77de7acbb79844422c8708c9f24d209827f1849312d2"
REQUIRED_START_BYTES = 17179869184  # 10 GiB scratch + 6 GiB headroom
MIN_RUNTIME_HEADROOM_BYTES = 6442450944  # 6 GiB
MAX_CGROUP_CURRENT_BYTES = 30064771072  # 28 GiB conservative stop threshold
MAX_WALLTIME_SECONDS = 18000
RETRY_TASK_ID = "TSC100_TI0001_COD9009647_2x2x2_REF"
RETRY_ATTEMPT_ID = "attempt-0002"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def run(cmd, *, cwd=REPO, env=None, check=True, capture=False):
    return subprocess.run(cmd, cwd=cwd, env=env, check=check, text=True,
                          stdout=subprocess.PIPE if capture else None,
                          stderr=subprocess.PIPE if capture else None)


def available_bytes(path: Path) -> int:
    st = os.statvfs(path)
    return st.f_bavail * st.f_frsize


def task_output_path(task_id: str, attempt_id: str | None) -> Path:
    task_root = BASE / "results" / task_id
    if attempt_id is None:
        return task_root
    if task_id != RETRY_TASK_ID or attempt_id != RETRY_ATTEMPT_ID:
        raise ValueError(
            f"Only {RETRY_ATTEMPT_ID} is permitted for {RETRY_TASK_ID}"
        )
    return task_root / attempt_id


def read_int(path: str) -> int:
    value = Path(path).read_text().strip()
    if value == "max":
        return 2**63 - 1
    return int(value)


def active_compute_processes():
    found = []
    for p in Path("/proc").iterdir():
        if not p.name.isdigit():
            continue
        try:
            status = (p / "status").read_text(errors="replace")
            comm = next(line.split("\t", 1)[1].strip() for line in status.splitlines() if line.startswith("Name:"))
            state = next(line.split()[1] for line in status.splitlines() if line.startswith("State:"))
            if state != "Z" and comm in {"gpaw", "mpiexec", "mpirun"}:
                found.append({"pid": int(p.name), "comm": comm, "state": state})
        except (OSError, StopIteration, ValueError):
            continue
    return found


def check_mpi(env, expected_versions):
    code = ("from gpaw.mpi import world; import ase,gpaw,socket; "
            "print(f'{world.rank} {world.size} {socket.gethostname()} {gpaw.__version__} {ase.__version__}')")
    p = run([shutil.which("mpiexec") or "mpiexec", "-n", "4", sys.executable, "-c", code],
            cwd=REPO, env=env, capture=True)
    rows = [line.split() for line in p.stdout.splitlines() if line.strip()]
    if len(rows) != 4 or {int(r[0]) for r in rows} != {0, 1, 2, 3}:
        raise RuntimeError(f"MPI probe did not return ranks 0..3: {p.stdout!r}")
    if any(int(r[1]) != 4 or r[3] != expected_versions["GPAW_version"]
           or r[4] != expected_versions["ASE_version"] for r in rows):
        raise RuntimeError(f"MPI ranks disagree on size or software versions: {rows}")
    if len({r[2] for r in rows}) != 1:
        raise RuntimeError(f"MPI ranks span different hosts: {rows}")
    return rows


def verify_manifest(base: Path):
    p = subprocess.run(["sha256sum", "-c", "SHA256SUMS.txt"], cwd=base,
                       text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if p.returncode:
        raise RuntimeError("Package SHA check failed:\n" + p.stdout)
    return p.stdout


def preflight(task_id: str, attempt_id: str | None,
              confirm_window: bool, expected_lock_sha: str):
    cfg = json.loads(CONFIG_PATH.read_text())
    lock = json.loads(LOCK_PATH.read_text())
    checks = {}
    checks["lock_authorized"] = lock.get("execution_authorized") is True
    checks["user_confirmed_five_hour_window"] = bool(confirm_window)
    checks["config_sha256"] = sha256(CONFIG_PATH) == CONFIG_SHA == lock["config_sha256"]
    checks["lock_sha256"] = sha256(LOCK_PATH) == expected_lock_sha
    checks["worker_sha256"] = sha256(WORKER_PATH) == lock["runner_sha256"]
    checks["supervisor_sha256"] = sha256(Path(__file__).resolve()) == lock["supervisor_sha256"]
    checks["package_checksums"] = bool(verify_manifest(BASE))
    checks["correct_git_branch"] = run(["git", "branch", "--show-current"], capture=True).stdout.strip() == BRANCH
    run(["git", "fetch", "--no-tags", "origin",
         f"+refs/heads/{BRANCH}:refs/remotes/origin/{BRANCH}"])
    local_head = run(["git", "rev-parse", "HEAD"], capture=True).stdout.strip()
    remote_head = run(["git", "rev-parse", f"refs/remotes/origin/{BRANCH}"], capture=True).stdout.strip()
    checks["branch_fresh_and_clean"] = local_head == remote_head and not run(
        ["git", "status", "--porcelain", "--untracked-files=no"], capture=True
    ).stdout.strip() and run(["git", "diff", "--cached", "--quiet"], check=False).returncode == 0
    checks["input_commit_is_ancestor"] = run(
        ["git", "merge-base", "--is-ancestor", lock["input_package_commit"], "HEAD"], check=False
    ).returncode == 0
    run(["git", "fetch", "--no-tags", "origin", lock["C_review_commit"]])
    c_review_path = "research/mace-v12-transfer/coordination/reports/window-c/tsc_surface_input_review_20261017/C_INPUT_REVIEW.md"
    c_decision_path = "research/mace-v12-transfer/coordination/reports/window-c/tsc_surface_input_review_20261017/C_INPUT_DECISION.json"
    c_review = subprocess.check_output(["git", "show", f"{lock['C_review_commit']}:{c_review_path}"], cwd=REPO)
    c_decision = subprocess.check_output(["git", "show", f"{lock['C_review_commit']}:{c_decision_path}"], cwd=REPO)
    checks["C_review_sha256"] = sha256_bytes(c_review) == lock["C_review_sha256"]
    checks["C_decision_sha256"] = sha256_bytes(c_decision) == lock["C_decision_sha256"]
    task = next(t for t in cfg["tasks"] if t["structure_id"] == task_id)
    input_path = BASE / task["input_path"]
    checks["input_sha256"] = sha256(input_path) == task["input_sha256"]
    atoms = read(input_path, format="extxyz")
    symbols = Counter(atoms.get_chemical_symbols())
    h = hashlib.sha256()
    h.update(" ".join(atoms.get_chemical_symbols()).encode())
    h.update(np.asarray(atoms.cell, dtype="<f8").tobytes(order="C"))
    h.update(np.asarray(np.round(atoms.positions, 8), dtype="<f8").tobytes(order="C"))
    h.update(np.asarray(atoms.pbc, dtype=np.uint8).tobytes())
    h.update(np.asarray(atoms.arrays["source_atom_id"], dtype="<i8").tobytes())
    geometry_digest = h.hexdigest()
    checks["geometry"] = (
        len(atoms) == task["atom_count"] == 100
        and dict(symbols) == task["composition"]
        and geometry_digest == task["geometry_sha256"]
        and atoms.pbc.tolist() == [True, True, True]
        and np.array_equal(atoms.arrays.get("source_atom_id"), np.arange(1, 101))
        and "energy" not in atoms.info and "free_energy" not in atoms.info
        and "forces" not in atoms.arrays and "stress" not in atoms.info
    )
    checks["GPAW_version"] = gpaw.__version__ == cfg["software"]["GPAW_version"]
    checks["ASE_version"] = ase.__version__ == cfg["software"]["ASE_version"]
    checks["gpaw_data_version"] = importlib.metadata.version("gpaw-data") == cfg["software"]["gpaw_data_version"]
    setup_path = Path(gpaw_data.datapath()).resolve()
    os.environ["GPAW_SETUP_PATH"] = str(setup_path)
    gpaw_executable = shutil.which("gpaw") or str(Path(sys.executable).with_name("gpaw"))
    checks["setup_path"] = (
        setup_path.is_dir()
        and shutil.which("mpiexec") is not None
        and Path(gpaw_executable).is_file()
        and os.access(gpaw_executable, os.X_OK)
    )
    for filename, expected in cfg["paw_sha256"].items():
        checks[f"PAW_{filename}"] = sha256(setup_path / filename) == expected
    cpu_max = Path("/sys/fs/cgroup/cpu.max").read_text().split()
    checks["cpu_quota_4"] = len(cpu_max) == 2 and int(cpu_max[0]) / int(cpu_max[1]) >= 4
    checks["memory_limit_32GiB"] = read_int("/sys/fs/cgroup/memory.max") >= cfg["resources"]["minimum_memory_limit_bytes"]
    checks["memory_below_stop_threshold"] = read_int("/sys/fs/cgroup/memory.current") < MAX_CGROUP_CURRENT_BYTES
    checks["resource_thresholds_match_lock"] = (
        cfg["resources"]["minimum_start_free_bytes"] == REQUIRED_START_BYTES
        and cfg["resources"]["minimum_runtime_headroom_bytes"] == MIN_RUNTIME_HEADROOM_BYTES
        and cfg["resources"]["maximum_cgroup_memory_current_bytes"] == MAX_CGROUP_CURRENT_BYTES
        and cfg["resources"]["per_structure_walltime_seconds"] == MAX_WALLTIME_SECONDS
        and cfg["resources"]["automatic_retry"] is False
    )
    free = available_bytes(Path("/workspace"))
    output_path = task_output_path(task_id, attempt_id)
    retry_policy = lock.get("manual_retry_policy", {})
    if attempt_id is not None:
        previous_attempt = BASE / retry_policy.get("preserve_existing_attempt_path", "INVALID")
        previous_manifest = previous_attempt / "SHA256SUMS.txt"
        previous_check = subprocess.run(
            ["sha256sum", "-c", "SHA256SUMS.txt"], cwd=previous_attempt,
            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
        ) if previous_attempt.is_dir() else None
        checks["manual_attempt_authorized_by_lock"] = (
            retry_policy.get("task_id") == task_id
            and retry_policy.get("attempt_id") == attempt_id
            and retry_policy.get("automatic_retry") is False
        )
        checks["retry_path_matches_lock"] = (
            output_path.relative_to(BASE).as_posix() == retry_policy.get("new_attempt_path")
        )
        checks["previous_attempt_manifest_sha256"] = (
            previous_manifest.is_file()
            and sha256(previous_manifest) == retry_policy.get("preserved_attempt_manifest_sha256")
        )
        checks["previous_attempt_archive_hashes"] = (
            previous_check is not None and previous_check.returncode == 0
        )
        checks["previous_attempt_commit_in_history"] = run(
            ["git", "merge-base", "--is-ancestor",
             retry_policy.get("previous_attempt_commit", "INVALID"), "HEAD"], check=False
        ).returncode == 0
    checks["initial_scratch_and_headroom"] = free >= REQUIRED_START_BYTES
    checks["output_path_is_new"] = not output_path.exists()
    checks["no_active_DFT_process"] = not active_compute_processes()
    env = os.environ.copy()
    env["GPAW_MPI_BACKEND"] = "cgpaw"
    env["OMP_NUM_THREADS"] = "1"
    env["GPAW_SETUP_PATH"] = str(setup_path)
    mpi_rows = check_mpi(env, cfg["software"])
    checks["MPI_4_rank_consistency"] = len(mpi_rows) == 4
    return checks, {"task": task, "config": cfg, "lock": lock,
                    "output_path": str(output_path), "attempt_id": attempt_id,
                    "workspace_free_bytes": free, "memory_current_bytes": read_int("/sys/fs/cgroup/memory.current"),
                    "memory_limit_bytes": read_int("/sys/fs/cgroup/memory.max"),
                    "cpu_max": Path("/sys/fs/cgroup/cpu.max").read_text().strip(),
                    "hostname": os.uname().nodename,
                    "python_version": sys.version,
                    "GPAW_version": gpaw.__version__, "ASE_version": ase.__version__,
                    "gpaw_data_version": importlib.metadata.version("gpaw-data"),
                    "GPAW_SETUP_PATH": str(setup_path),
                    "MPI_probe": mpi_rows}


def terminate_group(proc):
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        proc.wait(timeout=30)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.wait()


def sha_tree(root: Path):
    lines = []
    for path in sorted(p for p in root.rglob("*") if p.is_file() and p.name != "SHA256SUMS.txt"):
        rel = path.relative_to(root).as_posix()
        lines.append(f"{sha256(path)}  ./{rel}")
    (root / "SHA256SUMS.txt").write_text("\n".join(lines) + "\n")


def verify_remote_dir(remote_ref: str, repo_relpath: str):
    tar_bytes = subprocess.run(["git", "archive", "--format=tar", remote_ref, repo_relpath],
                               cwd=REPO, check=True, stdout=subprocess.PIPE).stdout
    with tarfile.open(fileobj=io.BytesIO(tar_bytes), mode="r:") as tf:
        prefix = repo_relpath.rstrip("/") + "/"
        manifest = tf.extractfile(prefix + "SHA256SUMS.txt").read().decode()
        for line in manifest.splitlines():
            if not line.strip():
                continue
            expected, rel = line.split(None, 1)
            rel = rel.lstrip("*")
            if rel.startswith("./"):
                rel = rel[2:]
            actual = sha256_bytes(tf.extractfile(prefix + rel).read())
            if actual != expected:
                raise RuntimeError(f"fresh-fetch SHA mismatch for {repo_relpath}/{rel}")
    return len([x for x in manifest.splitlines() if x.strip()])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--task-id", required=True)
    parser.add_argument("--attempt-id", choices=[RETRY_ATTEMPT_ID])
    parser.add_argument("--lock-sha256", required=True)
    parser.add_argument("--confirm-5h-window", action="store_true")
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()
    checks, details = preflight(args.task_id, args.attempt_id,
                                args.confirm_5h_window, args.lock_sha256)
    report = {"status": "PREFLIGHT_PASS" if all(checks.values()) else "PREFLIGHT_HOLD",
              "checks": checks, "details": details,
              "minimum_start_free_bytes": REQUIRED_START_BYTES,
              "minimum_runtime_headroom_bytes": MIN_RUNTIME_HEADROOM_BYTES,
              "maximum_cgroup_current_bytes": MAX_CGROUP_CURRENT_BYTES,
              "per_structure_walltime_seconds": MAX_WALLTIME_SECONDS,
              "checkpoint_policy": "disabled; preserve failure logs, no restart or automatic retry"}
    print(json.dumps(report, indent=2))
    if not all(checks.values()):
        return 2
    if args.preflight_only:
        return 0

    task = details["task"]
    cfg = details["config"]
    out = Path(details["output_path"])
    out.mkdir(parents=True, exist_ok=False)
    input_path = BASE / task["input_path"]
    shutil.copy2(input_path, out / "input.extxyz")
    (out / "RUN_MANIFEST.json").write_text(json.dumps({
        "structure_id": args.task_id,
        "attempt_id": args.attempt_id or "attempt-0001",
        "role": details["lock"]["role"],
        "input_sha256": task["input_sha256"], "config_sha256": CONFIG_SHA,
        "execution_lock_sha256": args.lock_sha256,
        "source_input_commit": details["lock"]["input_package_commit"],
        "C_review_commit": details["lock"]["C_review_commit"],
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "walltime_limit_seconds": MAX_WALLTIME_SECONDS,
        "checkpoint_enabled": False, "automatic_retry": False,
        "environment": {
            "hostname": details["hostname"], "python_version": details["python_version"],
            "GPAW_version": details["GPAW_version"], "ASE_version": details["ASE_version"],
            "gpaw_data_version": details["gpaw_data_version"],
            "GPAW_SETUP_PATH": details["GPAW_SETUP_PATH"],
            "GPAW_MPI_BACKEND": "cgpaw", "MPI_ranks": 4, "OMP_NUM_THREADS": 1,
            "cpu_max": details["cpu_max"],
            "memory_limit_bytes": details["memory_limit_bytes"],
            "memory_current_bytes_at_start": details["memory_current_bytes"],
            "workspace_free_bytes_at_start": details["workspace_free_bytes"],
            "paw_sha256": cfg["paw_sha256"],
        },
    }, indent=2) + "\n")
    scratch_dir = out / "scratch"
    scratch_dir.mkdir()
    env = os.environ.copy()
    env.update({"GPAW_MPI_BACKEND": "cgpaw", "OMP_NUM_THREADS": "1",
                "GPAW_SETUP_PATH": str(Path(gpaw_data.datapath()).resolve()),
                "TMPDIR": str(scratch_dir), "TMP": str(scratch_dir)})
    cmd = [shutil.which("mpiexec") or "mpiexec", "-n", "4",
           shutil.which("gpaw") or str(Path(sys.executable).with_name("gpaw")), "python",
           "--", str(WORKER_PATH), "--task-id", args.task_id,
           "--input", str(out / "input.extxyz"), "--run-dir", str(out),
           "--config", str(CONFIG_PATH), "--config-sha256", CONFIG_SHA,
           "--execution-lock", str(LOCK_PATH), "--execution-lock-sha256", args.lock_sha256]
    started = time.monotonic()
    stop_reason = None
    min_free_bytes = details["workspace_free_bytes"]
    max_memory_bytes = details["memory_current_bytes"]
    with (out / "stdout.log").open("wb") as stdout, (out / "stderr.log").open("wb") as stderr:
        proc = subprocess.Popen(cmd, cwd=REPO, env=env, stdout=stdout, stderr=stderr,
                                start_new_session=True)
        while proc.poll() is None:
            elapsed = time.monotonic() - started
            free = available_bytes(Path("/workspace"))
            mem = read_int("/sys/fs/cgroup/memory.current")
            min_free_bytes = min(min_free_bytes, free)
            max_memory_bytes = max(max_memory_bytes, mem)
            if elapsed >= MAX_WALLTIME_SECONDS:
                stop_reason = "WALLTIME_LIMIT_18000_SECONDS"
            elif free < MIN_RUNTIME_HEADROOM_BYTES:
                stop_reason = "DISK_HEADROOM_BELOW_6_GIB"
            elif mem >= MAX_CGROUP_CURRENT_BYTES:
                stop_reason = "CGROUP_MEMORY_AT_OR_ABOVE_28_GIB"
            if stop_reason:
                terminate_group(proc)
                break
            time.sleep(2)
        exit_code = proc.wait()
    result_exists = (out / "RESULT.json").exists()
    success = exit_code == 0 and result_exists and stop_reason is None
    status = {
        "structure_id": args.task_id,
        "status": "SUCCESS" if success else "FAILED",
        "exit_code": exit_code,
        "stop_reason": stop_reason,
        "elapsed_seconds": time.monotonic() - started,
        "result_written_only_after_scf_convergence": result_exists,
        "checkpoint_written": False,
        "automatic_retry": False,
        "minimum_workspace_free_bytes_observed": min_free_bytes,
        "maximum_cgroup_memory_current_bytes_observed": max_memory_bytes,
    }
    (out / "RUN_STATUS.json").write_text(json.dumps(status, indent=2) + "\n")
    sha_tree(out)

    rel = out.relative_to(REPO).as_posix()
    run(["git", "add", "--", rel])
    run(["git", "commit", "--only", "-m", f"[E] Archive surface DFT {args.task_id}", "--", rel])
    run(["git", "push", "origin", f"HEAD:refs/heads/{BRANCH}"])
    run(["git", "fetch", "--no-tags", "origin", f"+refs/heads/{BRANCH}:refs/remotes/origin/{BRANCH}"])
    head = run(["git", "rev-parse", "HEAD"], capture=True).stdout.strip()
    remote = run(["git", "rev-parse", f"refs/remotes/origin/{BRANCH}"], capture=True).stdout.strip()
    if head != remote:
        raise RuntimeError(f"Remote branch tip mismatch after push: local={head} remote={remote}")
    checked = verify_remote_dir(f"refs/remotes/origin/{BRANCH}", rel)
    print(json.dumps({"GITHUB_COMMIT_SHA": head, "REMOTE_SHA256_CHECKED_FILES": checked,
                      "TASK_STATUS": status["status"], "RESULT_PATH": rel}, indent=2))
    return 0 if success else 3


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"PREFLIGHT_OR_EXECUTION_ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        raise
