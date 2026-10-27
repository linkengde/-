#!/usr/bin/env python3
"""Safe sequential launcher with a hard 30 GiB memory-current stop threshold."""
from __future__ import annotations

import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[3]
CGROUP = Path("/sys/fs/cgroup")
MEMORY_LIMIT = 34359738368
MEMORY_STOP = 32212254720  # stop at 30 GiB, leaving 2 GiB below the cgroup limit
MINIMUM_DISK = 3221225472


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def memory_current() -> int:
    return int((CGROUP / "memory.current").read_text())


def active_dft() -> list[str]:
    matches = []
    for proc in Path("/proc").iterdir():
        if not proc.name.isdigit():
            continue
        try:
            cmd = (proc / "cmdline").read_bytes().decode(errors="ignore")
        except OSError:
            continue
        if "run_dft.py" in cmd and str(os.getpid()) not in proc.name:
            matches.append(proc.name)
    return matches


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: launch_one.py LABEL")
    label = sys.argv[1]
    manifest_path = ROOT / "input_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    record = next(item for item in manifest["records"] if item["label"] == label)
    assert manifest["launch_enabled"] and record["launch_enabled"] and record["status"] == "queued"
    source = ROOT / record["input"]
    assert source.is_file() and sha(source) == record["input_sha256"]
    assert int((CGROUP / "memory.max").read_text()) == MEMORY_LIMIT
    assert memory_current() < 3 * 1024**3, "memory usage is too high before launch"
    assert shutil.disk_usage(ROOT).free >= MINIMUM_DISK
    assert not active_dft(), "another GPAW calculation is active"

    setup_root = Path("/workspace/.venvs/gpaw-mpi/lib/python3.12/site-packages/gpaw_data")
    try:
        import gpaw_data
        setup_root = Path(gpaw_data.datapath())
    except Exception as exc:
        raise SystemExit(f"cannot resolve installed GPAW setups: {exc}")
    for name, digest in manifest["method"].get("PAW_setup_sha256", {}).items():
        path = setup_root / name
        assert path.is_file() and sha(path) == digest, f"PAW setup mismatch: {name}"

    target = ROOT / "calculations" / label
    scratch = ROOT / "calculations" / f".running_{label}"
    assert not target.exists() and not scratch.exists(), "existing run or archive found; preserve it"
    launcher_log = ROOT / "calculations" / f"{label}.launcher.log"
    ROOT.joinpath("calculations").mkdir(parents=True, exist_ok=True)
    lock_path = Path("/workspace/.mace-v19-dft.lock")
    with lock_path.open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        assert not active_dft(), "another GPAW calculation started after preflight"
        env = os.environ.copy()
        dft_deps = "/workspace/.local/gpaw-mpi-deps/usr/lib/x86_64-linux-gnu"
        dft_libs = "/workspace/.local/gpaw-libs/usr/lib/x86_64-linux-gnu"
        openblas_dir = "/workspace/.local/openblas/usr/lib/x86_64-linux-gnu/openblas-pthread"
        env.update({
            "GPAW_MPI_BACKEND": "cgpaw",
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
            "LD_LIBRARY_PATH": ":".join([dft_deps + "/openmpi/lib", dft_deps, dft_libs, openblas_dir]),
        })
        py = "/workspace/.venvs/gpaw-mpi/bin/python"
        command = ["mpirun", "--bind-to", "core", "--map-by", "core", "-n", "4", py, str(ROOT / "run_dft.py"), label, str(scratch)]
        start = time.time()
        memory_start = memory_current()
        max_mem = memory_start
        with launcher_log.open("w") as log:
            child = subprocess.Popen(command, cwd=REPO, env=env, stdout=log, stderr=subprocess.STDOUT)
            last_progress = -1
            stopped_for_memory = False
            while child.poll() is None:
                time.sleep(10)
                mem = memory_current()
                max_mem = max(max_mem, mem)
                if mem >= MEMORY_STOP:
                    child.terminate()
                    stopped_for_memory = True
                    try:
                        child.wait(timeout=30)
                    except subprocess.TimeoutExpired:
                        child.kill()
                    break
                progress_path = scratch / "progress.json"
                if progress_path.exists():
                    try:
                        progress = json.loads(progress_path.read_text())
                    except json.JSONDecodeError:
                        progress = {}
                    iteration = progress.get("iteration", 0)
                    if iteration != last_progress:
                        print(f"{label}: SCF {iteration}; cgroup memory {mem/1024**3:.2f} GiB; peak observed {max_mem/1024**3:.2f} GiB", flush=True)
                        last_progress = iteration
        return_code = child.wait()
        runtime = {
            "label": label,
            "return_code": return_code,
            "stopped_for_memory": stopped_for_memory,
            "elapsed_s": time.time() - start,
            "memory_start_bytes": memory_start,
            "memory_high_water_observed_bytes": max_mem,
            "cgroup_limit_bytes": MEMORY_LIMIT,
            "memory_stop_threshold_bytes": MEMORY_STOP,
        }
        (ROOT / "calculations" / f"{label}.runtime.json").write_text(json.dumps(runtime, indent=2) + "\n")
        if stopped_for_memory or return_code != 0:
            if scratch.exists():
                failed = ROOT / "calculations" / "failed_attempts" / label
                failed.parent.mkdir(exist_ok=True)
                if failed.exists():
                    raise SystemExit(f"failed-attempt archive already exists; preserve both: {failed}")
                scratch.rename(failed)
            raise SystemExit(f"calculation failed (return {return_code}, memory stop={stopped_for_memory}); no automatic retry")

    scratch.rename(target)
    verification = subprocess.check_output([py, str(ROOT / "verify_result.py"), label], cwd=REPO, text=True, env=env)
    verified = json.loads(verification)
    assert verified["status"] == "PASS"
    for path in target.glob("*.gpw"):
        path.unlink()
    (target / "verification.json").write_text(json.dumps(verified, indent=2) + "\n")
    (target / "SHA256SUMS.txt").write_text("".join(
        f"{sha(path)}  {path.name}\n" for path in sorted(target.iterdir()) if path.is_file() and path.name != "SHA256SUMS.txt"
    ))
    print(json.dumps(verified, indent=2), flush=True)


if __name__ == "__main__":
    main()
