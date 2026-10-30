#!/usr/bin/env python3
"""Sequential 4-rank preflight launcher with a 30 GiB cgroup stop guard."""
from __future__ import annotations

import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
import re

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[3]
CGROUP = Path("/sys/fs/cgroup")
LIMIT = 34359738368
STOP = 32212254720
MIN_DISK = 5 * 1024**3
PYTHON = "/workspace/.venvs/gpaw-mpi/bin/python"
LOCK = Path("/workspace/.mace-v19-dft.lock")

def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def memory_current() -> int:
    return int((CGROUP / "memory.current").read_text())

def active_dft() -> list[str]:
    found = []
    for proc in Path("/proc").iterdir():
        if not proc.name.isdigit():
            continue
        try:
            cmd = (proc / "cmdline").read_bytes().decode(errors="ignore")
        except OSError:
            continue
        if "run_preflight_label.py" in cmd and int(proc.name) != os.getpid():
            found.append(proc.name)
    return found

def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: launch_preflight.py EXACT_REGISTERED_LABEL")
    manifest = json.loads((ROOT / "input_manifest.json").read_text())
    label = sys.argv[1]
    job = next(row for row in manifest["preflight_jobs"] if row["label"] == label)
    assert job["launch_enabled"] is True
    source = ROOT / job["input"]
    assert source.is_file() and sha(source) == job["input_sha256"]
    assert int((CGROUP / "memory.max").read_text()) == LIMIT
    assert memory_current() < 3 * 1024**3, "memory use before launch exceeds 3 GiB"
    assert shutil.disk_usage(ROOT).free >= MIN_DISK, "less than 5 GiB free disk"
    assert not active_dft(), "another preflight DFT process is active"
    target = ROOT / "calculations" / label
    scratch = ROOT / "calculations" / f".running_{label}"
    assert not target.exists() and not scratch.exists(), "existing calculation/archive found; preserve it"
    (ROOT / "calculations").mkdir(parents=True, exist_ok=True)
    launcher_log = ROOT / "calculations" / f"{label}.launcher.log"
    env = os.environ.copy()
    deps = "/workspace/.local/gpaw-mpi-deps/usr/lib/x86_64-linux-gnu"
    env.update({"GPAW_MPI_BACKEND": "cgpaw", "OMP_NUM_THREADS": "1",
                "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
                "NUMEXPR_NUM_THREADS": "1"})
    env.pop("V23_RESTART_CHECKPOINT", None)
    env.pop("V23_RESTART_ITERATION_OFFSET", None)
    env["LD_LIBRARY_PATH"] = ":".join([deps + "/openmpi/lib", deps,
        "/workspace/.local/gpaw-libs/usr/lib/x86_64-linux-gnu",
        "/workspace/.local/openblas/usr/lib/x86_64-linux-gnu/openblas-pthread",
        env.get("LD_LIBRARY_PATH", "")])
    checkpoint_dir = ROOT / "calculations" / "latest_checkpoint"
    checkpoint = checkpoint_dir / "latest_state.gpw"
    checkpoint_meta_path = checkpoint_dir / "latest_checkpoint.json"
    if checkpoint.is_file() and checkpoint_meta_path.is_file():
        try:
            checkpoint_meta = json.loads(checkpoint_meta_path.read_text())
            if (checkpoint_meta.get("label") == label
                    and checkpoint_meta.get("input_sha256") == job["input_sha256"]
                    and sha(checkpoint) == checkpoint_meta.get("checkpoint_sha256")):
                env["V23_RESTART_CHECKPOINT"] = str(checkpoint)
                env["V23_RESTART_ITERATION_OFFSET"] = str(
                    int(checkpoint_meta.get("last_reported_scf_iteration") or 0))
                print(f"Resuming {label} from verified latest GPAW checkpoint at SCF≈"
                      f"{env['V23_RESTART_ITERATION_OFFSET']}", flush=True)
        except (OSError, ValueError, json.JSONDecodeError):
            print("Latest checkpoint metadata/hash did not verify; starting from the registered input.",
                  flush=True)
    method = manifest["method_candidate"]
    info = subprocess.check_output([str(Path(PYTHON).with_name("gpaw")), "info"], env=env, text=True)
    assert re.search(r"MPI enabled\s+yes", info) and re.search(r"scalapack\s+yes", info), "GPAW requires four-rank MPI and ScaLAPACK"
    probe = "import hashlib,importlib.metadata as m,gpaw_data,pathlib; p=pathlib.Path(gpaw_data.datapath()); print(m.version('gpaw'),m.version('ase'),m.version('gpaw-data')); [print(n,hashlib.sha256((p/n).read_bytes()).hexdigest()) for n in " + repr(list(method["PAW_setup_sha256"])) + "]"
    package_info = subprocess.check_output([PYTHON, "-c", probe], env=env, text=True).splitlines()
    assert package_info[0] == "26.7.0 3.29.0 1.2.1", package_info
    setup_hashes = dict(line.split() for line in package_info[1:])
    assert setup_hashes == method["PAW_setup_sha256"], "Installed PAW setup digest differs from registered method"
    with LOCK.open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        assert not active_dft(), "a DFT runner started after preflight"
        start = time.time()
        start_mem = memory_current()
        peak = start_mem
        stopped = False
        with launcher_log.open("w") as log:
            child = subprocess.Popen(["mpirun", "--bind-to", "core", "--map-by", "core", "-n", "4",
                PYTHON, str(ROOT / "run_preflight_label.py"), label, str(scratch)],
                cwd=REPO, env=env, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            last_iteration = -1
            while child.poll() is None:
                time.sleep(10)
                now = memory_current()
                peak = max(peak, now)
                if now >= STOP:
                    stopped = True
                    os.killpg(child.pid, signal.SIGTERM)
                    try:
                        child.wait(timeout=30)
                    except subprocess.TimeoutExpired:
                        os.killpg(child.pid, signal.SIGKILL)
                    break
                progress_path = scratch / "progress.json"
                if progress_path.exists():
                    try:
                        progress = json.loads(progress_path.read_text())
                    except (OSError, json.JSONDecodeError):
                        progress = {}
                    iteration = int(progress.get("iteration", 0))
                    if iteration != last_iteration:
                        print(f"{label}: SCF {iteration}; memory={now/1024**3:.2f} GiB; peak={peak/1024**3:.2f} GiB", flush=True)
                        last_iteration = iteration
            code = child.wait()
        runtime = {"label": label, "return_code": code, "stopped_for_memory": stopped,
                   "elapsed_s": time.time() - start, "memory_start_bytes": start_mem,
                   "memory_high_water_observed_bytes": peak, "cgroup_limit_bytes": LIMIT,
                   "memory_stop_threshold_bytes": STOP}
        (ROOT / "calculations" / f"{label}.runtime.json").write_text(json.dumps(runtime, indent=2) + "\n")
        if stopped or code != 0:
            if scratch.exists():
                # Preserve the latest GPAW state so a later invocation can resume it.
                (scratch / "state.next.gpw").unlink(missing_ok=True)
                failed_root = ROOT / "calculations" / "failed_attempts"
                failed_root.mkdir(exist_ok=True)
                failed = failed_root / f"{label}__attempt_{int(start)}"
                if failed.exists():
                    raise SystemExit(f"failed-attempt archive already exists: {failed}")
                scratch.rename(failed)
                shutil.copy2(ROOT / "calculations" / f"{label}.runtime.json", failed / "runtime.json")
                shutil.copy2(launcher_log, failed / "launcher.log")
            raise SystemExit(f"DFT preflight failed (exit={code}, memory_stop={stopped}); no automatic retry")
    # The archive contains only auditable text/label artifacts, never .gpw state.
    (scratch / "state.gpw").unlink(missing_ok=True)
    if not (scratch / "summary.json").is_file():
        raise SystemExit("runner exited without summary; preserve scratch for review")
    scratch.rename(target)
    inventory = "".join(f"{sha(p)}  {p.name}\n" for p in sorted(target.iterdir())
                        if p.is_file() and p.name != "SHA256SUMS.txt")
    (target / "SHA256SUMS.txt").write_text(inventory)
    print(json.dumps({"status": "RUN_COMPLETE_TEXT_ARCHIVE", "label": label,
                      "summary": json.loads((target / "summary.json").read_text()),
                      "runtime": runtime}, indent=2), flush=True)

if __name__ == "__main__":
    main()
