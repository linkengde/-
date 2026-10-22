#!/usr/bin/env python3
"""Run exactly one registered V19 independent-development label."""
import fcntl
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[3]
PYTHON = "/workspace/.venvs/gpaw-mpi/bin/python"

if len(sys.argv) != 2:
    raise SystemExit("usage: run_queue.py LABEL")
label = sys.argv[1]
sys.path.insert(0, str(ROOT))
from entry import load

manifest, record = load(label, launch=True)
lock = open("/workspace/.mace-v19-dft.lock", "a")
fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
for proc in Path("/proc").iterdir():
    if not proc.name.isdigit():
        continue
    try:
        names = [Path(arg.decode(errors="ignore")).name for arg in (proc / "cmdline").read_bytes().split(b"\0")]
    except OSError:
        continue
    if any(name in {"run_pw_reference.py", "run_pw_reference_draft.py"} for name in names):
        raise SystemExit("Another GPAW reference runner is active; do not duplicate")

out = Path("/workspace/mace_v19_independent_dev") / label
archive = REPO / record["result_directory"] / "calculations" / label
assert not out.exists() and not archive.exists(), "Existing run/archive: preserve and inspect"
assert shutil.disk_usage(ROOT).free >= manifest["minimum_free_disk_before_launch_bytes"], "Insufficient free disk"
env = os.environ.copy()
deps = "/workspace/.local/gpaw-mpi-deps/usr/lib/x86_64-linux-gnu"
env.update(
    GPAW_MPI_BACKEND="cgpaw",
    OMP_NUM_THREADS="1",
    OPENBLAS_NUM_THREADS="1",
    MKL_NUM_THREADS="1",
    NUMEXPR_NUM_THREADS="1",
    XDG_CACHE_HOME="/workspace/.cache",
    MPLCONFIGDIR="/workspace/.cache/matplotlib",
    FC_CACHEDIR="/workspace/.cache/fontconfig",
    OMPI_ALLOW_RUN_AS_ROOT="1",
    OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1",
)
env["LD_LIBRARY_PATH"] = ":".join(
    [deps + "/openmpi/lib", deps, "/workspace/.local/gpaw-libs/usr/lib/x86_64-linux-gnu", "/workspace/.local/openblas/usr/lib/x86_64-linux-gnu/openblas-pthread", env.get("LD_LIBRARY_PATH", "")]
)
info = subprocess.check_output([str(Path(PYTHON).with_name("gpaw")), "info"], env=env, text=True)
assert re.search(r"MPI enabled\s+yes", info) and re.search(r"scalapack\s+yes", info), info
mpi = subprocess.run(
    ["mpirun", "--bind-to", "core", "--map-by", "core", "-n", "4", PYTHON, "-c", "from gpaw.mpi import world; import _gpaw; assert world.size == 4 and hasattr(_gpaw, 'scalapack_diagonalize_dc'); print('4-rank GPAW+ScaLAPACK PASS', world.size)"],
    env=env, text=True, capture_output=True, check=True,
)
print(mpi.stdout, end="", flush=True)

sys.path.insert(0, str(REPO / "research/mace-v12-transfer/coordination"))
import sync_tasks
sync_tasks.progress("window-a", label, "running", 0, None, "Starting registered fixed-geometry independent-development PW-PBE label on four GPAW MPI ranks")
out.parent.mkdir(parents=True, exist_ok=True)
launcher = ROOT / "calculations" / f"{label}.launcher.log"
launcher.parent.mkdir(parents=True, exist_ok=True)
with launcher.open("w") as log:
    child = subprocess.Popen(["mpirun", "--bind-to", "core", "--map-by", "core", "-n", "4", PYTHON, str(ROOT / "run_pw_reference.py"), label, str((REPO / record["input"]).resolve()), str(out.resolve())], env=env, stdout=log, stderr=subprocess.STDOUT)
    last = 0
    while child.poll() is None:
        time.sleep(45)
        progress_file = out / "progress.json"
        if progress_file.exists():
            try:
                iteration = int(json.loads(progress_file.read_text()).get("iteration", 0))
            except json.JSONDecodeError:
                continue
            if iteration >= last + 10:
                sync_tasks.progress("window-a", label, "running", iteration, None, "Four-rank GPAW SCF active; checkpoint retained locally")
                last = iteration
    if child.returncode:
        sync_tasks.progress("window-a", label, "failed", last, None, f"GPAW runner exit {child.returncode}; preserve run and checkpoint")
        raise SystemExit(child.returncode)

assert archive.is_dir()
files = [p for p in archive.iterdir() if p.is_file()]
assert all(p.suffix != ".gpw" and p.stat().st_size < 100_000_000 for p in files)
(archive / "SHA256SUMS.txt").write_text("".join(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n" for p in sorted(files)))
verified = subprocess.check_output([PYTHON, str(ROOT / "verify_result.py"), label, str(archive)], env=env, text=True)
assert json.loads(verified)["status"] == "PASS"
files = [p for p in archive.iterdir() if p.is_file() and p.name not in {"SHA256SUMS.txt"}]
(archive / "SHA256SUMS.txt").write_text("".join(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n" for p in sorted(files)))
subprocess.run([PYTHON, str(ROOT / "verify_result.py"), label, str(archive)], env=env, check=True)
subprocess.run([sys.executable, str(REPO / "research/mace-v12-transfer/coordination/sync_tasks.py"), "publish", "window-a", "--label", label], check=True)
shutil.rmtree(out)
print(f"Verified and uploaded {label}; removed local checkpoint/run cache.", flush=True)
