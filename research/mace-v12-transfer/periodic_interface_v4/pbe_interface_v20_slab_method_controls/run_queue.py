"""Run one exact, A-owned slab k-mesh or smearing control."""
import fcntl
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

from entry import ROOT, REPO, load, sha

if len(sys.argv) != 2:
    raise SystemExit("usage: run_queue.py EXACT_LABEL")
label = sys.argv[1]
manifest, record = load(label, launch=True)

lock = open("/workspace/.mace-v19-dft.lock", "a")
fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
for proc in Path("/proc").iterdir():
    if not proc.name.isdigit():
        continue
    try:
        names = [Path(x.decode(errors="ignore")).name for x in (proc / "cmdline").read_bytes().split(b"\0")]
    except OSError:
        continue
    assert not any(n in ("run_pw_reference.py", "run_pw_reference_draft.py") for n in names), "Active GPAW runner: no duplicate"

out = Path("/workspace/mace_v20_slab_method_controls") / label
archive = ROOT / "calculations" / label
assert not out.exists() and not archive.exists(), "Existing run/archive: preserve and review"

sys.path.insert(0, str(REPO / "research/mace-v12-transfer/coordination"))
import sync_tasks as sync

sync.claim(record["owner_task"])
subprocess.run(
    [sys.executable, str(sync.HERE / "sync_tasks.py"), "progress", record["owner_task"], "--job", label,
     "--state", "running", "--iteration", "0", "--note", "Starting registered slab k-mesh/smearing control; four MPI ranks"],
    check=True,
)

env = os.environ.copy()
deps = "/workspace/.local/gpaw-mpi-deps/usr/lib/x86_64-linux-gnu"
env.update(
    GPAW_MPI_BACKEND="cgpaw", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
    NUMEXPR_NUM_THREADS="1", XDG_CACHE_HOME="/workspace/.cache", MPLCONFIGDIR="/workspace/.cache/matplotlib",
    FC_CACHEDIR="/workspace/.cache/fontconfig",
)
env["LD_LIBRARY_PATH"] = ":".join([
    deps + "/openmpi/lib", deps, "/workspace/.local/gpaw-libs/usr/lib/x86_64-linux-gnu",
    "/workspace/.local/openblas/usr/lib/x86_64-linux-gnu/openblas-pthread", env.get("LD_LIBRARY_PATH", ""),
])
python = "/workspace/.venvs/gpaw-mpi/bin/python"
out.parent.mkdir(exist_ok=True)
gpaw_info = subprocess.check_output([str(Path(python).with_name("gpaw")), "info"], env=env, text=True)
assert re.search(r"MPI enabled\s+yes", gpaw_info) and re.search(r"scalapack\s+yes", gpaw_info)

launcher = ROOT / "calculations" / f"{label}.launcher.log"
launcher.parent.mkdir(exist_ok=True)
with launcher.open("w") as log:
    child = subprocess.Popen(
        ["mpirun", "--bind-to", "core", "--map-by", "core", "-n", "4", python,
         str(ROOT / "run_pw_reference.py"), label, str(ROOT / record["input"]), str(out)],
        env=env, stdout=log, stderr=subprocess.STDOUT,
    )
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
                subprocess.run(
                    [sys.executable, str(sync.HERE / "sync_tasks.py"), "progress", record["owner_task"], "--job", label,
                     "--state", "running", "--iteration", str(iteration), "--note", "Slab k-mesh/smearing SCF active; four MPI ranks"],
                    check=True,
                )
                last = iteration

if child.returncode:
    subprocess.run(
        [sys.executable, str(sync.HERE / "sync_tasks.py"), "progress", record["owner_task"], "--job", label,
         "--state", "failed", "--iteration", str(last), "--note", f"Calculator exit {child.returncode}; preserve same run/checkpoint"],
        check=True,
    )
    raise SystemExit(child.returncode)

def checksum():
    (archive / "SHA256SUMS.txt").write_text("".join(
        f"{sha(path)}  {path.name}\n" for path in sorted(archive.iterdir())
        if path.is_file() and path.name != "SHA256SUMS.txt" and path.suffix != ".gpw"
    ))

checksum()
verification = json.loads(subprocess.check_output([python, str(ROOT / "verify_result.py"), label, str(archive)], env=env, text=True))
assert verification["status"] == "PASS"
(archive / "verification.json").write_text(json.dumps(verification, indent=2) + "\n")
checksum()
subprocess.run([python, str(ROOT / "verify_result.py"), label, str(archive)], env=env, check=True)
state = out / "state.gpw"
if state.exists():
    state.unlink()
subprocess.run([sys.executable, str(sync.HERE / "sync_tasks.py"), "publish", record["owner_task"], "--label", label], check=True)
print("Assigned slab numerical control verified and published", flush=True)
