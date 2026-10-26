#!/usr/bin/env python3
"""Run the separately registered B1 k-mesh task sequentially, pair by pair."""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[5]
TASK = json.loads((HERE / "task_manifest.json").read_text())
PYTHON = "/workspace/.venvs/gpaw-mpi/bin/python"
MPI_LOCK = Path("/workspace/.mace-v19-dft.lock")
OUT = HERE / TASK["outputs_directory"]
OUT.mkdir(parents=True, exist_ok=True)
(HERE / "comparisons").mkdir(exist_ok=True)
(HERE / "launcher_logs").mkdir(exist_ok=True)

def write_progress(**updates):
    path = HERE / "task_progress.json"
    current = json.loads(path.read_text()) if path.exists() else {
        "task_id": TASK["task_id"], "status": "queued", "started_epoch_s": time.time(),
        "labels": {}, "pairs": {}, "updated_epoch_s": time.time(),
    }
    current.update(updates)
    current["updated_epoch_s"] = time.time()
    path.write_text(json.dumps(current, indent=2) + "\n")
    return current

def env_for_gpaw():
    env = os.environ.copy()
    deps = "/workspace/.local/gpaw-mpi-deps/usr/lib/x86_64-linux-gnu"
    env.update({
        "GPAW_MPI_BACKEND":"cgpaw", "OMP_NUM_THREADS":"1", "OPENBLAS_NUM_THREADS":"1",
        "MKL_NUM_THREADS":"1", "NUMEXPR_NUM_THREADS":"1", "XDG_CACHE_HOME":"/workspace/.cache",
        "MPLCONFIGDIR":"/workspace/.cache/matplotlib", "FC_CACHEDIR":"/workspace/.cache/fontconfig",
    })
    env["LD_LIBRARY_PATH"] = ":".join([
        deps + "/openmpi/lib", deps, "/workspace/.local/gpaw-libs/usr/lib/x86_64-linux-gnu",
        "/workspace/.local/openblas/usr/lib/x86_64-linux-gnu/openblas-pthread", env.get("LD_LIBRARY_PATH", ""),
    ])
    return env

def ensure_run(record, conditional=False):
    label = record["label"]
    run_dir = OUT / label
    if run_dir.exists():
        # Resume only from an already complete, independently verified archive.
        subprocess.run([PYTHON, str(HERE / "verify_job.py"), label], check=True, capture_output=True, text=True)
        write_progress(labels={**read_progress().get("labels", {}), label:"complete_verified"})
        return
    cmd = ["mpirun", "--bind-to", "core", "--map-by", "core", "-n", "4", PYTHON, str(HERE / "run_one.py"), label]
    if conditional:
        cmd.append("--conditional")
    logfile = HERE / "launcher_logs" / f"{label}.launcher.log"
    write_progress(status="running", active_label=label,
                   labels={**read_progress().get("labels", {}), label:"running"})
    with logfile.open("w") as stream:
        result = subprocess.run(cmd, cwd=ROOT, env=env_for_gpaw(), stdout=stream, stderr=subprocess.STDOUT)
    if result.returncode:
        write_progress(status="failed", active_label=label,
                       labels={**read_progress().get("labels", {}), label:f"failed_exit_{result.returncode}"})
        raise SystemExit(f"{label} failed; preserve output and inspect {logfile}")
    verification = json.loads(subprocess.check_output([PYTHON, str(HERE / "verify_job.py"), label], cwd=ROOT, text=True))
    assert verification["status"] == "PASS"
    run_dir.joinpath("verification.json").write_text(json.dumps(verification, indent=2) + "\n")
    files = sorted(p for p in run_dir.iterdir() if p.is_file() and p.name not in ("SHA256SUMS.txt", "state.gpw"))
    run_dir.joinpath("SHA256SUMS.txt").write_text("".join(
        f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n" for p in files
    ))
    subprocess.run([PYTHON, str(HERE / "verify_job.py"), label], cwd=ROOT, check=True, capture_output=True, text=True)
    state = run_dir / "state.gpw"
    if state.exists():
        state.unlink()
    write_progress(labels={**read_progress().get("labels", {}), label:"complete_verified"}, active_label=None)

def read_progress():
    path = HERE / "task_progress.json"
    return json.loads(path.read_text()) if path.exists() else {}

def compare(coarse, fine, outfile):
    result = json.loads(subprocess.check_output([
        PYTHON, str(HERE / "compare_pair.py"), coarse, fine, "--output", outfile
    ], cwd=ROOT, text=True))
    write_progress(pairs={**read_progress().get("pairs", {}), outfile:result})
    return result

def main():
    if len(sys.argv) != 1:
        raise SystemExit("This registered task runs the full queue; no ad hoc labels are accepted.")
    subprocess.run([PYTHON, str(HERE / "verify_task.py")], cwd=ROOT, check=True)
    subprocess.run([PYTHON, str(BASE_PACK_VERIFIER)], cwd=ROOT, check=True)
    with MPI_LOCK.open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise SystemExit("Another registered GPAW task holds the shared four-rank DFT lock")
        write_progress(status="running", active_label=None)
        pairs = {}
        for record in TASK["records"]:
            pairs.setdefault(record["pair_id"], []).append(record)
        failures = []
        for pair_id, members in pairs.items():
            members = sorted(members, key=lambda r: 0 if r["pair_member"] == "coarse" else 1)
            for record in members:
                ensure_run(record)
            coarse, fine = members
            comparison = compare(coarse["label"], fine["label"], f"{pair_id}.json")
            if comparison["pass"] is False:
                failures.append(pair_id)
                if coarse["phase"] == "liquid_Ag":
                    conditional = next(r for r in TASK["conditional_records"] if r["compare_to_label"] == fine["label"])
                    trigger_path = HERE / "comparisons" / f"{conditional['pair_id']}.json"
                    trigger_path.write_text(json.dumps(comparison, indent=2) + "\n")
                    ensure_run(conditional, conditional=True)
                    conditional_comparison = compare(fine["label"], conditional["label"], f"{conditional['pair_id']}_2x2x2_vs_3x3x3.json")
                    if conditional_comparison["pass"] is False:
                        failures.append(conditional["pair_id"])
        status = "completed_with_mesh_threshold_failures" if failures else "completed_all_mesh_pairs_pass"
        write_progress(status=status, active_label=None, threshold_failures=failures, finished_epoch_s=time.time())
        print(json.dumps(read_progress(), indent=2))

BASE_PACK_VERIFIER = ROOT / "research/mace-v12-transfer/coordination/reports/window-c/b1_kmesh_execution_pack_20261024/verify_pack.py"

if __name__ == "__main__":
    main()
