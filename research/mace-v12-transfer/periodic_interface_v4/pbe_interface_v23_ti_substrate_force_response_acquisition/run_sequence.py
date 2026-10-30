#!/usr/bin/env python3
"""Recover or run the six v23 labels, with durable logs and no automatic retries."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[3]
PYTHON = "/workspace/.venvs/gpaw-mpi/bin/python"
CGROUP = Path("/sys/fs/cgroup")
STOP = 30 * 1024**3
RUNNER = "run_preflight_label.py"
LAUNCH = ROOT / "launch_preflight.py"
VERIFY = ROOT / "verify_preflight.py"
MANIFEST = json.loads((ROOT / "input_manifest.json").read_text())
JOBS = MANIFEST["preflight_jobs"]
PRIORITY = {"AgTi": 0, "AgSi": 1, "AgC": 2}
JOBS = sorted(JOBS, key=lambda x: (PRIORITY[x["label"].split("_", 1)[0]],
                                   0 if "_m0p02_" in x["label"] else 1))
STATUS_PATH = ROOT / "sequence_status.json"
CALCULATIONS = ROOT / "calculations"


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_status() -> dict:
    if STATUS_PATH.exists():
        return json.loads(STATUS_PATH.read_text())
    return {"status": "running", "started_utc": utc(), "completed_labels": [],
            "current_label": None, "automatic_retries": 0}


status = read_status()
status["status"] = "running"
status["supervisor_started_utc"] = utc()
status["automatic_retries"] = 0


def save_status() -> None:
    tmp = STATUS_PATH.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(status, indent=2) + "\n")
    tmp.replace(STATUS_PATH)


def memory_current() -> int:
    return int((CGROUP / "memory.current").read_text())


def memory_peak() -> int:
    p = CGROUP / "memory.peak"
    return int(p.read_text()) if p.exists() else memory_current()


def runner_pids(label: str) -> list[int]:
    found = []
    for proc in Path("/proc").iterdir():
        if not proc.name.isdigit():
            continue
        try:
            cmd = (proc / "cmdline").read_bytes().decode(errors="ignore")
            stat = (proc / "stat").read_text()
            state = stat.rsplit(")", 1)[1].strip().split()[0]
        except (OSError, IndexError):
            continue
        if RUNNER in cmd and label in cmd and state != "Z":
            found.append(int(proc.name))
    return found


def stop_runners(label: str) -> None:
    groups = set()
    for pid in runner_pids(label):
        try:
            groups.add(os.getpgid(pid))
        except ProcessLookupError:
            pass
    for pgid in groups:
        try:
            os.killpg(pgid, signal.SIGTERM)
        except ProcessLookupError:
            pass
    deadline = time.monotonic() + 30
    while time.monotonic() < deadline and runner_pids(label):
        time.sleep(1)
    for pgid in groups:
        try:
            os.killpg(pgid, signal.SIGKILL)
        except ProcessLookupError:
            pass


def write_sums(folder: Path) -> None:
    rows = []
    for path in sorted(folder.iterdir()):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            rows.append(f"{digest(path)}  {path.name}\n")
    (folder / "SHA256SUMS.txt").write_text("".join(rows))


def verify_archive(label: str) -> dict:
    folder = CALCULATIONS / label
    proc = subprocess.run([PYTHON, str(VERIFY), label], cwd=REPO,
                          capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"archive verifier failed for {label}: {proc.stderr or proc.stdout}")
    first = json.loads(proc.stdout)
    if first.get("status") != "PASS":
        raise RuntimeError(f"archive verifier did not pass for {label}")
    verification = folder / "verification.json"
    if not verification.exists():
        verification.write_text(json.dumps({
            "status": "PASS", "label": label, "checks": first["checks"],
            "verifier_sha256": digest(VERIFY), "verification_record_note":
                "Archive was verified, then this record was added and included in SHA256SUMS.txt; the archive was verified again."
        }, indent=2) + "\n")
        write_sums(folder)
    second = subprocess.run([PYTHON, str(VERIFY), label], cwd=REPO,
                            capture_output=True, text=True)
    if second.returncode != 0:
        raise RuntimeError(f"final archive verifier failed for {label}: {second.stderr or second.stdout}")
    result = json.loads(second.stdout)
    if result.get("status") != "PASS":
        raise RuntimeError(f"final archive verification status is not PASS for {label}")
    return result


def record_pass(label: str, result: dict) -> None:
    rows = status.setdefault("completed_labels", [])
    if not any(row.get("label") == label for row in rows):
        rows.append({"label": label, "archive_verification": "PASS",
                     "archive_sha256": result["archive_sha256"]})
    status["current_label"] = None
    save_status()


def preserve_failed_attempt(label: str, reason: str, started: float,
                            peak: int, stopped_for_memory: bool) -> None:
    scratch = CALCULATIONS / f".running_{label}"
    if scratch.exists():
        (scratch / "state.next.gpw").unlink(missing_ok=True)
        failed_root = CALCULATIONS / "failed_attempts"
        failed_root.mkdir(exist_ok=True)
        failed = failed_root / f"{label}__attempt_{int(started)}"
        if failed.exists():
            failed = failed_root / f"{label}__attempt_{int(started)}_{int(time.time())}"
        if scratch.exists():
            scratch.rename(failed)
        launcher_log = CALCULATIONS / f"{label}.launcher.log"
        if launcher_log.exists():
            (failed / "launcher.log").write_bytes(launcher_log.read_bytes())
        (failed / "runtime.json").write_text(json.dumps({
            "label": label, "return_code": None,
            "stopped_for_memory": stopped_for_memory,
            "elapsed_s": time.time() - started,
            "memory_high_water_observed_bytes": peak,
            "cgroup_limit_bytes": int((CGROUP / "memory.max").read_text()),
            "memory_stop_threshold_bytes": STOP,
            "recovery_supervisor_reason": reason,
        }, indent=2) + "\n")


def archive_recovered(label: str, peak: int) -> dict:
    scratch = CALCULATIONS / f".running_{label}"
    target = CALCULATIONS / label
    if target.exists():
        raise RuntimeError(f"refusing to overwrite existing archive: {target}")
    progress = json.loads((scratch / "progress.json").read_text())
    if progress.get("status") != "complete" or not (scratch / "summary.json").is_file():
        raise RuntimeError(f"running attempt lacks completed summary for {label}")
    (scratch / "state.gpw").unlink(missing_ok=True)
    started = float(progress.get("started_epoch_s", time.time()))
    runtime = {
        "label": label, "return_code": 0, "stopped_for_memory": False,
        "elapsed_s": time.time() - started,
        "memory_start_bytes": None,
        "memory_high_water_observed_bytes": peak,
        "cgroup_limit_bytes": int((CGROUP / "memory.max").read_text()),
        "memory_stop_threshold_bytes": STOP,
        "recovered_after_launcher_exit": True,
    }
    (scratch / "runtime.json").write_text(json.dumps(runtime, indent=2) + "\n")
    scratch.rename(target)
    write_sums(target)
    return verify_archive(label)


def recover_current(label: str) -> dict:
    scratch = CALCULATIONS / f".running_{label}"
    if not scratch.exists():
        target = CALCULATIONS / label
        if target.exists():
            return verify_archive(label)
        raise RuntimeError(f"no running scratch or archive exists for current label {label}")
    progress_path = scratch / "progress.json"
    started = float(json.loads(progress_path.read_text()).get("started_epoch_s", time.time()))
    peak = memory_peak()
    initial = memory_current()
    print(f"Recovering active DFT {label}; monitor memory guard at 30 GiB.", flush=True)
    while True:
        progress = json.loads(progress_path.read_text())
        current = memory_current()
        peak = max(peak, current, memory_peak())
        pids = runner_pids(label)
        if current >= STOP:
            stop_runners(label)
            preserve_failed_attempt(label, f"memory {current} crossed 30 GiB stop guard",
                                    started, peak, True)
            raise RuntimeError(f"memory guard stopped {label}")
        if progress.get("status") == "complete":
            deadline = time.monotonic() + 60
            while runner_pids(label) and time.monotonic() < deadline:
                time.sleep(2)
            if runner_pids(label):
                raise RuntimeError(f"SCF completed but MPI ranks did not exit for {label}")
            return archive_recovered(label, peak)
        if not pids:
            preserve_failed_attempt(label, "MPI ranks exited before a complete progress record",
                                    started, peak, False)
            raise RuntimeError(f"active DFT ended without a verified result for {label}")
        if int(time.time()) % 60 < 10:
            print(f"{label}: SCF={progress.get('iteration', 0)} (updates every 10 steps); "
                  f"memory={current/1024**3:.2f} GiB; peak={peak/1024**3:.2f} GiB", flush=True)
        time.sleep(10)


def run_label(label: str) -> dict:
    target = CALCULATIONS / label
    scratch = CALCULATIONS / f".running_{label}"
    if target.exists():
        return verify_archive(label)
    if scratch.exists():
        raise RuntimeError(f"unregistered partial scratch exists for {label}; preserving and stopping")
    status["current_label"] = label
    save_status()
    log_path = CALCULATIONS / f"{label}.supervisor.log"
    CALCULATIONS.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as log:
        proc = subprocess.Popen([PYTHON, str(LAUNCH), label], cwd=REPO,
                                stdout=log, stderr=subprocess.STDOUT,
                                start_new_session=True)
        peak = memory_peak()
        start = time.time()
        while proc.poll() is None:
            current = memory_current()
            peak = max(peak, current, memory_peak())
            if current >= STOP:
                stop_runners(label)
            time.sleep(10)
        code = proc.wait()
    if code != 0:
        status.update({"status": "failed_stop_no_retry", "failed_label": label,
                       "launcher_exit_code": code, "finished_utc": utc()})
        save_status()
        raise RuntimeError(f"launcher failed for {label}, exit={code}; no retry")
    result = verify_archive(label)
    record_pass(label, result)
    return result


def main() -> None:
    CALCULATIONS.mkdir(parents=True, exist_ok=True)
    current = status.get("current_label")
    completed = {row["label"] for row in status.get("completed_labels", [])}
    if current and current not in completed:
        print(f"Resuming registered current item: {current}", flush=True)
        result = recover_current(current)
        record_pass(current, result)
        completed.add(current)
        print(f"{current}: archive verification PASS", flush=True)
    for job in JOBS:
        label = job["label"]
        if label in completed:
            continue
        print(f"=== {label}: resource preflight and DFT ===", flush=True)
        result = run_label(label)
        completed.add(label)
        print(f"{label}: archive verification PASS", flush=True)
    status.update({"status": "all_six_archives_verified", "current_label": None,
                   "finished_utc": utc(), "automatic_retries": 0})
    save_status()
    print(json.dumps(status, indent=2), flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        status.setdefault("status", "failed_stop_no_retry")
        status["status"] = "failed_stop_no_retry"
        status["supervisor_error"] = str(exc)
        status["finished_utc"] = utc()
        save_status()
        print(f"STOPPED WITHOUT RETRY: {exc}", file=sys.stderr, flush=True)
        raise
