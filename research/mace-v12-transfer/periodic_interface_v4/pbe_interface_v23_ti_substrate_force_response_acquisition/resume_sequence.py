#!/usr/bin/env python3
"""Resume the first incomplete v23 label, using a matching latest checkpoint when available."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[3]
PYTHON = "/workspace/.venvs/gpaw-mpi/bin/python"
STATUS_PATH = ROOT / "sequence_status.json"
MANIFEST = json.loads((ROOT / "input_manifest.json").read_text())
JOBS = {job["label"]: job for job in MANIFEST["preflight_jobs"]}
CHECKPOINT_DIR = ROOT / "calculations" / "latest_checkpoint"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def running_processes() -> list[str]:
    found = []
    for proc in Path("/proc").iterdir():
        if not proc.name.isdigit() or int(proc.name) == os.getpid():
            continue
        try:
            cmd = (proc / "cmdline").read_bytes().decode(errors="ignore")
        except OSError:
            continue
        if "run_sequence.py" in cmd or "run_preflight_label.py" in cmd:
            found.append(proc.name)
    return found


def snapshot_failed_scratch(label: str) -> None:
    candidates = [ROOT / "calculations" / f".running_{label}"]
    candidates.extend(sorted((ROOT / "calculations" / "failed_attempts").glob(f"{label}__*/"),
                             key=lambda p: p.stat().st_mtime, reverse=True))
    for folder in candidates:
        source = folder / "state.gpw"
        if not source.is_file() or source.stat().st_size <= 0:
            continue
        job = JOBS[label]
        CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)
        target = CHECKPOINT_DIR / "latest_state.gpw"
        tmp = CHECKPOINT_DIR / "latest_state.gpw.tmp"
        shutil.copyfile(source, tmp)
        if sha(source) != sha(tmp):
            tmp.unlink(missing_ok=True)
            raise RuntimeError("failed-attempt checkpoint changed while copying; refusing to resume")
        checkpoint_hash = sha(tmp)
        os.replace(tmp, target)
        try:
            progress = json.loads((folder / "progress.json").read_text())
        except (OSError, json.JSONDecodeError):
            progress = {}
        metadata = {
            "status": "stable_checkpoint_snapshot", "label": label,
            "input_sha256": job["input_sha256"], "checkpoint_sha256": checkpoint_hash,
            "checkpoint_bytes": target.stat().st_size,
            "last_reported_scf_iteration": progress.get("iteration"),
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "resume_advice": "Recovered from an inactive failed-attempt directory.",
        }
        temp_meta = CHECKPOINT_DIR / "latest_checkpoint.json.tmp"
        temp_meta.write_text(json.dumps(metadata, indent=2) + "\n")
        os.replace(temp_meta, CHECKPOINT_DIR / "latest_checkpoint.json")
        print(f"Saved checkpoint from {folder} at iteration≈{progress.get('iteration')}", flush=True)
        return


def main() -> int:
    automatic = "--automatic" in sys.argv[1:]
    if running_processes():
        raise RuntimeError(f"refusing to resume while DFT/supervisor processes are active: {running_processes()}")
    status = json.loads(STATUS_PATH.read_text())
    label = status.get("failed_label") or status.get("current_label")
    if label not in JOBS:
        raise RuntimeError(f"no registered incomplete label to resume: {label}")
    automatic_counts = status.setdefault("automatic_resume_attempts_by_label", {})
    if automatic and automatic_counts.get(label, 0) >= 1:
        print(f"automatic resume limit already reached for {label}; latest checkpoint is preserved", flush=True)
        return 2

    scratch = ROOT / "calculations" / f".running_{label}"
    if scratch.exists():
        failed_root = ROOT / "calculations" / "failed_attempts"
        failed_root.mkdir(exist_ok=True)
        destination = failed_root / f"{label}__interrupted_{int(time.time())}"
        scratch.rename(destination)
    snapshot_failed_scratch(label)

    checkpoint = CHECKPOINT_DIR / "latest_state.gpw"
    meta_path = CHECKPOINT_DIR / "latest_checkpoint.json"
    if checkpoint.is_file() and meta_path.is_file():
        meta = json.loads(meta_path.read_text())
        valid = (meta.get("label") == label
                 and meta.get("input_sha256") == JOBS[label]["input_sha256"]
                 and meta.get("checkpoint_sha256") == sha(checkpoint))
    else:
        valid = False
    status.setdefault("resume_history", []).append({
        "label": label, "requested_utc": datetime.now(timezone.utc).isoformat(),
        "checkpoint_used": valid, "checkpoint_sha256": meta.get("checkpoint_sha256") if valid else None,
        "automatic": automatic,
    })
    if automatic:
        automatic_counts[label] = automatic_counts.get(label, 0) + 1
    status["status"] = "running"
    status["current_label"] = None
    for key in ("failed_label", "launcher_exit_code", "finished_utc"):
        status.pop(key, None)
    status["resume_from_label"] = label
    tmp = STATUS_PATH.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(status, indent=2) + "\n")
    tmp.replace(STATUS_PATH)
    print(f"Resuming sequence at {label}; matching checkpoint={'yes' if valid else 'no, fresh start'}",
          flush=True)
    proc = subprocess.run([PYTHON, str(ROOT / "run_sequence.py")], cwd=REPO)
    return proc.returncode


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"resume refused: {exc}", file=sys.stderr, flush=True)
        raise
