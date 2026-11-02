#!/usr/bin/env python3
"""Keep one stable, checksummed copy of the newest active GPAW checkpoint."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import time
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parent
CALCULATIONS = ROOT / "calculations"
DEST = CALCULATIONS / "latest_checkpoint"
STATUS_PATH = ROOT / "sequence_status.json"
MANIFEST = json.loads((ROOT / "input_manifest.json").read_text())
JOBS = {job["label"]: job for job in MANIFEST["preflight_jobs"]}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_json(path: Path, value: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2) + "\n")
    tmp.replace(path)


def main() -> None:
    DEST.mkdir(parents=True, exist_ok=True)
    stable_key = None
    stable_count = 0
    saved_key = None
    while True:
        try:
            status = json.loads(STATUS_PATH.read_text())
        except (OSError, json.JSONDecodeError):
            time.sleep(10)
            continue
        label = status.get("current_label")
        if status.get("status") == "failed_stop_no_retry":
            print(f"sequence stopped; leaving the newest verified checkpoint at {DEST}", flush=True)
            return
        if label not in JOBS:
            if status.get("status") in {"all_six_archives_verified", "failed_stop_no_retry"}:
                print(f"checkpoint keeper stopped; last snapshot retained at {DEST}", flush=True)
                return
            stable_key, stable_count = None, 0
            time.sleep(10)
            continue

        source = CALCULATIONS / f".running_{label}" / "state.gpw"
        try:
            st = source.stat()
            key = (label, st.st_size, st.st_mtime_ns)
        except OSError:
            stable_key, stable_count = None, 0
            time.sleep(10)
            continue
        if st.st_size <= 0:
            stable_key, stable_count = None, 0
            time.sleep(10)
            continue
        if key == stable_key:
            stable_count += 1
        else:
            stable_key, stable_count = key, 0
        if stable_count >= 2 and key != saved_key:
            tmp = DEST / "latest_state.gpw.tmp"
            target = DEST / "latest_state.gpw"
            try:
                before = source.stat()
                shutil.copyfile(source, tmp)
                after = source.stat()
                if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                    tmp.unlink(missing_ok=True)
                    stable_key, stable_count = None, 0
                    time.sleep(10)
                    continue
                source_hash = digest(source)
                copy_hash = digest(tmp)
                if source_hash != copy_hash:
                    tmp.unlink(missing_ok=True)
                    stable_key, stable_count = None, 0
                    time.sleep(10)
                    continue
                os.replace(tmp, target)
                progress_path = CALCULATIONS / f".running_{label}" / "progress.json"
                try:
                    progress = json.loads(progress_path.read_text())
                except (OSError, json.JSONDecodeError):
                    progress = {}
                meta = {
                    "status": "stable_checkpoint_snapshot",
                    "label": label,
                    "input_sha256": JOBS[label]["input_sha256"],
                    "checkpoint_sha256": copy_hash,
                    "checkpoint_bytes": target.stat().st_size,
                    "last_reported_scf_iteration": progress.get("iteration"),
                    "created_utc": datetime.now(timezone.utc).isoformat(),
                    "source_file_mtime_ns": after.st_mtime_ns,
                    "resume_advice": "Use resume_sequence.py; it checks label and input hashes before reuse.",
                }
                atomic_json(DEST / "latest_checkpoint.json", meta)
                saved_key = key
                print(f"saved latest checkpoint for {label}: "
                      f"SCF≈{meta['last_reported_scf_iteration']}, {target.stat().st_size} bytes",
                      flush=True)
            except OSError as exc:
                tmp.unlink(missing_ok=True)
                print(f"checkpoint snapshot deferred: {exc}", flush=True)
        time.sleep(10)


if __name__ == "__main__":
    main()
