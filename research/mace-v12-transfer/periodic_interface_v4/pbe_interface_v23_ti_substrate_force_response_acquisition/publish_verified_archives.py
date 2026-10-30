#!/usr/bin/env python3
"""Commit and push each v23 result only after its archive verifier passes."""
from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[3]
PYTHON = "/workspace/.venvs/gpaw-mpi/bin/python"
VERIFY = ROOT / "verify_preflight.py"
MANIFEST = json.loads((ROOT / "input_manifest.json").read_text())
STATUS = ROOT / "sequence_status.json"


def git(*args: str, check: bool = True) -> subprocess.CompletedProcess[str]:
    proc = subprocess.run(["git", "-C", str(REPO), *args], text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if check and proc.returncode:
        raise RuntimeError(f"git {' '.join(args)} failed ({proc.returncode}): {proc.stdout}")
    return proc


def verify(label: str) -> bool:
    proc = subprocess.run([PYTHON, str(VERIFY), label], cwd=REPO,
                          text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if proc.returncode:
        print(f"{label}: archive not ready for publication: {proc.stdout.strip()}", flush=True)
        return False
    result = json.loads(proc.stdout)
    if result.get("status") != "PASS":
        print(f"{label}: verifier did not return PASS", flush=True)
        return False
    return True


def publish(label: str) -> None:
    folder = ROOT / "calculations" / label
    if not (folder / "verification.json").is_file() or not (folder / "SHA256SUMS.txt").is_file():
        raise RuntimeError(f"{label}: missing verification or checksum file")
    rel = folder.relative_to(REPO).as_posix()
    if git("diff", "--cached", "--quiet", check=False).returncode != 0:
        raise RuntimeError("Refusing to publish while another change is already staged")
    git("add", "--sparse", "--", rel)
    staged = git("diff", "--cached", "--name-only").stdout.splitlines()
    if not staged or any(not (p == rel or p.startswith(rel + "/")) for p in staged):
        raise RuntimeError(f"Unexpected staged paths while publishing {label}: {staged}")
    git("commit", "-m", f"Add verified v23 DFT result for {label}")
    while True:
        proc = git("push", "origin", "HEAD:main", check=False)
        print(proc.stdout.strip(), flush=True)
        if proc.returncode == 0:
            break
        print(f"{label}: push failed; keeping the verified local commit and retrying in 60 s",
              flush=True)
        time.sleep(60)


def main() -> None:
    jobs = MANIFEST["preflight_jobs"]
    published = set()
    pending_push: str | None = None
    while True:
        status = json.loads(STATUS.read_text())
        completed = {row["label"] for row in status.get("completed_labels", [])
                     if row.get("archive_verification") == "PASS"}
        if pending_push is not None:
            publish(pending_push)
            print(f"{pending_push}: verified result pushed to main", flush=True)
            published.add(pending_push)
            pending_push = None
        for job in jobs:
            label = job["label"]
            if label in completed and label not in published:
                if verify(label):
                    pending_push = label
                    break
        if status.get("status") in {"all_six_archives_verified", "failed_stop_no_retry"}:
            if pending_push is None and completed.issubset(published):
                print(f"publisher finished with status={status['status']}; "
                      f"published {len(published)} verified result archives", flush=True)
                return
        time.sleep(20)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"publisher stopped: {exc}", file=sys.stderr, flush=True)
        raise
