#!/usr/bin/env python3
"""Make one automatic resume attempt after an interrupted v23 DFT sequence."""
from __future__ import annotations

import json
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[3]
PYTHON = "/workspace/.venvs/gpaw-mpi/bin/python"
STATUS_PATH = ROOT / "sequence_status.json"
RESUME = ROOT / "resume_sequence.py"


def active(script_fragment: str) -> bool:
    for proc in Path("/proc").iterdir():
        if not proc.name.isdigit():
            continue
        try:
            cmd = (proc / "cmdline").read_bytes().decode(errors="ignore")
            stat = (proc / "stat").read_text()
            state = stat.rsplit(")", 1)[1].strip().split()[0]
        except (OSError, IndexError):
            continue
        if script_fragment in cmd and state != "Z":
            return True
    return False


def main() -> None:
    stale_observations = 0
    while True:
        try:
            status = json.loads(STATUS_PATH.read_text())
        except (OSError, json.JSONDecodeError):
            time.sleep(20)
            continue
        if status.get("status") == "all_six_archives_verified":
            print("all v23 labels verified; recovery monitor finished", flush=True)
            return
        if status.get("status") == "failed_stop_no_retry":
            label = status.get("failed_label") or status.get("current_label")
            count = status.get("automatic_resume_attempts_by_label", {}).get(label, 0)
            if count >= 1:
                print(f"one automatic resume already used for {label}; preserving checkpoint for manual review",
                      flush=True)
                return
            if active("run_sequence.py") or active("run_preflight_label.py"):
                time.sleep(20)
                continue
            print(f"attempting one resume for interrupted label {label}", flush=True)
            proc = subprocess.run([PYTHON, str(RESUME), "--automatic"], cwd=REPO)
            print(f"resume process returned {proc.returncode}", flush=True)
            stale_observations = 0
            continue
        if status.get("status") == "running" and status.get("current_label"):
            if not active("run_sequence.py") and not active("run_preflight_label.py"):
                stale_observations += 1
                if stale_observations >= 3:
                    print("sequence supervisor disappeared with an active label; attempting one resume",
                          flush=True)
                    proc = subprocess.run([PYTHON, str(RESUME), "--automatic"], cwd=REPO)
                    print(f"resume process returned {proc.returncode}", flush=True)
                    stale_observations = 0
            else:
                stale_observations = 0
        else:
            stale_observations = 0
        time.sleep(20)


if __name__ == "__main__":
    main()
