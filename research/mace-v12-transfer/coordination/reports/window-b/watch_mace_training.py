#!/usr/bin/env python3
"""Publish coarse-grained MACE epoch progress while a parent training command runs."""
import argparse
import os
from pathlib import Path
import re
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
REPO = Path(subprocess.check_output(["git", "rev-parse", "--show-toplevel"], cwd=HERE, text=True).strip())
SYNC = REPO / "research/mace-v12-transfer/coordination/sync_tasks.py"
EPOCH = re.compile(r"\bEpoch\s+(\d+):")
parser = argparse.ArgumentParser()
parser.add_argument("--parent-pid", type=int, required=True)
parser.add_argument("--interval", type=int, default=45)
args = parser.parse_args()
log = REPO / "research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v14_interface_energy/train.stdout"
last_published = -1
while True:
    try:
        os.kill(args.parent_pid, 0)
    except OSError:
        break
    try:
        contents = log.read_text(errors="replace")
    except OSError:
        contents = ""
    epochs = [int(x) for x in EPOCH.findall(contents)]
    current = max(epochs, default=-1)
    if current >= 0 and current >= last_published + 4:
        completed = min(current + 1, 80)
        subprocess.run([
            sys.executable, str(SYNC), "progress", "window-b",
            "--job", "v14_training", "--state", "running", "--iteration", str(completed),
            "--note", f"MACE CPU training active; epoch {completed}/80.",
        ], cwd=REPO, check=False)
        last_published = current
    time.sleep(args.interval)
