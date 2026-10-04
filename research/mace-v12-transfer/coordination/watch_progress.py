#!/usr/bin/env python3
"""Publish small per-window SCF progress snapshots while a label queue runs."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time


HERE = Path(__file__).resolve().parent
REPO = Path(subprocess.check_output(["git", "rev-parse", "--show-toplevel"], cwd=HERE, text=True).strip())
SYNC = HERE / "sync_tasks.py"
ITER_RE = re.compile(r"\|iter:\s*(\d+)\|")


def read_json(path):
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None


def iteration_from_log(path):
    try:
        lines = path.read_text(errors="replace").splitlines()
    except OSError:
        return 0
    for line in reversed(lines[-100:]):
        match = ITER_RE.search(line)
        if match:
            return int(match.group(1))
    return 0


def parent_is_alive(pid):
    if pid is None:
        return True
    try:
        os.kill(pid, 0)
        return True
    except OSError:
        return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task", choices=("window-a", "window-b"))
    parser.add_argument("run_root", type=Path)
    parser.add_argument("archive_root", type=Path)
    parser.add_argument("labels", nargs="+")
    parser.add_argument("--step", type=int, default=10, help="Publish after each N additional SCF iterations.")
    parser.add_argument("--interval", type=int, default=45, help="Seconds between local progress checks.")
    parser.add_argument("--parent-pid", type=int)
    args = parser.parse_args()
    if args.step < 1 or args.interval < 10:
        raise SystemExit("Use a positive iteration step and an interval of at least 10 seconds.")
    run_root = args.run_root.resolve()
    archive_root = args.archive_root.resolve()
    print(f"Watching {args.task} queue; publishing every {args.step} SCF iterations.", flush=True)
    last_published = {}
    terminal_published = set()
    queued_published = set()

    while parent_is_alive(args.parent_pid):
        active = []
        terminal = 0
        for label in args.labels:
            run_dir = run_root / label
            archive_dir = archive_root / label
            run_progress = read_json(run_dir / "progress.json")
            archived_progress = read_json(archive_dir / "progress.json")
            status_record = run_progress or archived_progress
            if status_record is None:
                if label not in queued_published:
                    command = [
                        sys.executable, str(SYNC), "progress", args.task,
                        "--job", label, "--state", "queued", "--iteration", "0",
                        "--note", "Assigned label is queued; waiting for its SCF run.",
                    ]
                    result = subprocess.run(command, cwd=REPO, text=True, capture_output=True)
                    if result.returncode == 0:
                        queued_published.add(label)
                        print(f"Published {label} queued.", flush=True)
                    else:
                        print(f"Queued progress sync failed for {label}: {result.stderr.strip() or result.stdout.strip()}", flush=True)
                continue
            state = status_record.get("status", "unknown")
            if state == "running":
                log_iteration = iteration_from_log(run_dir / "gpaw.log")
                iteration = max(int(status_record.get("iteration", 0)), log_iteration)
                elapsed = status_record.get("elapsed_s")
                if elapsed is None and status_record.get("started_utc_epoch_s") is not None:
                    elapsed = max(0.0, time.time() - float(status_record["started_utc_epoch_s"]))
                active.append((label, iteration, elapsed, status_record))
            elif state in ("complete", "not_converged", "failed"):
                terminal += 1
                progress_state = "completed" if state == "complete" else state
                if label not in terminal_published:
                    iteration = int(status_record.get("iteration", 0))
                    elapsed = status_record.get("elapsed_s")
                    note = (
                        "SCF converged; compact archive verification/publication is pending."
                        if state == "complete"
                        else f"SCF reached terminal state {state}; inspect the local log before resuming."
                    )
                    command = [
                        sys.executable, str(SYNC), "progress", args.task,
                        "--job", label, "--state", progress_state,
                        "--iteration", str(iteration), "--note", note,
                    ]
                    if elapsed is not None:
                        command.extend(("--elapsed-s", str(float(elapsed))))
                    result = subprocess.run(command, cwd=REPO, text=True, capture_output=True)
                    if result.returncode == 0:
                        terminal_published.add(label)
                        print(f"Published {label} terminal SCF state {progress_state}.", flush=True)
                    else:
                        print(f"Terminal progress sync failed for {label}: {result.stderr.strip() or result.stdout.strip()}", flush=True)

        if not active and terminal == len(args.labels):
            print("All assigned labels reached a terminal SCF state; stopping progress watcher.", flush=True)
            return

        for label, iteration, elapsed, record in active:
            previous = last_published.get(label, -args.step)
            if iteration < previous + args.step and label in last_published:
                continue
            ranks = record.get("mpi_ranks", "unknown")
            note = f"SCF active; {ranks} MPI ranks"
            command = [
                sys.executable, str(SYNC), "progress", args.task,
                "--job", label, "--state", "running", "--iteration", str(iteration),
                "--note", note,
            ]
            if elapsed is not None:
                command.extend(("--elapsed-s", str(float(elapsed))))
            result = subprocess.run(command, cwd=REPO, text=True, capture_output=True)
            if result.returncode == 0:
                last_published[label] = iteration
                print(f"Published {label} iteration {iteration}.", flush=True)
            else:
                print(f"Progress sync failed for {label} iteration {iteration}: {result.stderr.strip() or result.stdout.strip()}", flush=True)

        time.sleep(args.interval)

    print("Runner process ended; stopping progress watcher.", flush=True)


if __name__ == "__main__":
    main()
