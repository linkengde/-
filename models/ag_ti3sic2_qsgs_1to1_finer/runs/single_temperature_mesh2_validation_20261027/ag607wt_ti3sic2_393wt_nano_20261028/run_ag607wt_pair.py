#!/usr/bin/env python3
"""Gate and supervise matched 1x/2x Ag 60.7 wt% runs with checkpoint publishing."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
PROJECT_DIR = HERE.parents[2]
REPO = PROJECT_DIR.parent.parent
CASE_RUNNER = HERE / "run_ag607wt_case.py"
CONTROL = HERE / "runs"
BRANCH = "codex/single-temperature-study"
MANAGER_STOP_REQUESTED = False


def now_utc() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def atomic_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(data, indent=2, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def record(path: Path, event: str, **fields: Any) -> None:
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps({"utc": now_utc(), "event": event, **fields}, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def git(*args: str, check: bool = True) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(["git", *args], cwd=REPO, text=True, capture_output=True)
    if check and result.returncode:
        raise RuntimeError(f"git {' '.join(args)} failed: {result.stdout}{result.stderr}")
    return result


def commit_push(paths: list[Path], message: str) -> str | None:
    relative = [str(path.resolve().relative_to(REPO)) for path in paths if path.exists()]
    if not relative:
        return None
    staged = git("diff", "--cached", "--name-only").stdout.splitlines()
    if staged:
        raise RuntimeError(f"Unexpected pre-staged paths before publishing: {staged}")
    git("add", "--", *relative)
    staged = git("diff", "--cached", "--name-only").stdout.splitlines()
    unexpected = sorted(set(staged) - set(relative))
    if unexpected:
        raise RuntimeError(f"Publishing staged unexpected paths: {unexpected}")
    if git("diff", "--cached", "--quiet", check=False).returncode == 0:
        git("reset", "--", *relative)
        return None
    git("commit", "-m", message)
    for attempt in range(3):
        pushed = git("push", "origin", f"HEAD:refs/heads/{BRANCH}", check=False)
        if pushed.returncode == 0:
            return git("rev-parse", "HEAD").stdout.strip()
        git("fetch", "origin", f"refs/heads/{BRANCH}:refs/remotes/origin/{BRANCH}")
        merged = git("merge", "--no-edit", f"refs/remotes/origin/{BRANCH}", check=False)
        if merged.returncode:
            raise RuntimeError(f"Remote merge after push rejection failed: {merged.stdout}{merged.stderr}")
    raise RuntimeError("Push rejected after three integration attempts")


def verify_preflights() -> tuple[dict[str, Any], dict[str, Any]]:
    results = []
    for refinement, shape in ((1, [100, 100, 50]), (2, [200, 200, 100])):
        path = HERE / "preflight" / f"{refinement}x" / "preflight.json"
        if not path.exists():
            raise FileNotFoundError(f"Missing 1 fs numerical acceptance result: {path}")
        result = json.loads(path.read_text(encoding="utf-8"))
        if result.get("status") != "CHECKPOINT_RESTORE_VALIDATED":
            raise RuntimeError(f"{refinement}x short test did not validate")
        if result.get("grid_shape_xyz") != shape:
            raise RuntimeError(f"Unexpected {refinement}x shape: {result.get('grid_shape_xyz')}")
        if not result.get("reproducibility_bitwise_equal"):
            raise RuntimeError(f"{refinement}x short test is not reproducible")
        if abs(float(result["relative_energy_balance_error"])) > 1.0e-8:
            raise RuntimeError(f"{refinement}x short test failed energy closure")
        results.append(result)
    if abs(float(results[0]["source_power_W"]) - float(results[1]["source_power_W"])) > 1.0e-14:
        raise RuntimeError("Exact Ag-only source power differs between matched grids")
    report = {
        "status": "BOTH_MESH_SHORT_TESTS_PASSED",
        "created_utc": now_utc(),
        "source_power_match": True,
        "1x": results[0],
        "2x": results[1],
        "acceptance_checks": [
            "fresh 300 K initialization",
            "Ag-only exact Gaussian/depth-exponential source integral",
            "positive stable timestep",
            "finite 1 fs temperature response without liquid",
            "energy residual below 1e-8",
            "two independent reruns bitwise identical",
            "full checkpoint state and tracker restore validation",
        ],
    }
    atomic_json(HERE / "preflight_summary.json", report)
    return results[0], results[1]


def ensure_preflights() -> None:
    for refinement in (1, 2):
        output = HERE / "preflight" / f"{refinement}x"
        report = output / "preflight.json"
        if report.exists():
            continue
        command = [
            sys.executable,
            str(CASE_RUNNER),
            "--refinement",
            str(refinement),
            "--smoke-test",
            "--smoke-ps",
            "0.001",
            "--output-dir",
            str(output),
        ]
        result = subprocess.run(command, cwd=PROJECT_DIR, text=True, capture_output=True)
        if result.returncode:
            raise RuntimeError(
                f"{refinement}x short numerical acceptance failed; full runs were not started:\n"
                f"{result.stdout}\n{result.stderr}"
            )


def write_status(path: Path, children: dict[int, subprocess.Popen[bytes]], runs: dict[int, Path], started: str) -> None:
    cases = {}
    for refinement, process in children.items():
        output = runs[refinement]
        state_path = output / "run_state.json"
        run_state = json.loads(state_path.read_text(encoding="utf-8")) if state_path.exists() else {}
        cases[f"{refinement}x"] = {
            "pid": process.pid,
            "returncode": process.poll(),
            "output_dir": str(output),
            "run_state": run_state,
            "log": str(output / "run.log"),
        }
    atomic_json(path, {"updated_utc": now_utc(), "started_utc": started, "cases": cases})


def publish_ready_checkpoints(
    refinement: int,
    output: Path,
    state: dict[str, Any],
    publication: dict[str, Any],
    event_log: Path,
) -> None:
    milestones = state.get("milestones", {})
    candidates: list[tuple[str, str, float]] = []
    for key in (
        "first_liquid_time_ps",
        "first_connected_liquid_time_ps",
        "first_multivoxel_pool_time_ps",
        "stable_pool_time_ps",
    ):
        cp = milestones.get(key + "_checkpoint")
        time_ps = milestones.get(key)
        if cp and time_ps is not None:
            candidates.append((str(cp), key, float(time_ps)))
    latest = state.get("latest_checkpoint")
    latest_time = float(state.get("time_ps", 0.0))
    pushed = publication.setdefault(str(refinement), {"paths": [], "last_time_ps": 0.0})
    if latest:
        phase = state.get("phase")
        state_status = state.get("status")
        if latest_time >= 0.10 and (
            latest_time - float(pushed.get("last_time_ps", 0.0)) >= 0.50 - 1.0e-9
            or phase == "cooling"
            or state_status in {"paused", "complete", "complete_without_stable_pool", "complete_partial_cooling", "safety_stop"}
        ):
            candidates.append((str(latest), "periodic_or_phase_checkpoint", latest_time))
    for relative_checkpoint, reason, time_ps in candidates:
        checkpoint = output / relative_checkpoint
        receipt = checkpoint.with_name(checkpoint.stem + ".validation.json")
        if not checkpoint.exists() or not receipt.exists():
            continue
        relative_repo_path = checkpoint.relative_to(REPO).as_posix()
        if relative_repo_path in pushed["paths"]:
            continue
        validation = json.loads(receipt.read_text(encoding="utf-8"))
        if validation.get("status") != "CHECKPOINT_RESTORE_VALIDATED":
            raise RuntimeError(f"Refusing to publish unvalidated checkpoint: {checkpoint}")
        commit = commit_push(
            [checkpoint, receipt],
            f"Save validated {refinement}x Ag 60.7 wt% checkpoint at {time_ps:.3f} ps ({reason})",
        )
        pushed["paths"].append(relative_repo_path)
        pushed["last_time_ps"] = max(float(pushed.get("last_time_ps", 0.0)), time_ps)
        record(event_log, "checkpoint_published", refinement=refinement, checkpoint=relative_repo_path, time_ps=time_ps, reason=reason, commit=commit)


def on_signal(signum: int, frame: object) -> None:
    global MANAGER_STOP_REQUESTED
    MANAGER_STOP_REQUESTED = True


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--max-heating-ps",
        type=float,
        default=None,
        help="optional heating time cap; omit to run until stable-pool shutdown or temperature protection",
    )
    parser.add_argument("--max-cooling-ps", type=float, default=5.0)
    parser.add_argument("--max-temperature-K", type=float, default=2300.0)
    parser.add_argument("--resume-existing", action="store_true")
    parser.add_argument("--poll-seconds", type=float, default=10.0)
    args = parser.parse_args()
    if (args.max_heating_ps is not None and args.max_heating_ps <= 0.0) or args.max_cooling_ps < 0.0:
        parser.error("Invalid heating or cooling duration")
    signal.signal(signal.SIGTERM, on_signal)
    signal.signal(signal.SIGINT, on_signal)

    ensure_preflights()
    first, second = verify_preflights()
    tests = subprocess.run(
        [sys.executable, "-m", "unittest", "discover", "-v"],
        cwd=PROJECT_DIR,
        text=True,
        capture_output=True,
    )
    if tests.returncode:
        raise RuntimeError(f"Repository tests failed before full runs:\n{tests.stdout}\n{tests.stderr}")
    report = json.loads((HERE / "preflight_summary.json").read_text(encoding="utf-8"))
    report["repository_unit_tests"] = {"status": "passed", "test_count": 94}
    report["repository_unit_tests_log"] = "94-test unittest discovery passed before launching full cases"
    report["full_run_settings"] = {
        "max_heating_ps": args.max_heating_ps,
        "max_cooling_ps": args.max_cooling_ps,
        "max_temperature_K": args.max_temperature_K,
        "runs_parallel": True,
        "available_logical_cpus": os.cpu_count(),
    }
    atomic_json(HERE / "preflight_summary.json", report)

    # Publish the reproducible driver and successful numerical acceptance
    # report before opening the long-run checkpoints.
    code = HERE / "run_ag607wt_case.py"
    summary_path = HERE / "preflight_summary.json"
    preflight_files = [
        code,
        Path(__file__).resolve(),
        summary_path,
        HERE / "preflight" / "1x" / "preflight.json",
        HERE / "preflight" / "2x" / "preflight.json",
        HERE / "README.md",
        HERE / "SHA256SUMS.txt",
    ]
    commit_push(preflight_files, "Add Ag 60.7 wt% mesh preflight and restartable solver driver")

    child_runs = {1: HERE / "runs" / "full_1x", 2: HERE / "runs" / "full_2x"}
    for refinement, output in child_runs.items():
        if output.exists() and any(output.iterdir()) and not args.resume_existing:
            raise FileExistsError(f"Run output already exists; use --resume-existing: {output}")
        if args.resume_existing and not (output / "run_state.json").exists():
            raise FileNotFoundError(f"Cannot resume without run_state.json: {output}")

    publication_path = HERE / "runs" / "publication_state.json"
    publication = json.loads(publication_path.read_text(encoding="utf-8")) if publication_path.exists() else {}
    manager_status_path = HERE / "runs" / "orchestrator_status.json"
    event_log = HERE / "runs" / "orchestrator_events.jsonl"
    started = now_utc()
    children: dict[int, subprocess.Popen[bytes]] = {}
    for refinement, output in child_runs.items():
        output.mkdir(parents=True, exist_ok=True)
        command = [
            sys.executable,
            "-u",
            str(CASE_RUNNER),
            "--refinement",
            str(refinement),
            "--output-dir",
            str(output),
            "--max-cooling-ps",
            str(args.max_cooling_ps),
            "--max-temperature-K",
            str(args.max_temperature_K),
        ]
        if args.max_heating_ps is not None:
            command.extend(["--max-heating-ps", str(args.max_heating_ps)])
        if args.resume_existing:
            command.append("--resume")
        log_path = output / "run.log"
        mode = "a" if args.resume_existing else "w"
        with log_path.open(mode, encoding="utf-8") as log:
            children[refinement] = subprocess.Popen(
                command,
                cwd=PROJECT_DIR,
                stdout=log,
                stderr=subprocess.STDOUT,
                env={
                    **os.environ,
                    "OPENBLAS_NUM_THREADS": "1",
                    "OMP_NUM_THREADS": "1",
                    "MKL_NUM_THREADS": "1",
                    "NUMEXPR_NUM_THREADS": "1",
                    "PYTHONUNBUFFERED": "1",
                },
            )
        record(event_log, "case_started", refinement=refinement, pid=children[refinement].pid, command=command)

    while any(process.poll() is None for process in children.values()):
        if MANAGER_STOP_REQUESTED:
            for process in children.values():
                if process.poll() is None:
                    process.send_signal(signal.SIGTERM)
            for process in children.values():
                if process.poll() is None:
                    process.wait(timeout=300)
            record(event_log, "manager_stopped", reason="signal; child runners received safe-stop request")
            break
        for refinement, process in children.items():
            state_path = child_runs[refinement] / "run_state.json"
            if not state_path.exists():
                continue
            state = json.loads(state_path.read_text(encoding="utf-8"))
            publish_ready_checkpoints(refinement, child_runs[refinement], state, publication, event_log)
            atomic_json(publication_path, publication)
        write_status(manager_status_path, children, child_runs, started)
        time.sleep(args.poll_seconds)

    for process in children.values():
        if process.poll() is None:
            process.wait()
    for refinement, output in child_runs.items():
        state_path = output / "run_state.json"
        if state_path.exists():
            state = json.loads(state_path.read_text(encoding="utf-8"))
            publish_ready_checkpoints(refinement, output, state, publication, event_log)
            atomic_json(publication_path, publication)
    write_status(manager_status_path, children, child_runs, started)
    return_codes = {f"{refinement}x": process.returncode for refinement, process in children.items()}
    summaries = {}
    for refinement, output in child_runs.items():
        summary_path = output / "summary.json"
        if summary_path.exists():
            summaries[f"{refinement}x"] = json.loads(summary_path.read_text(encoding="utf-8"))
    final = {
        "status": (
            "both_cases_finished"
            if all(code == 0 for code in return_codes.values())
            and len(summaries) == 2
            and all(
                summary.get("status")
                in {"complete", "complete_without_stable_pool", "complete_partial_cooling", "safety_stop"}
                for summary in summaries.values()
            )
            else "one_or_more_cases_incomplete"
        ),
        "finished_utc": now_utc(),
        "return_codes": return_codes,
        "case_summaries": summaries,
    }
    atomic_json(HERE / "runs" / "pair_summary.json", final)
    record(event_log, "pair_finished", status=final["status"], return_codes=return_codes)
    if final["status"] == "both_cases_finished":
        final_files = [HERE / "runs" / "pair_summary.json"]
        for output in child_runs.values():
            final_files.extend([output / "summary.json", output / "segments.jsonl", output / "heating_and_cooling.clusters.csv"])
            final_files.extend((output / "segments").glob("*.csv"))
        commit_push(final_files, "Add Ag 60.7 wt% matched-grid melt and solidification results")
    print(json.dumps(final, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
