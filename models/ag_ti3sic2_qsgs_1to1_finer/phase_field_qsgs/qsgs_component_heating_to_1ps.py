#!/usr/bin/env python3
"""Continue the native-voxel QSGS Ag heating pilot in 50 fs stages to 1 ps.

The first 0.1--0.15 ps stage is run separately. This driver waits for that
stage's checkpoint, then resumes sequentially. Intermediate checkpoints are
kept under /tmp to avoid committing a large 3-D field for every 50 fs stage;
the final 1 ps fields and energy ledger are written in the repository.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

PHASE_DIR = Path(__file__).resolve().parent
DIAG_DIR = PHASE_DIR / "nested_grid_diagnostics"
STAGE_RUNNER = PHASE_DIR / "qsgs_component_heating_extend.py"
DEFAULT_FIRST_JSON = DIAG_DIR / "qsgs_component_native_heating_to_0p15ps.json"
DEFAULT_FIRST_FIELDS = DIAG_DIR / "qsgs_component_native_heating_to_0p15ps_fields.npz"
DEFAULT_PROGRESS = DIAG_DIR / "qsgs_component_native_heating_to_1ps_progress.json"
DEFAULT_LOG = DIAG_DIR / "qsgs_component_native_heating_to_1ps_orchestrator.log"
INTERMEDIATE_DIR = Path("/tmp/qsgs_component_heating_to_1ps")
FINAL_STEM = "qsgs_component_native_heating_to_1ps"


def _tag(time_ps: float) -> str:
    return f"{time_ps:.2f}".replace(".", "p")


def _write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    tmp.replace(path)


def _pid_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    return True


def _append_new_samples(target: list[dict], incoming: list[dict]) -> None:
    last_time = float(target[-1]["time_ps"]) if target else -1.0
    for row in incoming:
        if float(row["time_ps"]) > last_time + 1.0e-9:
            target.append(row)
            last_time = float(row["time_ps"])


def _pool_gate(samples: list[dict]) -> dict:
    def _has_any_liquid(row: dict) -> bool:
        return bool(row.get("initial_melting_detected", False)) or float(
            row.get("max_ag_liquid_fraction", row.get("ag_liquid_volume_fraction", 0.0))
        ) > 0.0

    initial_melting = any(_has_any_liquid(row) for row in samples)
    minimum_nodes = 50
    minimum_volume = 1.44e3
    sample_interval_ps = 0.025
    valid_runs: list[list[dict]] = []
    active_run: list[dict] = []
    for row in samples:
        cluster = row.get("liquid_connectivity_at_f_l_ge_0p5", {})
        count = int(cluster.get("largest_cluster_node_count", 0))
        volume = float(cluster.get("largest_cluster_volume_A3", 0.0))
        min_volume = float(row.get("minimum_pool_volume_A3", minimum_volume))
        energy_ok = float(row.get("energy_residual_relative_to_source", 0.0)) <= 1.0e-10
        qualifies = count >= minimum_nodes and volume >= max(minimum_volume, min_volume) and energy_ok
        time = float(row["time_ps"])
        contiguous = bool(active_run) and abs(
            time - float(active_run[-1]["time_ps"]) - sample_interval_ps
        ) <= 1.0e-6
        overlap = int(row.get("largest_cluster_overlap_with_previous_sample", {}).get("shared_node_count", 0))
        same_region = bool(active_run) and overlap > 0
        if qualifies and (not active_run or (contiguous and same_region)):
            active_run.append(row)
        else:
            if active_run:
                valid_runs.append(active_run)
            active_run = [row] if qualifies else []
    if active_run:
        valid_runs.append(active_run)

    def _gate_for_duration(required_ps: float) -> dict:
        qualifying = []
        for run in valid_runs:
            if len(run) < 2:
                continue
            run_duration = float(run[-1]["time_ps"]) - float(run[0]["time_ps"])
            if run_duration + 1.0e-9 < required_ps:
                continue
            end_time = float(run[0]["time_ps"]) + required_ps
            window = [row for row in run if float(row["time_ps"]) <= end_time + 1.0e-9]
            volumes = [
                float(row["liquid_connectivity_at_f_l_ge_0p5"]["largest_cluster_volume_A3"])
                for row in window
            ]
            overlap_rows = window[1:]
            jaccard = [
                float(row["largest_cluster_overlap_with_previous_sample"]["jaccard_fraction"])
                for row in overlap_rows
                if row.get("largest_cluster_overlap_with_previous_sample", {}).get("jaccard_fraction") is not None
            ]
            qualifying.append({
                "start_time_ps": float(window[0]["time_ps"]),
                "first_pass_time_ps": end_time,
                "end_sample_time_ps": float(window[-1]["time_ps"]),
                "observed_duration_ps": float(window[-1]["time_ps"]) - float(window[0]["time_ps"]),
                "sample_count": len(window),
                "minimum_node_count": min(
                    int(row["liquid_connectivity_at_f_l_ge_0p5"]["largest_cluster_node_count"])
                    for row in window
                ),
                "minimum_cluster_volume_A3": min(volumes),
                "maximum_cluster_volume_A3": max(volumes),
                "volume_peak_to_peak_fraction_of_mean": float(
                    (max(volumes) - min(volumes)) / max(float(np.mean(volumes)), 1.0e-30)
                ),
                "minimum_adjacent_jaccard_fraction": min(jaccard) if jaccard else None,
                "latest_cluster_metrics": window[-1]["liquid_connectivity_at_f_l_ge_0p5"],
            })
        return {"passed": bool(qualifying), "first_qualifying_window": qualifying[0] if qualifying else None}

    short_005 = _gate_for_duration(0.05)
    short_03 = _gate_for_duration(0.3)
    enhanced_05 = _gate_for_duration(0.5)
    return {
        "initial_melting_detected": initial_melting,
        "initial_melting_sample_times_ps": [
            float(row["time_ps"]) for row in samples if _has_any_liquid(row)
        ],
        "preliminary_pool_gate_met": short_005["passed"],
        "preliminary_local_pool_gate": {
            "liquid_threshold_fl": 0.5,
            "minimum_nodes": minimum_nodes,
            "approximate_minimum_volume_A3": minimum_volume,
            "consecutive_samples": 3,
            "sample_interval_ps": sample_interval_ps,
            "required_adjacent_overlap_nodes": ">0",
            **short_005,
        },
        "short_time_persistence_0p3ps": short_03,
        "enhanced_persistence_0p5ps": enhanced_05,
        "candidate_runs": [
            {
                "start_time_ps": float(run[0]["time_ps"]),
                "end_time_ps": float(run[-1]["time_ps"]),
                "duration_ps": float(run[-1]["time_ps"]) - float(run[0]["time_ps"]),
                "sample_count": len(run),
            }
            for run in valid_runs
        ],
        "locality_metrics": "For a qualifying cluster, inspect source-center distance, periodic bbox/extents and Ag-volume fraction; no extra numeric locality cutoff was supplied.",
    }


def _wait_for_first_stage(json_path: Path, fields_path: Path, pid: int | None, log) -> dict:
    while not json_path.exists():
        if pid is not None and not _pid_alive(pid):
            raise RuntimeError(
                f"The active 0.1--0.15 ps process (PID {pid}) exited before writing {json_path}"
            )
        log.write(f"Waiting for active 0.1--0.15 ps stage; pid={pid}\n")
        log.flush()
        time.sleep(60)
    result = json.loads(json_path.read_text(encoding="utf-8"))
    if result.get("status") != "NATIVE_QSGS_AG_HEATING_PILOT_COMPLETED":
        raise RuntimeError(f"Unexpected first-stage status in {json_path}")
    if abs(float(result["heating"]["end_time_ps"]) - 0.15) > 1.0e-9:
        raise RuntimeError("First stage did not end at 0.15 ps")
    if not fields_path.exists():
        raise FileNotFoundError(fields_path)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first-json", type=Path, default=DEFAULT_FIRST_JSON)
    parser.add_argument("--first-fields", type=Path, default=DEFAULT_FIRST_FIELDS)
    parser.add_argument("--wait-pid", type=int, default=None)
    parser.add_argument("--start-time-ps", type=float, default=0.15)
    parser.add_argument("--stop-time-ps", type=float, default=1.0)
    parser.add_argument("--stage-ps", type=float, default=0.05)
    parser.add_argument("--dt-fs", type=float, default=0.3125)
    parser.add_argument("--progress", type=Path, default=DEFAULT_PROGRESS)
    parser.add_argument("--log", type=Path, default=DEFAULT_LOG)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    if args.stop_time_ps > 1.0 + 1.0e-12:
        parser.error("This task is capped at 1 ps")
    if args.stop_time_ps <= args.start_time_ps or args.stage_ps <= 0 or args.dt_fs <= 0:
        parser.error("Invalid time interval or step")
    interval_steps = args.stage_ps * 1000.0 / args.dt_fs
    if abs(interval_steps - round(interval_steps)) > 1.0e-9:
        parser.error("Stage duration must be an integer multiple of dt")

    stages = []
    current = args.start_time_ps
    while current < args.stop_time_ps - 1.0e-12:
        end = min(args.stop_time_ps, current + args.stage_ps)
        stages.append((round(current, 10), round(end, 10)))
        current = end
    if args.dry_run:
        print(json.dumps({"stages": stages, "dt_fs": args.dt_fs}, indent=2))
        return

    args.first_json = args.first_json.resolve()
    args.first_fields = args.first_fields.resolve()
    args.start_time_ps = round(args.start_time_ps, 10)
    args.log.parent.mkdir(parents=True, exist_ok=True)
    args.progress.parent.mkdir(parents=True, exist_ok=True)
    INTERMEDIATE_DIR.mkdir(parents=True, exist_ok=True)
    state = {
        "status": "WAITING_FOR_0P15PS_CHECKPOINT",
        "start_time_ps": args.start_time_ps,
        "stop_time_ps": args.stop_time_ps,
        "stage_ps": args.stage_ps,
        "dt_fs": args.dt_fs,
        "completed_stages": [],
        "samples": [],
        "stable_pool_gate": {
            "project_candidate_rule": "largest connected liquid Ag cluster persists >=3 ps; volume variation <=5% over last 2 ps",
            "evaluable_before_1ps": False,
        },
    }
    _write_json(args.progress, state)

    with args.log.open("a", encoding="utf-8", buffering=1) as log:
        log.write("Waiting for the already-running 0.1--0.15 ps stage to finish.\n")
        first_result = _wait_for_first_stage(args.first_json, args.first_fields, args.wait_pid, log)
        current_json = args.first_json
        current_fields = args.first_fields
        current_time = args.start_time_ps
        if abs(float(first_result["heating"]["end_time_ps"]) - current_time) > 1.0e-9:
            raise RuntimeError("Initial stage end time does not match requested start")
        state["status"] = "HEATING_IN_STAGES"
        state["first_stage_end"] = first_result["samples"][-1]
        _append_new_samples(state["samples"], first_result["samples"])
        state["melting_and_pool_assessment"] = _pool_gate(state["samples"])
        _write_json(args.progress, state)
        log.write(f"First stage complete at {current_time:.2f} ps; continuing to {args.stop_time_ps:.2f} ps.\n")

        env = os.environ.copy()
        env.update({"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"})
        for start, end in stages:
            if abs(start - current_time) > 1.0e-9:
                raise RuntimeError(f"Stage chain discontinuity: expected {current_time}, got {start}")
            final_stage = abs(end - args.stop_time_ps) <= 1.0e-9
            if final_stage:
                output = DIAG_DIR / f"{FINAL_STEM}.json"
                stage_log = DIAG_DIR / f"{FINAL_STEM}.log"
                driver_log = DIAG_DIR / f"{FINAL_STEM}_driver.log"
            else:
                tag = _tag(end)
                output = INTERMEDIATE_DIR / f"heating_to_{tag}ps.json"
                stage_log = INTERMEDIATE_DIR / f"heating_to_{tag}ps.log"
                driver_log = INTERMEDIATE_DIR / f"heating_to_{tag}ps_driver.log"
            command = [
                sys.executable, str(STAGE_RUNNER),
                "--checkpoint", str(current_fields),
                "--checkpoint-json", str(current_json),
                "--checkpoint-time-ps", f"{start:.10g}",
                "--end-time-ps", f"{end:.10g}",
                "--dt-fs", f"{args.dt_fs:.10g}",
                "--output", str(output),
                "--log", str(stage_log),
            ]
            log.write(f"Launching stage {start:.2f}->{end:.2f} ps.\n")
            with driver_log.open("w", encoding="utf-8") as driver:
                completed = subprocess.run(command, cwd=PHASE_DIR.parents[2], env=env, stdout=driver, stderr=subprocess.STDOUT)
            if completed.returncode != 0:
                raise RuntimeError(f"Stage {start:.2f}->{end:.2f} failed; see {driver_log}")
            result = json.loads(output.read_text(encoding="utf-8"))
            if result.get("status") != "NATIVE_QSGS_AG_HEATING_PILOT_COMPLETED":
                raise RuntimeError(f"Stage {start:.2f}->{end:.2f} did not report completion")
            final = result["samples"][-1]
            if abs(float(final["time_ps"]) - end) > 1.0e-8:
                raise RuntimeError(f"Stage ended at {final['time_ps']}, expected {end}")
            if float(final["energy_residual_relative_to_source"]) > 1.0e-10:
                raise RuntimeError(f"Energy balance failed at {end} ps")
            if float(final["electron_temperature_max_K"]) > 20000.0:
                raise RuntimeError(f"Electron temperature exceeded table limit at {end} ps")
            fields = output.with_name(f"{output.stem}_fields.npz")
            if not fields.exists():
                raise FileNotFoundError(fields)

            record = {
                "start_time_ps": start,
                "end_time_ps": end,
                "status": "completed",
                "max_electron_temperature_K": final["electron_temperature_max_K"],
                "max_lattice_temperature_K": final["lattice_temperature_max_K"],
                "ag_liquid_volume_fraction": final["ag_liquid_volume_fraction"],
                "largest_liquid_cluster_volume_A3": final["liquid_connectivity_at_f_l_ge_0p5"]["largest_cluster_volume_A3"],
                "energy_residual_relative_to_source": final["energy_residual_relative_to_source"],
            }
            state["completed_stages"].append(record)
            _append_new_samples(state["samples"], result["samples"])
            state["melting_and_pool_assessment"] = _pool_gate(state["samples"])
            state["current_time_ps"] = end
            state["latest_sample"] = final
            state["status"] = "COMPLETED_AT_1PS" if final_stage else "HEATING_IN_STAGES"
            _write_json(args.progress, state)
            log.write(json.dumps(record, ensure_ascii=False) + "\n")
            log.write(f"Completed {end:.2f} ps: Tl,max={final['lattice_temperature_max_K']:.6f} K; liquid={final['ag_liquid_volume_fraction']:.8g}.\n")

            # Keep only the latest temporary restart. The 0.1 ps and final 1 ps
            # snapshots remain as durable repository checkpoints.
            if current_fields.parent == INTERMEDIATE_DIR and current_fields != fields:
                current_fields.unlink(missing_ok=True)
            if current_json.parent == INTERMEDIATE_DIR and current_json != output:
                current_json.unlink(missing_ok=True)
            if current_fields.parent == INTERMEDIATE_DIR and current_fields != fields:
                current_fields.with_name(current_fields.stem.removesuffix("_fields") + ".log").unlink(missing_ok=True)
            current_json, current_fields, current_time = output, fields, end

        state["task_result"] = {
            "stop_reason": "1 ps hard limit reached",
            "heat_source_on_at_stop": True,
            "liquid_observed": state["melting_and_pool_assessment"]["initial_melting_detected"],
            "preliminary_pool_size_and_persistence_gate_met": state["melting_and_pool_assessment"]["preliminary_pool_gate_met"],
            "stable_pool_claimed": False,
            "three_stage_classification": {
                "initial_melting": state["melting_and_pool_assessment"]["initial_melting_detected"],
                "preliminary_local_pool_size_and_persistence": state["melting_and_pool_assessment"]["preliminary_pool_gate_met"],
                "short_time_persistence_0p3ps": state["melting_and_pool_assessment"]["short_time_persistence_0p3ps"]["passed"],
                "enhanced_persistence_0p5ps": state["melting_and_pool_assessment"]["enhanced_persistence_0p5ps"]["passed"],
                "stable_pool": False,
                "locality_review": "Use saved cluster bbox/extents, distance to heat-source center, Ag-volume fraction and source-power fraction; no unsupported numeric locality cutoff is imposed.",
            },
            "note": "The existing 3 ps persistence / 2 ps volume-window candidate criterion cannot be demonstrated within a 1 ps cap.",
        }
        _write_json(args.progress, state)
        log.write("Reached the 1 ps task limit; no further heating was scheduled.\n")


if __name__ == "__main__":
    main()
