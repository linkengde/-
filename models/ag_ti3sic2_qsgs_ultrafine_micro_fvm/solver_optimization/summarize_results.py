#!/usr/bin/env python3
"""Compile run summaries and step logs into one comparison CSV."""

from __future__ import annotations

import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
CHECKPOINTS = HERE / "results" / "checkpoints"
OUT = HERE / "results" / "OPTIMIZATION_COMPARISON.csv"


def read_jsonl(path: Path):
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def compile_rows():
    rows = []
    for folder in sorted(CHECKPOINTS.iterdir()):
        if not folder.is_dir() or folder.name == "shared_inputs":
            continue
        summary_path = folder / "run_summary.json"
        if not summary_path.exists():
            continue
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        configuration = summary.get("run_configuration", {})
        steps = [item.get("accepted_step", {}) for item in read_jsonl(folder / "step_log.jsonl")]
        hot = summary.get("heated_endpoint_metrics") or {}
        solver_step_wall = sum(float(item.get("step_wall_seconds", 0.0)) for item in steps)
        row = {
            "run_id": folder.name,
            "status": summary.get("status"),
            "case": summary.get("case"),
            "n": summary.get("resolution"),
            "solver": summary.get("solver"),
            "linear_policy": summary.get("linear_policy"),
            "preconditioner": summary.get("preconditioner", "fft"),
            "dt_factor": configuration.get("dt_factor", ""),
            "domain_edge_m": configuration.get("domain_edge_m", ""),
            "source_rate_W_m3_on_Ag_volume": configuration.get("source_rate_W_m3_on_Ag_volume", ""),
            "threads": summary.get("threads"),
            "accepted_steps": len(steps),
            "retries": summary.get("retry_count_this_invocation"),
            "newton_iterations": sum(int(item.get("nonlinear_iterations", 0)) for item in steps),
            "pcg_iterations": sum(int(item.get("pcg_iterations_total", 0)) for item in steps),
            "solver_step_wall_seconds": solver_step_wall,
            "wall_seconds": summary.get("wall_seconds_total_estimated", summary.get("wall_seconds_this_invocation")),
            "cpu_seconds": summary.get("cpu_seconds_total_estimated", summary.get("cpu_seconds_this_invocation")),
            "whole_run_timing_reliable": summary.get("whole_run_timing_reliable", True),
            "average_cpu_cores": summary.get("average_cpu_cores_used"),
            "peak_rss_mib": summary.get("peak_rss_mib"),
            "physical_time_s": summary.get("physical_time_s"),
            "phase_elapsed_s": summary.get("phase_elapsed_s"),
            "max_temperature_K": max((float(item.get("max_temperature_K", 0.0)) for item in steps), default=0.0),
            "ag_max_liquid_fraction_hot": hot.get("Ag_max_liquid_fraction"),
            "ag_bulk_liquid_fraction_hot": hot.get("Ag_bulk_liquid_fraction"),
            "ag_liquid_weighted_volume_m3_hot": hot.get("Ag_liquid_weighted_volume_m3"),
            "normalized_energy_residual": summary.get("normalized_energy_residual"),
            "jit_warmup_seconds": summary.get("jit_warmup_seconds"),
            "amg_setup_seconds": "",
            "amg_levels": "",
            "notes": summary.get("timing_note", "step wall sums accepted-step solver time; whole-run wall includes checkpoint/log IO"),
        }
        rows.append(row)

    baseline_by_n = {
        row["n"]: row["solver_step_wall_seconds"]
        for row in rows if row["solver"] == "baseline" and row["status"] == "COMPLETE"
    }
    for row in rows:
        baseline = baseline_by_n.get(row["n"])
        row["solver_step_speedup_vs_completed_baseline"] = (
            baseline / row["solver_step_wall_seconds"]
            if baseline and row["solver_step_wall_seconds"] > 0 else ""
        )

    amg_path = HERE / "results" / "AMG_PROBE.json"
    if amg_path.exists():
        for item in json.loads(amg_path.read_text(encoding="utf-8")):
            rows.append({
                "run_id": f"n50_hot_single_step_{item['preconditioner']}",
                "status": "PROBE_COMPLETE", "case": "hot_state_numerical_pulse",
                "n": item["resolution"], "solver": "numba", "linear_policy": item["linear_policy"],
                "preconditioner": item["preconditioner"], "threads": 4,
                "dt_factor": "hot-state probe", "domain_edge_m": "", "source_rate_W_m3_on_Ag_volume": "",
                "accepted_steps": 1, "retries": 0,
                "newton_iterations": item["newton_iterations"], "pcg_iterations": item["pcg_iterations"],
                "solver_step_wall_seconds": item["step_wall_seconds"],
                "wall_seconds": item["step_wall_seconds"], "cpu_seconds": "",
                "whole_run_timing_reliable": "single-step probe",
                "average_cpu_cores": "", "peak_rss_mib": "",
                "physical_time_s": item["dt_s"], "phase_elapsed_s": item["dt_s"],
                "max_temperature_K": item["max_temperature_K"],
                "ag_max_liquid_fraction_hot": item["max_ag_liquid_fraction"],
                "ag_bulk_liquid_fraction_hot": item["ag_bulk_liquid_fraction"],
                "ag_liquid_weighted_volume_m3_hot": "",
                "normalized_energy_residual": item["normalized_energy_residual"],
                "jit_warmup_seconds": item.get("jit_warmup_seconds_excluded"),
                "amg_setup_seconds": item.get("amg_setup_seconds_sum", 0.0),
                "amg_levels": item.get("amg_levels", 0),
                "solver_step_speedup_vs_completed_baseline": "",
                "notes": "single preheated-state step; not a complete thermal cycle",
            })
    return rows


def main():
    rows = compile_rows()
    if not rows:
        raise SystemExit("No run summaries found")
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with OUT.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} rows to {OUT}")


if __name__ == "__main__":
    main()
