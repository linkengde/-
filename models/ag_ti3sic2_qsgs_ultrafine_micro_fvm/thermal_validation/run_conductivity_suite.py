#!/usr/bin/env python3
"""顺序运行资源受控的性能、固定几何、形貌和界面热阻算例。"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time


HERE = Path(__file__).resolve().parent
CASE_SCRIPT = HERE / "run_conductivity_case.py"
RESULTS = HERE / "results" / "cases"
SIZES = ("ultrafine_500nm", "micro_10um")
SEEDS = (1105, 1136, 1182)


def case_key(grid_n: int, seed: int, size: str, resistance: float, kind: str) -> str:
    return f"{kind}_n{grid_n}_seed{seed}_{size}_R{resistance:.0e}"


def run_case(
    grid_n: int,
    seed: int,
    size: str,
    resistance: float,
    kind: str,
    snapshot: bool = False,
) -> None:
    key = case_key(grid_n, seed, size, resistance, kind)
    output = RESULTS / f"{key}.json"
    if output.is_file():
        print(json.dumps({"status": "reused", "case_id": key}), flush=True)
        return
    cmd = [
        sys.executable, str(CASE_SCRIPT),
        "--grid-n", str(grid_n), "--seed", str(seed),
        "--size-name", size, "--interface-resistance", f"{resistance:.17g}",
        "--kind", kind,
    ]
    if snapshot:
        cmd.append("--snapshot")
    env = os.environ.copy()
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    started = time.perf_counter()
    completed = subprocess.run(cmd, check=False, capture_output=True, text=True, env=env)
    if completed.stdout:
        print(completed.stdout.rstrip(), flush=True)
    if completed.returncode != 0:
        if completed.stderr:
            print(completed.stderr.rstrip(), file=sys.stderr, flush=True)
        raise RuntimeError(f"case failed: {key}; exit={completed.returncode}")
    if not output.is_file():
        raise RuntimeError(f"case completed without result file: {key}")
    print(json.dumps({"status": "saved", "case_id": key,
                      "suite_wall_seconds": time.perf_counter() - started}), flush=True)


def performance() -> None:
    for n in (50, 100, 150):
        run_case(n, 1105, "ultrafine_500nm", 0.0, "performance", snapshot=(n == 150))


def fixed_geometry_convergence() -> None:
    for resistance in (0.0, 1.0e-9):
        run_case(150, 1105, "ultrafine_500nm", resistance, "qsgs_morphology")
        run_case(300, 1105, "ultrafine_500nm", resistance, "fixed_nested_300")


def qsgs_morphology() -> None:
    for n in (150, 200, 300):
        for seed in SEEDS:
            if n == 150 and seed == 1105:
                # Reuse the measured performance geometry's identical R=0 case.
                continue
            run_case(n, seed, "ultrafine_500nm", 0.0, "qsgs_morphology")


def scale_invariance() -> None:
    # Same binary mask, 10 μm and 200 μm, R=0; the second solve checks unit scaling.
    run_case(150, 1105, "micro_10um", 0.0, "qsgs_morphology")


def finite_interface() -> None:
    # R=1e-9 m² K/W is a numerical sensitivity value, not measured material data.
    for size in SIZES:
        for seed in SEEDS:
            run_case(150, seed, size, 1.0e-9, "qsgs_morphology")
    # Three-point sensitivity curve for one representative seed at both physical scales.
    for size in SIZES:
        for resistance in (1.0e-10, 1.0e-8):
            run_case(150, 1105, size, resistance, "qsgs_morphology")


PHASES = {
    "performance": performance,
    "fixed-convergence": fixed_geometry_convergence,
    "qsgs-morphology": qsgs_morphology,
    "scale-invariance": scale_invariance,
    "finite-interface": finite_interface,
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", choices=(*PHASES, "all"), default="all")
    args = parser.parse_args()
    phases = tuple(PHASES) if args.phase == "all" else (args.phase,)
    suite_started = time.perf_counter()
    for phase in phases:
        print(json.dumps({"phase": phase, "status": "starting"}), flush=True)
        PHASES[phase]()
        print(json.dumps({"phase": phase, "status": "completed"}), flush=True)
    summary = {"status": "SUITE_PHASES_COMPLETED", "phases": list(phases),
               "wall_seconds": time.perf_counter() - suite_started}
    (HERE / "results").mkdir(parents=True, exist_ok=True)
    (HERE / "results" / "suite_run.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
