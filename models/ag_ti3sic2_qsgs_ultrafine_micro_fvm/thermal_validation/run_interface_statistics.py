#!/usr/bin/env python3
"""顺序刷新界面法向热流统计，保留在各算例JSON中。"""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import time


HERE = Path(__file__).resolve().parent
RUNNER = HERE / "run_conductivity_case.py"
CASES = []
for size in ("ultrafine_500nm", "micro_10um"):
    for seed in (1105, 1136, 1182):
        CASES.append((size, seed, 0.0))
        CASES.append((size, seed, 1.0e-9))
    for resistance in (1.0e-10, 1.0e-8):
        CASES.append((size, 1105, resistance))


def main() -> None:
    total_start = time.perf_counter()
    env = os.environ.copy()
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    for index, (size, seed, resistance) in enumerate(CASES, start=1):
        command = [
            sys.executable, str(RUNNER), "--grid-n", "150", "--seed", str(seed),
            "--size-name", size, "--interface-resistance", f"{resistance:.17g}",
            "--kind", "qsgs_morphology",
        ]
        start = time.perf_counter()
        result = subprocess.run(command, check=False, capture_output=True, text=True, env=env)
        if result.stdout:
            print(f"{index}/{len(CASES)} {result.stdout.rstrip()}", flush=True)
        if result.returncode:
            if result.stderr:
                print(result.stderr, file=sys.stderr, flush=True)
            raise SystemExit(result.returncode)
        print(f"case_seconds={time.perf_counter()-start:.2f}", flush=True)
    print(f"completed_interface_cases={len(CASES)} total_seconds={time.perf_counter()-total_start:.2f}", flush=True)


if __name__ == "__main__":
    main()
