#!/usr/bin/env python3
"""Run the existing heating driver with wall-time/RSS and optional cProfile."""

from __future__ import annotations

import argparse
import cProfile
import json
import pstats
import resource
import runpy
import sys
import time
from pathlib import Path


def main() -> None:
    if "--" not in sys.argv:
        raise SystemExit("usage: profile_ag70_short_run.py --measurement FILE [--profile FILE] -- DRIVER.py DRIVER_ARGS...")
    divider = sys.argv.index("--")
    parser = argparse.ArgumentParser()
    parser.add_argument("--measurement", type=Path, required=True)
    parser.add_argument("--profile", type=Path)
    runner_args = parser.parse_args(sys.argv[1:divider])
    child_args = sys.argv[divider + 1:]
    if not child_args:
        raise SystemExit("missing heating driver command")

    driver = Path(child_args[0]).resolve()
    sys.argv = child_args
    profiler = cProfile.Profile() if runner_args.profile else None
    started = time.perf_counter()
    try:
        if profiler:
            profiler.enable()
        runpy.run_path(str(driver), run_name="__main__")
    finally:
        if profiler:
            profiler.disable()
        elapsed = time.perf_counter() - started
        usage = resource.getrusage(resource.RUSAGE_SELF)
        measurement = {
            "driver": str(driver),
            "driver_arguments": child_args[1:],
            "elapsed_wall_seconds": elapsed,
            "max_rss_kib": usage.ru_maxrss,
            "profiled": profiler is not None,
        }
        runner_args.measurement.parent.mkdir(parents=True, exist_ok=True)
        runner_args.measurement.write_text(json.dumps(measurement, indent=2) + "\n", encoding="utf-8")
        if profiler:
            assert runner_args.profile is not None
            runner_args.profile.parent.mkdir(parents=True, exist_ok=True)
            profiler.dump_stats(str(runner_args.profile))
            summary_path = runner_args.profile.with_suffix(".top_functions.txt")
            with summary_path.open("w", encoding="utf-8") as handle:
                pstats.Stats(profiler, stream=handle).strip_dirs().sort_stats("cumulative").print_stats(40)
        print(json.dumps(measurement, indent=2), flush=True)


if __name__ == "__main__":
    main()
