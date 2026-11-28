#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""汇总10 ps/5 ps局部凝固事件细化与完整循环的凝固时间差。"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"


def read(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(path)
    return json.loads(path.read_text(encoding="utf-8"))


def relative(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), 1e-30)


def main() -> None:
    base = read(RESULTS / "pankratz_v2_n100_dt1" / "summary.json")
    half = read(RESULTS / "pankratz_v2_n100_dt0.5" / "summary.json")
    quarter = read(RESULTS / "pankratz_v2_n100_dt0.25" / "summary.json")
    event10 = read(RESULTS / "pankratz_v2_n100_freezing_event_refinement" / "summary.json")
    event5 = read(RESULTS / "pankratz_v2_n100_freezing_event_cap_5ps" / "summary.json")
    t10 = float(event10["cooling_time_to_complete_freeze_s"])
    t5 = float(event5["cooling_time_to_complete_freeze_s"])
    payload = {
        "cases": {"uniform_dt_1": base, "uniform_dt_0.5": half, "uniform_dt_0.25": quarter,
                  "liquid_event_cap_10ps": event10, "liquid_event_cap_5ps": event5},
        "relative_change_5ps_vs_10ps_cooling_freeze_time": relative(t10, t5),
        "relative_change_event_5ps_vs_uniform_dt_0.25": relative(
            float(quarter["cooling_time_to_complete_freeze_s"]), t5),
        "relative_change_event_5ps_vs_uniform_dt_1": relative(
            float(base["cooling_time_to_complete_freeze_s"]), t5),
        "event_time_convergence_3percent": relative(t10, t5) <= 0.03,
        "interpretation": (
            "10 ps and 5 ps capped cooling cases start from the identical saved n=100 heated field; "
            "the source energy and initial enthalpy reference are carried into the cooling-only restart. "
            "This is local event refinement; it does not change the result that uniform dt factors 1, 0.5, "
            "and 0.25 do not meet the complete-freeze-time 3% pairwise target."
        ),
    }
    out = RESULTS / "FREEZING_EVENT_REFINEMENT.json"
    out.write_text(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                   encoding="utf-8")
    print(json.dumps({k: v for k, v in payload.items() if k != "cases"},
                     ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
