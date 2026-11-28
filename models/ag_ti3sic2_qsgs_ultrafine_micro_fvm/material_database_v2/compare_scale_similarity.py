#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""比较10 μm与200 μm同无量纲热源计算的场与热循环时间尺度。"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
BASE = RESULTS / "same_load_50" / "pankratz_v2_n50_heated_state.npz"
SCALED = RESULTS / "scaled_n50_micro200um_similarity" / "heated_state.npz"


def main() -> None:
    with np.load(BASE, allow_pickle=False) as base, np.load(SCALED, allow_pickle=False) as scaled:
        metrics = {}
        for key in ("temperature_K", "ag_fraction", "liquid_fraction"):
            a = np.asarray(base[key], dtype=np.float64)
            b = np.asarray(scaled[key], dtype=np.float64)
            metrics[key] = {
                "max_absolute_difference": float(np.max(np.abs(a - b))),
                "max_relative_difference_to_baseline": float(
                    np.max(np.abs(a - b)) / max(float(np.max(np.abs(a))), 1e-300)),
            }
        base_time = float(base["time_s"])
        scaled_time = float(scaled["time_s"])
        edge_base = float(np.asarray(base["spacing_m"])[0] * base["temperature_K"].shape[0])
        edge_scaled = float(np.asarray(scaled["spacing_m"])[0] * scaled["temperature_K"].shape[0])
    reference = json.loads((RESULTS / "same_load_50" / "pankratz_v2_summary.json").read_text())
    large = json.loads((RESULTS / "scaled_n50_micro200um_similarity" / "summary.json").read_text())
    spatial_scale = edge_scaled / edge_base
    temporal_expected = spatial_scale**2
    ratios = {
        "first_melt_time": large["first_melt_time_s"] / reference["first_melt_time_s"],
        "complete_freeze_cooling_time": (
            large["cooling_time_to_complete_freeze_s"]
            / reference["cooling_time_to_complete_freeze_s"]),
        "heated_snapshot_time": scaled_time / base_time,
    }
    result = {
        "status": "COMPLETE",
        "base_edge_m": edge_base,
        "scaled_edge_m": edge_scaled,
        "length_scale": spatial_scale,
        "expected_time_scale_L_squared": temporal_expected,
        "observed_time_ratios": ratios,
        "field_comparison_at_heating_end": metrics,
        "base_summary": "same_load_50/pankratz_v2_summary.json",
        "scaled_summary": "scaled_n50_micro200um_similarity/summary.json",
        "interpretation": "Numerical dimensionless-similarity check only; not equal physical arc loading or contact performance.",
    }
    (RESULTS / "SCALE_SIMILARITY_COMPARISON.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
