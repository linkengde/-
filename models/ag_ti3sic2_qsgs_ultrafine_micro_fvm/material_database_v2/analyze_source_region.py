#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从已保存温度场分别统计全局Ag熔化与局部热源区Ag熔化。"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy import ndimage

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"


def measure(label: str, field_path: Path, source_path: Path) -> dict:
    with np.load(field_path, allow_pickle=False) as data:
        ag = np.asarray(data["ag_fraction"], dtype=np.float64)
        liquid = np.asarray(data["liquid_fraction"], dtype=np.float64)
        spacing = np.asarray(data["spacing_m"], dtype=np.float64)
        temperature = np.asarray(data["temperature_K"], dtype=np.float64)
        time_s = float(data["time_s"])
    source = np.asarray(np.load(source_path, allow_pickle=False), dtype=np.float64)
    cell_volume = float(np.prod(spacing))
    source_region = source > 0.0
    ag_volume = float(np.sum(ag) * cell_volume)
    liq_volume = float(np.sum(ag * liquid) * cell_volume)
    source_ag_volume = float(np.sum(ag[source_region]) * cell_volume)
    source_liq_volume = float(np.sum((ag * liquid)[source_region]) * cell_volume)
    source_ag_weight = float(np.sum(ag[source_region]))
    source_local_fraction = (float(np.sum((ag * liquid)[source_region])) / source_ag_weight
                             if source_ag_weight else 0.0)
    source_max_fraction = float(np.max(liquid[source_region & (ag > 1e-14)], initial=0.0))
    # “明确液态区”按Ag体积分数≥0.5且液相率≥0.5定义，并用26邻接判断连通。
    # 这只是数值分辨下的连通性指标，不代替实验尺度验证。
    melt_mask = (ag >= 0.5) & (liquid >= 0.5)
    labels, component_count = ndimage.label(melt_mask, structure=np.ones((3, 3, 3), dtype=np.uint8))
    sizes = np.bincount(labels.ravel())[1:] if component_count else np.zeros(0, dtype=np.int64)
    largest_cells = int(sizes.max(initial=0))
    largest_label = int(np.argmax(sizes) + 1) if largest_cells else 0
    largest_mask = labels == largest_label if largest_label else np.zeros_like(melt_mask)
    return {
        "case": label, "field_file": str(field_path.relative_to(HERE)),
        "source_mask_file": str(source_path.relative_to(HERE)),
        "time_s": time_s, "resolution": int(ag.shape[0]),
        "domain_edge_m": float(ag.shape[0] * spacing[0]),
        "Ag_total_volume_m3": ag_volume,
        "Ag_bulk_liquid_fraction": float(np.sum(ag * liquid) / np.sum(ag)),
        "Ag_liquid_volume_m3": liq_volume,
        "Ag_max_temperature_K": float(np.max(temperature[ag > 1e-14], initial=300.0)),
        "Ti3SiC2_max_temperature_K": float(np.max(temperature[(1.0 - ag) > 1e-14], initial=300.0)),
        "global_min_temperature_K": float(temperature.min()),
        "source_region_cell_count": int(np.count_nonzero(source_region)),
        "source_region_Ag_volume_m3": source_ag_volume,
        "source_region_Ag_liquid_fraction": source_local_fraction,
        "source_region_Ag_liquid_volume_m3": source_liq_volume,
        "source_region_Ag_max_liquid_fraction": source_max_fraction,
        "source_region_peak_temperature_K": float(np.max(temperature[source_region], initial=300.0)),
        "melt_threshold_Ag_fraction": 0.5,
        "melt_threshold_liquid_fraction": 0.5,
        "melt_component_connectivity": "26-neighbor",
        "melt_component_count": int(component_count),
        "largest_melt_component_cells": largest_cells,
        "largest_melt_component_geometric_volume_m3": float(largest_cells * cell_volume),
        "largest_melt_component_Ag_liquid_weighted_volume_m3": float(
            np.sum((ag * liquid)[largest_mask]) * cell_volume),
        "melt_component_interpretation": "thresholded control volumes; numerical continuity only",
        "source_power_W": float(np.sum(source) * cell_volume),
        "source_total_energy_J_to_snapshot": float(np.sum(source) * cell_volume * time_s),
    }


def main() -> None:
    cases = []
    same = RESULTS / "same_load_50"
    for tag in ("old_candidate", "pankratz_v2"):
        cases.append(measure(
            f"50³_same_load_{tag}", same / f"{tag}_n50_heated_state.npz",
            same / f"checkpoint_{tag}" / "source_field_W_m3.npy"))
    for factor in ("1", "0.5", "0.25"):
        case_dir = RESULTS / f"pankratz_v2_n100_dt{factor}"
        if not (case_dir / "heated_state.npz").exists():
            continue
        cases.append(measure(
            f"100³_Pankratz_v2_dt{factor}", case_dir / "heated_state.npz",
            case_dir / "checkpoint" / "source_field_W_m3.npy"))
    probe150 = RESULTS / "pankratz_v2_n150_first_melt_probe"
    if (probe150 / "first_melt_state.npz").exists():
        cases.append(measure(
            "150³_Pankratz_v2_first_melt_probe", probe150 / "first_melt_state.npz",
            probe150 / "checkpoint" / "source_field_W_m3.npy"))
    for label in ("micro200um_similarity", "exploratory_target1400K"):
        case_dir = RESULTS / f"scaled_n50_{label}"
        if (case_dir / "heated_state.npz").exists():
            cases.append(measure(
                f"50³_{label}", case_dir / "heated_state.npz",
                case_dir / "checkpoint" / "source_field_W_m3.npy"))
    (RESULTS / "LOCAL_SOURCE_REGION_METRICS.json").write_text(
        json.dumps(cases, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    if cases:
        with (RESULTS / "LOCAL_SOURCE_REGION_METRICS.csv").open(
                "w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(cases[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(cases)
    print(json.dumps(cases, ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
