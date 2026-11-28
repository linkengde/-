#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""以相同50³数值载荷测试Ag糊状区宽度0.5、1及2 K。"""
from __future__ import annotations

import json
from pathlib import Path

from compare_50_old_vs_pankratz import run_material_set

HERE = Path(__file__).resolve().parent
OUT = HERE / "results" / "same_load_50"


def main():
    summaries = []
    for width in (0.5, 1.0, 2.0):
        tag = f"Pankratz_v2_{width:g}K"
        summaries.append(run_material_set(tag, True, False, 600.0, width))
    results = {
        "cases": summaries,
        "relative_change_vs_1K": {
            metric: {
                f"{row['Ag_mushy_interval_width_K']:g}K": abs(float(row[metric]) - float(summaries[1][metric]))
                    / max(abs(float(summaries[1][metric])), 1e-30)
                for row in summaries if row is not summaries[1]
            }
            for metric in ("maximum_temperature_K", "Ag_bulk_liquid_fraction_at_heating_end",
                           "Ag_liquid_volume_at_heating_end_m3", "latent_energy_at_heating_end_J")
        },
        "interpretation": "Numerical mushy-width sensitivity under an unchanged local validation pulse; not an arc load.",
    }
    (OUT / "mushy_width_sensitivity.json").write_text(
        json.dumps(results, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(results, ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
