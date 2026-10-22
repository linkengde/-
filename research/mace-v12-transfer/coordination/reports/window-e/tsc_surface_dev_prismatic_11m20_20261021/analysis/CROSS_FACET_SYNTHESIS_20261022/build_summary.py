#!/usr/bin/env python3
"""Build a descriptive cross-facet summary from two archived run0 summaries."""
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
SOURCES = {
    "(10-10) 96 atom": ROOT / "tsc_surface_dev_prismatic_20261020/analysis/DEV_TSCBOP_RUN0_PAIR_20261021/SUMMARY.json",
    "(11-20) 120 atom": ROOT / "tsc_surface_dev_prismatic_11m20_20261021/analysis/TSCBOP_RUN0_PAIR_20261022/SUMMARY.json",
}

rows = []
for facet, path in SOURCES.items():
    data = json.loads(path.read_text())
    for geom_name in ("REF", "RESP_TIZ020"):
        geom = data["structures"][geom_name]
        for group in ("all", "by_element"):
            values = geom[group]
            if group == "all":
                rows.append({
                    "facet": facet, "quantity": "absolute_force_error", "geometry": geom_name,
                    "group": "all", "n": values["atom_count"],
                    "component_rmse_eV_A": values["component_rmse_eV_A"],
                    "vector_rmse_eV_A": values["vector_rmse_eV_A"],
                    "p95_vector_error_eV_A": values["p95_vector_error_eV_A"],
                    "max_vector_error_eV_A": values["max_vector_error_eV_A"],
                })
            else:
                for element, stats in values.items():
                    rows.append({
                        "facet": facet, "quantity": "absolute_force_error", "geometry": geom_name,
                        "group": element, "n": stats["atom_count"],
                        "component_rmse_eV_A": stats["component_rmse_eV_A"],
                        "vector_rmse_eV_A": stats["vector_rmse_eV_A"],
                        "p95_vector_error_eV_A": stats["p95_vector_error_eV_A"],
                        "max_vector_error_eV_A": stats["max_vector_error_eV_A"],
                    })
    response = data["response_error"]
    local_key = "within_3p3A_displaced" if "within_3p3A_displaced" in response else "within_3p3A_of_displaced_Ti"
    for region, stats in (("all", response["all"]), ("local_3p3A", response[local_key])):
        rows.append({
            "facet": facet, "quantity": "differential_response_error", "geometry": "RESP_minus_REF",
            "group": region, "n": stats["atom_count"],
            "component_rmse_eV_A": stats["component_rmse_eV_A"],
            "vector_rmse_eV_A": stats["vector_rmse_eV_A"],
            "p95_vector_error_eV_A": stats["p95_vector_error_eV_A"],
            "max_vector_error_eV_A": stats["max_vector_error_eV_A"],
        })

fields = ["facet", "quantity", "geometry", "group", "n", "component_rmse_eV_A",
          "vector_rmse_eV_A", "p95_vector_error_eV_A", "max_vector_error_eV_A"]
with (OUT / "CROSS_FACET_FORCE_SUMMARY.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

print(f"wrote {len(rows)} metric rows to {OUT / 'CROSS_FACET_FORCE_SUMMARY.csv'}")
