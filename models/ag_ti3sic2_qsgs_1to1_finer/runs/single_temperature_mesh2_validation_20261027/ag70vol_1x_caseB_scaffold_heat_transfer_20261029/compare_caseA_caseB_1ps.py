#!/usr/bin/env python3
"""Compare existing Case A and new Case B Ag70 1x states at 1 ps."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np


OUT_DIR = Path(__file__).resolve().parent
RUN_DIR = OUT_DIR.parent
PROJECT_DIR = OUT_DIR.parents[2]
PHASE_DIR = PROJECT_DIR / "phase_field_qsgs"
GEOMETRY = RUN_DIR / "qsgs_morphology_ag70vol.npz"
CASE_A_CHECKPOINT = RUN_DIR / "ag70vol_resume_1x/checkpoints/caseA_ag70_1x_final_1ps.npz"
CASE_B_CHECKPOINT = OUT_DIR / "checkpoints/caseB_ag70_1x_heat_final_1ps.npz"
sys.path.insert(0, str(PHASE_DIR))

from qsgs_enthalpy_solver import QSGSEnthalpyGrid  # noqa: E402
from qsgs_liquid_tracking import LiquidClusterTracker  # noqa: E402
from run_ag70_1x_caseB import energy_breakdown  # noqa: E402


def checkpoint_row(boundary: str, path: Path) -> dict[str, object]:
    model = QSGSEnthalpyGrid(
        None,
        boundary,
        native_geometry_npz=GEOMETRY,
        source_geometry_npz=GEOMETRY,
        mesh_refinement=1,
    )
    tracker = LiquidClusterTracker(model.spacing_A, model.size_A, liquid_threshold=0.5)
    model.load_checkpoint(path, tracker)
    ag_mask = model.ag_fraction > 0.0
    sc_mask = model.sc_fraction > 0.0
    ag_mean = float(np.sum(model.ag_fraction * model.state.temperature_K) / model.ag_fraction.sum())
    row: dict[str, object] = {
        "model": f"Case {boundary}",
        "boundary_case": boundary,
        "time_ps": model.time_s * 1.0e12,
        "geometry_sha256": QSGSEnthalpyGrid._sha256_file(GEOMETRY),
        "ag_peak_temperature_K": float(model.state.temperature_K[ag_mask].max()),
        "ag_mean_temperature_K": ag_mean,
        "ti3sic2_peak_temperature_K": float(model.state.temperature_K[sc_mask].max()) if boundary == "B" else None,
        "ti3sic2_mean_temperature_K": float(np.sum(model.sc_fraction * model.state.temperature_K) / model.sc_fraction.sum()) if boundary == "B" else None,
        "ag_liquid_volume_fraction": model._liquid_ag_volume_fraction(),
        "ag_latent_enthalpy_J": model._ag_latent_energy_J(),
        "system_total_enthalpy_J": float(model.enthalpy_density_J_m3.sum() * model.cell_volume_m3),
        "cumulative_source_energy_J": float(model.injected_energy_J),
        "energy_balance_residual_J": float(model.enthalpy_density_J_m3.sum() * model.cell_volume_m3 - model.injected_energy_J),
        "checkpoint_sha256": QSGSEnthalpyGrid._sha256_file(path),
    }
    if boundary == "B":
        row["ag_sensible_enthalpy_J"] = energy_breakdown(model)["ag_sensible_enthalpy_J"]
        row["ti3sic2_stored_enthalpy_J"] = energy_breakdown(model)["ti3sic2_stored_enthalpy_J"]
    else:
        row["ag_sensible_enthalpy_J"] = row["system_total_enthalpy_J"] - row["ag_latent_enthalpy_J"]
        row["ti3sic2_stored_enthalpy_J"] = 0.0
    return row


def main() -> None:
    if not CASE_A_CHECKPOINT.is_file():
        raise FileNotFoundError(CASE_A_CHECKPOINT)
    if not CASE_B_CHECKPOINT.is_file():
        raise FileNotFoundError(CASE_B_CHECKPOINT)
    rows = [checkpoint_row("A", CASE_A_CHECKPOINT), checkpoint_row("B", CASE_B_CHECKPOINT)]
    if not np.isclose(float(rows[0]["time_ps"]), float(rows[1]["time_ps"]), rtol=0.0, atol=1.0e-10):
        raise ValueError("Case A and Case B checkpoints are not at the same time")
    peak_delta = float(rows[1]["ag_peak_temperature_K"]) - float(rows[0]["ag_peak_temperature_K"])
    mean_delta = float(rows[1]["ag_mean_temperature_K"]) - float(rows[0]["ag_mean_temperature_K"])
    result = {
        "status": "MATCHED_TIME_CASE_A_CASE_B_COMPARISON",
        "time_ps": float(rows[0]["time_ps"]),
        "geometry_hashes_match": rows[0]["geometry_sha256"] == rows[1]["geometry_sha256"],
        "ag_peak_temperature_delta_B_minus_A_K": peak_delta,
        "ag_peak_temperature_reduction_percent_vs_A": -100.0 * peak_delta / float(rows[0]["ag_peak_temperature_K"]),
        "ag_mean_temperature_delta_B_minus_A_K": mean_delta,
        "both_models_have_zero_liquid": all(float(row["ag_liquid_volume_fraction"]) == 0.0 for row in rows),
        "models": rows,
        "interpretation": "At matched 1 ps, Case B stores heat in Ti3SiC2 and conducts it through both phases; this comparison changes the thermal boundary case, not the QSGS geometry or source parameters.",
        "limitations": [
            "candidate properties remain incompletely source-verified",
            "Case B assumes ideal interface contact and insulated external boundaries",
            "this is a thermal response comparison, not melt-flow or mechanical weld simulation",
        ],
    }
    (OUT_DIR / "caseA_caseB_matched_1ps.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    with (OUT_DIR / "caseA_caseB_matched_1ps.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
