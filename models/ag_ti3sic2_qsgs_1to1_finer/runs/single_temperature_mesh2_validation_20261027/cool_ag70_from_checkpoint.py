#!/usr/bin/env python3
"""Continue a saved Ag70 Case A state with the source switched off."""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

RUN_DIR = Path(__file__).resolve().parent
PROJECT_DIR = RUN_DIR.parents[1]
PHASE_DIR = PROJECT_DIR / "phase_field_qsgs"
GEOMETRY = RUN_DIR / "qsgs_morphology_ag70vol.npz"
sys.path.insert(0, str(PHASE_DIR))
from qsgs_enthalpy_solver import QSGSEnthalpyGrid, write_cluster_csv, write_csv  # noqa: E402
from qsgs_liquid_tracking import LiquidClusterTracker  # noqa: E402


def tag(time_ps: float) -> str:
    return f"{time_ps:.6f}".rstrip("0").rstrip(".").replace(".", "p") + "ps"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("--refinement", type=int, choices=(1, 2), default=1)
    parser.add_argument("--cooling-limit-ps", type=float, default=0.5)
    parser.add_argument("--segment-ps", type=float, default=0.02)
    parser.add_argument("--sample-ps", type=float, default=0.002)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--restore-only", action="store_true")
    args = parser.parse_args()
    if min(args.cooling_limit_ps, args.segment_ps, args.sample_ps) <= 0:
        parser.error("cooling limit, segment and sample intervals must be positive")
    if not args.checkpoint.is_file():
        parser.error(f"Missing checkpoint: {args.checkpoint}")
    output_dir = args.output_dir.resolve()
    checkpoints = output_dir / "checkpoints"
    segments = output_dir / "segments"
    checkpoints.mkdir(parents=True, exist_ok=True)
    segments.mkdir(parents=True, exist_ok=True)
    model = QSGSEnthalpyGrid(
        None, "A", native_geometry_npz=GEOMETRY, source_geometry_npz=GEOMETRY,
        mesh_refinement=args.refinement,
    )
    tracker = LiquidClusterTracker(model.spacing_A, model.size_A, liquid_threshold=0.5)
    restart = model.load_checkpoint(args.checkpoint.resolve(), tracker)
    if model._liquid_ag_volume_fraction() <= 0.0 and not np.any(model.state.ag_liquid_fraction > 0.0):
        raise ValueError("The source checkpoint contains no liquid Ag; a cooling branch is not applicable")
    if args.restore_only:
        print(json.dumps({
            "status": "COOLING_CHECKPOINT_RESTORE_VALIDATED",
            "time_ps": model.time_s * 1e12,
            "peak_temperature_K": float(model.state.temperature_K.max()),
            "ag_liquid_volume_fraction": model._liquid_ag_volume_fraction(),
            "maximum_ag_liquid_fraction": float(model.state.ag_liquid_fraction.max()),
            "injected_energy_J": restart["injected_energy_J"],
            "energy_balance_error_J": restart["energy_balance_error_J"],
            "tracker_rows": len(tracker.history),
        }, indent=2))
        return

    started = time.perf_counter()
    initial_time_ps = model.time_s * 1e12
    initial_energy_J = float(model.enthalpy_density_J_m3.sum() * model.cell_volume_m3)
    final_time_limit = initial_time_ps + args.cooling_limit_ps
    all_rows: list[dict[str, object]] = []
    all_solid_time_ps = None
    stop_reason = "cooling_time_limit_reached"
    stage = 0
    while model.time_s * 1e12 < final_time_limit - 1e-12:
        current_ps = model.time_s * 1e12
        end_ps = min(current_ps + args.segment_ps, final_time_limit)
        stage += 1
        rows = model.run(
            end_ps - current_ps,
            source_on=False,
            sample_interval_ps=args.sample_ps,
            track_liquid_connectivity=True,
            liquid_threshold=0.5,
            progress=True,
            tracker=tracker,
            checkpoint_times_ps=(end_ps,),
            checkpoint_dir=checkpoints,
            checkpoint_prefix=f"caseA_ag70_{args.refinement}x_cool",
        )
        all_rows.extend(rows)
        for row in rows:
            if float(row["max_ag_liquid_fraction"]) == 0.0 and all_solid_time_ps is None:
                all_solid_time_ps = float(row["time_ps"])
        total_energy = float(model.enthalpy_density_J_m3.sum() * model.cell_volume_m3)
        error = total_energy - model.injected_energy_J
        rel = error / model.injected_energy_J if model.injected_energy_J else 0.0
        if not np.isfinite(error) or abs(rel) > 1.0e-8:
            raise ArithmeticError(f"Cooling energy ledger failed: relative error {rel:.3e}")
        if not np.all(np.isfinite(model.state.temperature_K)) or not np.all(np.isfinite(model.state.ag_liquid_fraction)):
            raise FloatingPointError("Non-finite cooling state")
        checkpoint = checkpoints / f"caseA_ag70_{args.refinement}x_cool_{tag(end_ps)}.npz"
        if not checkpoint.exists():
            raise FileNotFoundError(f"Missing cooling checkpoint {checkpoint}")
        write_csv(segments / f"cooling_{tag(end_ps)}.csv", rows)
        print(json.dumps({
            "stage": stage,
            "time_ps": end_ps,
            "peak_temperature_K": float(model.state.temperature_K.max()),
            "ag_liquid_volume_fraction": model._liquid_ag_volume_fraction(),
            "max_ag_liquid_fraction": float(model.state.ag_liquid_fraction.max()),
            "connected_cluster_count": int(rows[-1].get("connected_liquid_cluster_count", 0)),
            "largest_cluster_volume_A3": float(rows[-1].get("largest_liquid_cluster_volume_A3", 0.0)),
            "all_solid_sample_time_ps": all_solid_time_ps,
            "energy_balance_error_J": error,
            "relative_energy_balance_error": rel,
            "checkpoint": str(checkpoint),
        }, sort_keys=True), flush=True)
        if all_solid_time_ps is not None:
            stop_reason = "all_ag_liquid_fraction_zero_at_sampled_time"
            break

    final_time_ps = model.time_s * 1e12
    final_checkpoint = checkpoints / f"caseA_ag70_{args.refinement}x_cool_final_{tag(final_time_ps)}.npz"
    model.save_checkpoint(final_checkpoint, tracker)
    if all_rows:
        write_csv(output_dir / "cooling.csv", all_rows)
    if tracker.history:
        write_cluster_csv(output_dir / "cooling.clusters.csv", tracker.history)
    final_energy = float(model.enthalpy_density_J_m3.sum() * model.cell_volume_m3)
    summary = {
        "status": "AG70VOL_CASE_A_SOURCE_OFF_COOLING_COMPLETED",
        "refinement": args.refinement,
        "source_heat_checkpoint": str(args.checkpoint.resolve()),
        "final_checkpoint": str(final_checkpoint.relative_to(output_dir)),
        "time_heat_off_ps": initial_time_ps,
        "time_end_ps": final_time_ps,
        "cooling_elapsed_ps": final_time_ps - initial_time_ps,
        "first_all_solid_sample_time_ps": all_solid_time_ps,
        "solidification_duration_from_heat_off_ps": None if all_solid_time_ps is None else all_solid_time_ps - initial_time_ps,
        "peak_temperature_final_K": float(model.state.temperature_K.max()),
        "ag_liquid_volume_fraction_final": model._liquid_ag_volume_fraction(),
        "max_ag_liquid_fraction_final": float(model.state.ag_liquid_fraction.max()),
        "initial_total_enthalpy_J": initial_energy_J,
        "final_total_enthalpy_J": final_energy,
        "final_energy_balance_error_J": final_energy - model.injected_energy_J,
        "relative_energy_balance_error": (final_energy - model.injected_energy_J) / model.injected_energy_J if model.injected_energy_J else 0.0,
        "elapsed_wall_time_s": time.perf_counter() - started,
        "stop_reason": stop_reason,
        "interpretation": "Insulated Case A source-off run: internal energy redistribution and latent-to-sensible conversion only; not heat loss to the environment.",
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
