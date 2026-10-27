#!/usr/bin/env python3
"""Check 70:30 3x3 central-source integration and a 1 fs energy step."""

from __future__ import annotations

import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
sys.path.insert(0, str(PROJECT / "phase_field_qsgs"))
from qsgs_enthalpy_solver import QSGSEnthalpyGrid  # noqa: E402

GEOMETRY = HERE / "qsgs_morphology_ag70vol.npz"


def build(tile_xy: int) -> QSGSEnthalpyGrid:
    return QSGSEnthalpyGrid(
        None,
        "A",
        native_geometry_npz=GEOMETRY,
        source_geometry_npz=GEOMETRY,
        tile_xy=tile_xy,
        mesh_refinement=1,
    )


def main() -> None:
    started = time.perf_counter()
    base = build(1)
    base_power = base.actual_source_power_W
    del base
    tile = build(3)
    source_center = list(tile.source_center_xy_A)
    power_difference_relative = (tile.actual_source_power_W - base_power) / base_power
    dt_s = min(tile.stable_timestep_s(), 1.0e-15)
    target_s = 1.0e-15
    while tile.time_s < target_s - 1.0e-30:
        tile.advance(min(dt_s, target_s - tile.time_s), source_on=True)
    total_energy = float(tile.enthalpy_density_J_m3.sum() * tile.cell_volume_m3)
    ledger_error = total_energy - tile.injected_energy_J
    relative_error = ledger_error / tile.injected_energy_J if tile.injected_energy_J else 0.0
    report = {
        "status": "AG70_TILE3_CENTRAL_SOURCE_PREFLIGHT_COMPLETED",
        "geometry_sha256": QSGSEnthalpyGrid._sha256_file(GEOMETRY),
        "tile_xy": 3,
        "grid_shape_xyz": list(tile.temperature_shape),
        "source_center_xy_A": source_center,
        "source_location": "single Gaussian centered in the central repeated QSGS cell; no sources in the other eight cells",
        "source_power_1x_W": base_power,
        "source_power_3x_W": tile.actual_source_power_W,
        "relative_source_power_difference_3x_vs_1x": power_difference_relative,
        "maximum_test_timestep_fs": dt_s * 1.0e15,
        "test_duration_fs": tile.time_s * 1.0e15,
        "injected_energy_J": tile.injected_energy_J,
        "total_enthalpy_J": total_energy,
        "energy_balance_error_J": ledger_error,
        "relative_energy_balance_error": relative_error,
        "peak_temperature_after_test_K": float(tile.state.temperature_K.max()),
        "peak_rss_mib_including_prior_python_peak": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
        "elapsed_wall_time_s": time.perf_counter() - started,
        "grid_cell_count_ratio_to_1x": int(np.prod(tile.temperature_shape) / (100 * 100 * 50)),
        "runtime_estimate_note": "A full 3x3 transient is not launched by this preflight; estimate only after reading the completed 1x runtime and available compute budget.",
    }
    output = HERE / "tile3_ag70_preflight.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    if abs(power_difference_relative) > 1.0e-10:
        raise ArithmeticError("3x3 central source power does not match the 1x source")
    if abs(relative_error) > 1.0e-10:
        raise ArithmeticError("3x3 one-step energy balance failed")


if __name__ == "__main__":
    main()
