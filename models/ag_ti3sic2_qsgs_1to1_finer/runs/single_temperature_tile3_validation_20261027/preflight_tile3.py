"""Validate a 3x3 Case A tile and run a short central-cell comparison."""

from __future__ import annotations

import json
import resource
import sys
import time
from pathlib import Path

import numpy as np


RUN_DIR = Path(__file__).resolve().parent
PROJECT_DIR = RUN_DIR.parents[1]
PHASE_DIR = PROJECT_DIR / "phase_field_qsgs"
sys.path.insert(0, str(PHASE_DIR))
from qsgs_enthalpy_solver import QSGSEnthalpyGrid, ThermalState  # noqa: E402


def main() -> None:
    geometry = PROJECT_DIR / "qsgs_morphology.npz"
    with np.load(geometry, allow_pickle=False) as archive:
        source_ag = ~np.asarray(archive["skeleton"], dtype=bool)
        size_A = np.asarray(archive["size_A"], dtype=float)
    shape = np.asarray(source_ag.shape)

    wall0 = time.perf_counter()
    base = QSGSEnthalpyGrid(None, "A", native_geometry_npz=geometry, tile_xy=1)
    base_init_s = time.perf_counter() - wall0
    wall1 = time.perf_counter()
    tiled = QSGSEnthalpyGrid(None, "A", native_geometry_npz=geometry, tile_xy=3)
    tiled_init_s = time.perf_counter() - wall1

    expected_ag = np.tile(base.ag_fraction, (3, 3, 1))
    geometry_exact = bool(np.array_equal(tiled.ag_fraction, expected_ag))
    if not geometry_exact:
        raise RuntimeError("3x3 geometry is not an exact lateral repetition of the original QSGS mask")

    middle = (slice(int(shape[0]), int(2 * shape[0])), slice(int(shape[1]), int(2 * shape[1])), slice(None))
    base_source = base.source_cell_power_W
    central_source = tiled.source_cell_power_W[middle]
    source_central_abs_diff_W = float(np.max(np.abs(central_source - base_source)))
    source_central_relative_l1 = float(
        np.sum(np.abs(central_source - base_source)) / max(np.sum(np.abs(base_source)), 1.0e-300)
    )
    central_source_power_W = float(central_source.sum())
    total_source_power_W = float(tiled.source_cell_power_W.sum())
    outer_source_power_W = total_source_power_W - central_source_power_W
    source_center_xy_A = tiled.source_center_xy_A
    expected_center_xy_A = (1.5 * float(size_A[0]), 1.5 * float(size_A[1]))

    base_stable = base.stable_timestep_s()
    tiled_stable = tiled.stable_timestep_s()
    if not np.isclose(base_stable, tiled_stable, rtol=1.0e-13, atol=0.0):
        raise RuntimeError("Local explicit stability limit changed under exact tiling")

    # A smooth nonuniform field checks that the periodic x/y, insulated-z
    # diffusion operator still conserves total energy on the expanded domain.
    nx, ny, nz = tiled.temperature_shape
    x = (np.arange(nx) + 0.5) / nx
    y = (np.arange(ny) + 0.5) / ny
    temp = 300.0 + 25.0 * np.sin(2.0 * np.pi * x[:, None, None]) * np.cos(
        2.0 * np.pi * y[None, :, None]
    )
    temp = np.broadcast_to(temp, tiled.temperature_shape).copy()
    diffusion_state = ThermalState(temp, np.zeros(tiled.temperature_shape, dtype=float))
    conductivity = tiled.conductivity_W_m_K(diffusion_state)
    diffusion = tiled.diffusion_power_density_W_m3(temp, conductivity)
    diffusion_integral_W = float(diffusion.sum() * tiled.cell_volume_m3)
    diffusion_scale_W = float(np.abs(diffusion).sum() * tiled.cell_volume_m3)
    diffusion_relative_residual = abs(diffusion_integral_W) / max(diffusion_scale_W, 1.0e-300)

    # One short heating interval on both domains measures central-cell aliasing
    # and the actual cost before committing to a multi-hour candidate run.
    sample_ps = 0.005
    wall2 = time.perf_counter()
    base_rows = base.run(sample_ps, source_on=True, sample_interval_ps=sample_ps)
    base_run_s = time.perf_counter() - wall2
    wall3 = time.perf_counter()
    tiled_rows = tiled.run(sample_ps, source_on=True, sample_interval_ps=sample_ps)
    tiled_run_s = time.perf_counter() - wall3

    base_temperature = base.state.temperature_K
    central_temperature = tiled.state.temperature_K[middle]
    central_ag = base.ag_fraction > 0.0
    temperature_delta = central_temperature - base_temperature
    temperature_scale = float(np.sum(central_ag * np.abs(base_temperature - 300.0)))
    central_l1_over_rise = float(
        np.sum(central_ag * np.abs(temperature_delta)) / max(temperature_scale, 1.0e-300)
    )
    central_max_delta_K = float(np.max(np.abs(temperature_delta[central_ag])))
    central_peak_delta_K = float(
        np.max(central_temperature[central_ag]) - np.max(base_temperature[central_ag])
    )
    base_energy = float(base.enthalpy_density_J_m3.sum() * base.cell_volume_m3)
    tiled_energy = float(tiled.enthalpy_density_J_m3.sum() * tiled.cell_volume_m3)
    base_energy_error = base_energy - base.injected_energy_J
    tiled_energy_error = tiled_energy - tiled.injected_energy_J

    report = {
        "status": "3X3_TILE_PREFLIGHT_PASSED",
        "geometry_input": "qsgs_morphology.npz",
        "geometry_sha256": QSGSEnthalpyGrid._sha256_file(geometry),
        "base_shape_xyz": base.temperature_shape,
        "tile3_shape_xyz": tiled.temperature_shape,
        "tile3_geometry_exact_lateral_repeat": geometry_exact,
        "base_spacing_A": base.spacing_A.tolist(),
        "tile3_spacing_A": tiled.spacing_A.tolist(),
        "boundary_case": "A; x/y periodic, z insulated; scaffold has zero thermal capacity and conductivity.",
        "source_integration_method": tiled.source_integration_method,
        "source_center_xy_A": source_center_xy_A,
        "expected_central_cell_source_center_xy_A": expected_center_xy_A,
        "base_source_power_W": base.actual_source_power_W,
        "tile3_total_source_power_W": tiled.actual_source_power_W,
        "tile3_source_power_in_central_cell_W": central_source_power_W,
        "tile3_source_power_outside_central_cell_W": outer_source_power_W,
        "tile3_source_power_outside_central_cell_fraction": outer_source_power_W / total_source_power_W,
        "central_source_map_max_abs_difference_W_per_cell": source_central_abs_diff_W,
        "central_source_map_relative_L1_difference": source_central_relative_l1,
        "base_initial_stable_dt_fs": base_stable * 1.0e15,
        "tile3_initial_stable_dt_fs": tiled_stable * 1.0e15,
        "one_step_operator_integrated_diffusion_power_W": diffusion_integral_W,
        "one_step_operator_relative_energy_residual": diffusion_relative_residual,
        "short_heating_interval_ps": sample_ps,
        "base_initialization_wall_s": base_init_s,
        "tile3_initialization_wall_s": tiled_init_s,
        "base_short_heating_wall_s": base_run_s,
        "tile3_short_heating_wall_s": tiled_run_s,
        "tile3_short_heating_relative_wall_cost": tiled_run_s / max(base_run_s, 1.0e-12),
        "base_short_run_energy_error_J": base_energy_error,
        "tile3_short_run_energy_error_J": tiled_energy_error,
        "base_short_run_injected_energy_J": base.injected_energy_J,
        "tile3_short_run_injected_energy_J": tiled.injected_energy_J,
        "central_cell_peak_temperature_base_K": float(np.max(base_temperature[central_ag])),
        "central_cell_peak_temperature_tile3_K": float(np.max(central_temperature[central_ag])),
        "central_cell_peak_temperature_difference_K": central_peak_delta_K,
        "central_cell_max_abs_temperature_difference_K": central_max_delta_K,
        "central_cell_Ag_weighted_temperature_L1_over_base_rise": central_l1_over_rise,
        "central_cell_Ag_liquid_volume_fraction_base": base._liquid_ag_volume_fraction(),
        "central_cell_Ag_liquid_volume_fraction_tile3": float(
            np.sum(base.ag_fraction * tiled.state.ag_liquid_fraction[middle]) / base.ag_fraction.sum()
        ),
        "whole_tile3_Ag_liquid_volume_fraction": tiled._liquid_ag_volume_fraction(),
        "max_rss_kib": int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
        "sample_rows_base": len(base_rows),
        "sample_rows_tile3": len(tiled_rows),
        "interpretation": "Short-run mapping, source, diffusion conservation, and central-cell comparison only; this does not test hotspot/pool finite-size convergence after substantial lateral heat spread.",
    }
    (RUN_DIR / "tile3_preflight.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
