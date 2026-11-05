#!/usr/bin/env python3
"""Cross-check P1H pressure loading and matched-time surface response.

This audit reads retained P1H fixed-grid binary snapshots and run CSV files.
The binary handoff is a projection of the octree fields, so the reconstructed
surface map is a diagnostic visualization; the solver's conservative CSVs
remain authoritative for integrated material volumes.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import struct
from pathlib import Path

import numpy as np

HEADER = struct.Struct("<8s6I5d")
SOURCE_X_NM = 9.58832
SOURCE_Y_NM = 9.65943950771472
SCENARIOS = {
    "zero_pressure": "production_pair/marangoni_on_zero_pressure",
    "gaussian_1p8MPa": "production_pair/marangoni_on_gaussian_1p8MPa",
    "gaussian_16p7MPa": "pressure_sensitivity_pi0p1/marangoni_on_gaussian_16p7MPa",
    "gaussian_167MPa": "pressure_sensitivity_pi1/marangoni_on_gaussian_167MPa",
}


def read_frame(path: Path, expected_magic: bytes, expected_fields: int):
    with gzip.open(path, "rb") as stream:
        raw = stream.read()
    if len(raw) < HEADER.size:
        raise ValueError(f"short frame: {path}")
    (magic, version, nx, ny, nz, nfields, case_flag, time_ps, dx_m,
     x0_m, y0_m, z0_m) = HEADER.unpack_from(raw)
    if magic != expected_magic or version != 1 or nfields != expected_fields:
        raise ValueError(f"unexpected P1H header in {path}")
    expected_bytes = HEADER.size + 4*nfields*nx*ny*nz
    if len(raw) != expected_bytes:
        raise ValueError(f"incomplete frame: {path}")
    fields = np.frombuffer(raw, dtype="<f4", offset=HEADER.size)
    fields = fields.reshape(nfields, nx, ny, nz)
    if not np.isfinite(fields).all():
        raise ValueError(f"non-finite field value in {path}")
    metadata = dict(shape=(nx, ny, nz), time_ps=time_ps, dx_m=dx_m,
                    origin_m=(x0_m, y0_m, z0_m), case_flag=case_flag)
    return fields, metadata


def read_csv(path: Path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    order = np.argsort(values)
    sorted_values = values[order]
    cumulative = np.cumsum(weights[order])
    if not cumulative.size or cumulative[-1] <= 0.:
        return float("nan")
    return float(sorted_values[min(np.searchsorted(cumulative, q*cumulative[-1]),
                                   len(sorted_values) - 1)])


def nearest_row(rows, time_ps: float):
    row = min(rows, key=lambda item: abs(float(item["time_ps"]) - time_ps))
    if abs(float(row["time_ps"]) - time_ps) > 1.e-5:
        raise ValueError(f"no matching time near {time_ps:.9f} ps")
    return row


def projected_force_z(vof: np.ndarray, p_effective: np.ndarray, dx_m: float):
    """Mirror the solver's z-face pressure-force sum on the projected grid."""
    current_interface = (vof > 1.e-8) & (vof < 1. - 1.e-8)
    lower_vof = np.roll(vof, 1, axis=2)
    lower_interface = np.roll(current_interface, 1, axis=2)
    delta_f = vof - lower_vof
    use_face = (delta_f != 0.) & (current_interface | lower_interface)
    use_face[:, :, 0] = False  # no periodic wrap in z
    lower_pressure = np.roll(p_effective, 1, axis=2)
    p_face = np.where(
        current_interface & lower_interface,
        0.5*(p_effective + lower_pressure),
        np.where(current_interface, p_effective, lower_pressure),
    )
    return float(np.sum(np.where(use_face, p_face*delta_f/dx_m, 0.)) * dx_m**3)


def surface_height_change(initial_vof: np.ndarray, current_vof: np.ndarray,
                          metadata: dict):
    nx, ny, nz = current_vof.shape
    dx_nm = metadata["dx_m"]*1.e9
    z0_nm = metadata["origin_m"][2]*1.e9
    z_faces_nm = z0_nm + (np.arange(nz) + 1.)*dx_nm
    initial_height = np.max(np.where(initial_vof > 0.5,
                                     z_faces_nm[None, None, :], -np.inf), axis=2)
    current_height = np.max(np.where(current_vof > 0.5,
                                     z_faces_nm[None, None, :], -np.inf), axis=2)
    change = initial_height - current_height
    x = metadata["origin_m"][0]*1.e9 + (np.arange(nx) + 0.5)*dx_nm
    y = metadata["origin_m"][1]*1.e9 + (np.arange(ny) + 0.5)*dx_nm
    xx, yy = np.meshgrid(x, y, indexing="ij")
    radius = np.hypot(xx - SOURCE_X_NM, yy - SOURCE_Y_NM)
    valid = np.isfinite(initial_height) & np.isfinite(current_height)
    return change, valid, radius


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True,
                        help="root containing the retained four-scenario P1H runs")
    parser.add_argument("--comparison-csv", type=Path, required=True,
                        help="committed arc_pressure_exploration_comparison.csv")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    comparison = {r["scenario"]: r for r in read_csv(args.comparison_csv)}
    rows_out = []
    baseline_change = None
    for scenario, relative_dir in SCENARIOS.items():
        run_dir = args.run_root / relative_dir
        main_path = run_dir / "p1h_marangoni_on_regular_t000.504856_ps.bin.gz"
        pressure_path = run_dir / "p1h_pressure_marangoni_on_regular_t000.504856_ps.bin.gz"
        initial_path = run_dir / "p1h_marangoni_on_initial_t000.000000_ps.bin.gz"
        main_fields, meta = read_frame(main_path, b"P1HSNAP1", 20)
        pressure_fields, pressure_meta = read_frame(pressure_path, b"P1HPRS01", 4)
        initial_fields, initial_meta = read_frame(initial_path, b"P1HSNAP1", 20)
        compare_keys = ("shape", "time_ps", "dx_m", "origin_m")
        if (any(meta[key] != pressure_meta[key] for key in compare_keys) or
            any(meta[key] != initial_meta[key] for key in ("shape", "dx_m", "origin_m"))):
            raise ValueError(f"snapshot metadata mismatch in {scenario}")

        force = projected_force_z(main_fields[0], pressure_fields[2], meta["dx_m"])
        pressure_history = read_csv(run_dir / "p1h_arc_pressure_history.csv")
        history_at_time = nearest_row(pressure_history, meta["time_ps"])
        logged_force = float(history_at_time["integrated_pressure_force_z_N"])
        force_relative_error = (abs(force - logged_force)/abs(logged_force)
                                if logged_force else abs(force - logged_force))
        if force_relative_error > 1.e-6:
            raise ValueError(f"projected and logged pressure forces disagree in {scenario}: "
                             f"{force_relative_error:.3g}")
        configured_pressure = float(history_at_time["configured_peak_Pa"])
        equivalent_area_nm2 = (abs(force)/configured_pressure*1.e18
                               if configured_pressure > 0. else 0.)

        history = read_csv(run_dir / "p1h_history.csv")
        thermal_row = nearest_row(history, meta["time_ps"])
        surface_change, valid, radius = surface_height_change(
            initial_fields[0], main_fields[0], meta)
        center = valid & (radius <= 5.)
        ring = valid & (radius > 5.) & (radius <= 8.)
        if not center.any() or not ring.any():
            raise ValueError("empty source or annular region in surface map")
        if scenario == "zero_pressure":
            baseline_change = surface_change.copy()
        center_response = float(np.mean(surface_change[center]))
        ring_response = float(np.mean(surface_change[ring]))
        delta_center_vs_zero = (0. if baseline_change is None else
                                float(np.mean((surface_change - baseline_change)[center])))
        delta_ring_vs_zero = (0. if baseline_change is None else
                              float(np.mean((surface_change - baseline_change)[ring])))

        # Track whether the conservatively advected Ag occupancy actually
        # migrates radially, and quantify its relation to the fixed skeleton.
        dx_nm = meta["dx_m"]*1.e9
        cell_volume_nm3 = dx_nm**3
        ag_now = np.maximum(main_fields[2], 0.)
        ag_initial = np.maximum(initial_fields[2], 0.)
        liquid_ag = np.maximum(main_fields[4], 0.)
        skeleton = np.maximum(main_fields[3], 0.)
        center_xy = radius <= 5.
        annulus_xy = (radius > 5.) & (radius <= 8.)
        outer_xy = radius > 8.
        def regional_volume(field, region):
            return float(np.sum(field[region, :])*cell_volume_nm3)
        ag_total_volume = float(np.sum(ag_now)*cell_volume_nm3)
        ag_initial_volume = float(np.sum(ag_initial)*cell_volume_nm3)
        ag_center_change = regional_volume(ag_now - ag_initial, center_xy)
        ag_annulus_change = regional_volume(ag_now - ag_initial, annulus_xy)
        ag_outer_change = regional_volume(ag_now - ag_initial, outer_xy)
        liquid_total_volume = float(np.sum(liquid_ag)*cell_volume_nm3)
        moved_equivalent_volume = float(0.5*np.sum(np.abs(ag_now - ag_initial))*cell_volume_nm3)
        skeleton_volume = float(np.sum(skeleton)*cell_volume_nm3)

        # Compare the discrete pressure and capillary potentials where the
        # pressure mask is active and pressure is at least 10% of its maximum.
        effective_pressure = pressure_fields[2]
        capillary_potential = pressure_fields[3] - effective_pressure
        local_pressure_max = float(np.max(effective_pressure))
        stress_mask = ((pressure_fields[1] > 0.) &
                       (effective_pressure >= 0.1*local_pressure_max) &
                       (local_pressure_max > 0.))
        stress_pressure = effective_pressure[stress_mask]
        stress_capillary = np.abs(capillary_potential[stress_mask])
        stress_curvature = np.abs(main_fields[11][stress_mask])
        if stress_pressure.size:
            pressure_capillary_ratio = stress_pressure/(stress_capillary + 1.e-30)
            pressure_median_MPa = float(np.median(stress_pressure)/1.e6)
            capillary_median_MPa = float(np.median(stress_capillary)/1.e6)
            capillary_p90_MPa = float(np.quantile(stress_capillary,0.9)/1.e6)
            capillary_max_MPa = float(np.max(stress_capillary)/1.e6)
            pressure_capillary_ratio_median = float(np.median(pressure_capillary_ratio))
            pressure_capillary_ratio_p10 = float(np.quantile(pressure_capillary_ratio,0.1))
            curvature_median_per_nm = float(np.median(stress_curvature)/1.e9)
            sigma_weighted_median_N_m = float(np.median(main_fields[13][stress_mask]))
        else:
            pressure_median_MPa = capillary_median_MPa = capillary_p90_MPa = float("nan")
            capillary_max_MPa = pressure_capillary_ratio_median = float("nan")
            pressure_capillary_ratio_p10 = curvature_median_per_nm = float("nan")
            sigma_weighted_median_N_m = float("nan")

        liquid_weights = liquid_ag.ravel()
        liquid_drag = main_fields[18].ravel()
        liquid_speed = np.sqrt(main_fields[7]**2 + main_fields[8]**2 + main_fields[9]**2).ravel()
        liquid_total_weight = float(liquid_weights.sum())
        if liquid_total_weight > 0.:
            drag_q50 = weighted_quantile(liquid_drag, liquid_weights, 0.50)
            drag_q90 = weighted_quantile(liquid_drag, liquid_weights, 0.90)
            drag_q99 = weighted_quantile(liquid_drag, liquid_weights, 0.99)
            liquid_drag_frac_1e9 = float(liquid_weights[liquid_drag >= 1.e9].sum()/liquid_total_weight)
            liquid_drag_frac_1e10 = float(liquid_weights[liquid_drag >= 1.e10].sum()/liquid_total_weight)
            liquid_drag_frac_1e11 = float(liquid_weights[liquid_drag >= 1.e11].sum()/liquid_total_weight)
            liquid_mean_speed = float(np.dot(liquid_speed, liquid_weights)/liquid_total_weight)
            sc_cells = skeleton > 1.e-6
            near_sc = sc_cells.copy()
            for axis in (0, 1):
                near_sc |= np.roll(sc_cells, 1, axis=axis) | np.roll(sc_cells, -1, axis=axis)
            for shift in (-1, 1):
                shifted = np.zeros_like(sc_cells)
                if shift < 0:
                    shifted[:, :, :-1] = sc_cells[:, :, 1:]
                else:
                    shifted[:, :, 1:] = sc_cells[:, :, :-1]
                near_sc |= shifted
            near_mask = near_sc.ravel() & (liquid_weights > 0.)
            far_mask = (~near_sc).ravel() & (liquid_weights > 0.)
            near_weight = float(liquid_weights[near_mask].sum())
            far_weight = float(liquid_weights[far_mask].sum())
            near_liquid_speed = (float(np.dot(liquid_speed[near_mask], liquid_weights[near_mask])/near_weight)
                                 if near_weight else float("nan"))
            far_liquid_speed = (float(np.dot(liquid_speed[far_mask], liquid_weights[far_mask])/far_weight)
                                if far_weight else float("nan"))
            near_liquid_drag = (float(np.dot(liquid_drag[near_mask], liquid_weights[near_mask])/near_weight)
                                if near_weight else float("nan"))
            far_liquid_drag = (float(np.dot(liquid_drag[far_mask], liquid_weights[far_mask])/far_weight)
                               if far_weight else float("nan"))
        else:
            drag_q50 = drag_q90 = drag_q99 = float("nan")
            liquid_drag_frac_1e9 = liquid_drag_frac_1e10 = liquid_drag_frac_1e11 = float("nan")
            liquid_mean_speed = near_liquid_speed = far_liquid_speed = float("nan")
            near_liquid_drag = far_liquid_drag = float("nan")
            near_weight = far_weight = 0.

        run_pressure = np.array([float(r["max_effective_interface_Pa"])
                                 for r in pressure_history])
        run_force = np.array([float(r["integrated_pressure_force_z_N"])
                              for r in pressure_history])
        run_time_s = np.array([float(r["time_ps"])*1.e-12 for r in pressure_history])
        peak_force_index = int(np.argmax(np.abs(run_force)))
        peak_effective_pressure = float(np.max(run_pressure))
        impulse = float(np.trapezoid(run_force, run_time_s))
        active = run_pressure > max(configured_pressure, 1.)*1.e-6
        active_force = run_force[active]
        active_load_all_down = (bool(np.all(active_force < 0.))
                                if active_force.size else "")
        if active_force.size and not active_load_all_down:
            raise ValueError(f"upward net pressure force found while load is active in {scenario}")

        summary = comparison[scenario]
        rows_out.append({
            "scenario": scenario,
            "configured_peak_pressure_Pa": configured_pressure,
            "common_snapshot_time_ps": meta["time_ps"],
            "snapshot_max_raw_pressure_Pa": float(np.max(pressure_fields[0])),
            "snapshot_max_interface_pressure_Pa": float(np.max(pressure_fields[2])),
            "snapshot_pressure_support_cells": int(np.count_nonzero(pressure_fields[2] > 0.)),
            "snapshot_force_z_reconstructed_N": force,
            "snapshot_force_z_logged_N": logged_force,
            "snapshot_force_relative_difference": force_relative_error,
            "snapshot_equivalent_projected_area_nm2": equivalent_area_nm2,
            "source_radius_mean_surface_change_nm": center_response,
            "annulus_5_to_8nm_mean_surface_change_nm": ring_response,
            "source_radius_response_vs_zero_nm": delta_center_vs_zero,
            "annulus_response_vs_zero_nm": delta_ring_vs_zero,
            "source_radius_max_downward_change_nm": float(np.max(surface_change[center])),
            "source_radius_downward_columns": int(np.count_nonzero(surface_change[center] > 0.)),
            "max_velocity_at_common_time_m_s": float(thermal_row["max_velocity_m_s"]),
            "liquid_Ag_volume_at_common_time_nm3": float(thermal_row["liquid_Ag_volume_nm3"]),
            "projected_Ag_total_volume_nm3": ag_total_volume,
            "projected_Ag_volume_drift_pct": 100.*(ag_total_volume/ag_initial_volume - 1.) if ag_initial_volume else float("nan"),
            "Ag_marker_moved_equivalent_volume_nm3": moved_equivalent_volume,
            "Ag_marker_moved_equivalent_fraction": moved_equivalent_volume/ag_initial_volume if ag_initial_volume else float("nan"),
            "Ag_center_volume_change_nm3": ag_center_change,
            "Ag_annulus_volume_change_nm3": ag_annulus_change,
            "Ag_outer_volume_change_nm3": ag_outer_change,
            "projected_liquid_Ag_volume_nm3": liquid_total_volume,
            "projected_skeleton_volume_nm3": skeleton_volume,
            "pressure_load_cell_count_above_10pct_local_peak": int(stress_pressure.size),
            "load_cell_median_effective_pressure_MPa": pressure_median_MPa,
            "load_cell_median_abs_capillary_potential_MPa": capillary_median_MPa,
            "load_cell_p90_abs_capillary_potential_MPa": capillary_p90_MPa,
            "load_cell_max_abs_capillary_potential_MPa": capillary_max_MPa,
            "load_cell_median_pressure_to_capillary_ratio": pressure_capillary_ratio_median,
            "load_cell_p10_pressure_to_capillary_ratio": pressure_capillary_ratio_p10,
            "load_cell_median_abs_curvature_per_nm": curvature_median_per_nm,
            "load_cell_median_weighted_surface_tension_N_m": sigma_weighted_median_N_m,
            "liquid_Ag_weighted_drag_q50_s-1": drag_q50,
            "liquid_Ag_weighted_drag_q90_s-1": drag_q90,
            "liquid_Ag_weighted_drag_q99_s-1": drag_q99,
            "liquid_Ag_fraction_drag_ge_1e9": liquid_drag_frac_1e9,
            "liquid_Ag_fraction_drag_ge_1e10": liquid_drag_frac_1e10,
            "liquid_Ag_fraction_drag_ge_1e11": liquid_drag_frac_1e11,
            "liquid_Ag_weighted_mean_speed_m_s": liquid_mean_speed,
            "liquid_Ag_weighted_mean_speed_near_skeleton_m_s": near_liquid_speed,
            "liquid_Ag_weighted_mean_speed_away_skeleton_m_s": far_liquid_speed,
            "liquid_Ag_weighted_mean_drag_near_skeleton_s-1": near_liquid_drag,
            "liquid_Ag_weighted_mean_drag_away_skeleton_s-1": far_liquid_drag,
            "peak_effective_pressure_over_run_Pa": peak_effective_pressure,
            "peak_force_z_over_run_N": float(run_force[peak_force_index]),
            "peak_force_time_ps": float(pressure_history[peak_force_index]["time_ps"]),
            "peak_force_equivalent_area_nm2": abs(float(run_force[peak_force_index]))/configured_pressure*1.e18 if configured_pressure else 0.,
            "integrated_pressure_impulse_N_s": impulse,
            "active_pressure_history_samples": int(active.sum()),
            "all_active_pressure_forces_downward": active_load_all_down,
            "heat_off_center_surface_metric_nm": float(summary["center_surface_metric_at_heat_off_nm"]),
            "solidified_center_surface_metric_nm": float(summary["center_surface_metric_after_solidification_nm"]),
            "solidified_mean_surface_metric_nm": float(summary["mean_surface_metric_after_solidification_nm"]),
            "maximum_connected_pool_nm3": float(summary["maximum_connected_pool_nm3"]),
            "maximum_velocity_over_run_m_s": float(summary["max_velocity_m_per_s"]),
        })

    if len(rows_out) != 4 or rows_out[0]["snapshot_force_z_logged_N"] != 0.:
        raise ValueError("unexpected scenario set or baseline force")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows_out[0]),
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows_out)
    print(f"wrote {args.output}")
    for r in rows_out:
        print(r["scenario"],
              f"Fz={r['snapshot_force_z_reconstructed_N']:.8e} N",
              f"history_delta={r['snapshot_force_relative_difference']:.3e}",
              f"center={r['source_radius_mean_surface_change_nm']:.6f} nm",
              f"center_vs_zero={r['source_radius_response_vs_zero_nm']:.6f} nm",
              f"Vmax={r['max_velocity_at_common_time_m_s']:.5f} m/s")


if __name__ == "__main__":
    main()
