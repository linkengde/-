#!/usr/bin/env python3
"""Summarize existing P1H pressure-on, peak, pressure-off, and solidified fields."""
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
P167_DIR = "pressure_sensitivity_pi1/marangoni_on_gaussian_167MPa"


def read_frame(path: Path, magic_expected: bytes, fields_expected: int):
    with gzip.open(path, "rb") as stream:
        raw = stream.read()
    magic, version, nx, ny, nz, nf, case_flag, time_ps, dx_m, x0, y0, z0 = HEADER.unpack_from(raw)
    if magic != magic_expected or version != 1 or nf != fields_expected:
        raise ValueError(f"invalid P1H field header: {path}")
    expected = HEADER.size + 4*nf*nx*ny*nz
    if len(raw) != expected:
        raise ValueError(f"truncated P1H field file: {path}")
    data = np.frombuffer(raw, dtype="<f4", offset=HEADER.size).reshape(nf, nx, ny, nz)
    if not np.isfinite(data).all():
        raise ValueError(f"non-finite field values: {path}")
    return data, dict(time_ps=time_ps, dx_m=dx_m, origin_m=(x0, y0, z0),
                      shape=(nx, ny, nz), case_flag=case_flag)


def rows(path: Path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def surface_metrics(initial_vof, vof, meta):
    nx, ny, nz = vof.shape
    dx = meta["dx_m"]*1.e9
    x0, y0, z0 = (v*1.e9 for v in meta["origin_m"])
    z_top = z0 + (np.arange(nz) + 1.)*dx
    h0 = np.max(np.where(initial_vof > 0.5, z_top[None, None, :], -np.inf), axis=2)
    h = np.max(np.where(vof > 0.5, z_top[None, None, :], -np.inf), axis=2)
    change = h0 - h
    x = x0 + (np.arange(nx) + 0.5)*dx
    y = y0 + (np.arange(ny) + 0.5)*dx
    xx, yy = np.meshgrid(x, y, indexing="ij")
    r = np.hypot(xx-SOURCE_X_NM, yy-SOURCE_Y_NM)
    valid = np.isfinite(h0) & np.isfinite(h)
    center = valid & (r <= 5.)
    ring = valid & (r > 5.) & (r <= 8.)
    return {
        "center_mean_surface_change_nm": float(np.mean(change[center])),
        "ring_5_to_8nm_mean_surface_change_nm": float(np.mean(change[ring])),
        "center_max_downward_surface_change_nm": float(np.max(change[center])),
        "center_columns_below_initial_surface": int(np.count_nonzero(change[center] > 0.)),
    }, r


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    run_dir = args.run_root / P167_DIR
    pressure_history = rows(run_dir / "p1h_arc_pressure_history.csv")
    thermal_history = rows(run_dir / "p1h_history.csv")
    first_active = next(r for r in pressure_history
                        if float(r["max_effective_interface_Pa"]) > 1.)
    last_active = next(r for r in reversed(pressure_history)
                       if float(r["max_effective_interface_Pa"]) > 1.)
    first_zero_after_active = next(r for r in pressure_history
                                   if float(r["time_ps"]) > float(last_active["time_ps"])
                                   and float(r["max_effective_interface_Pa"]) <= 1.)

    events = [
        ("pressure_first_effective", float(first_active["time_ps"]),
         "regular_t000.263509", "nearest saved full field; first-load time has history only"),
        ("maximum_pool_and_heat_off", 0.797446025151503,
         "heat_off_t000.797446", "saved at event"),
        ("pressure_fully_off", 0.817446025151503,
         "regular_t001.014300", "first saved full field after pressure reached zero"),
        ("fully_solidified", 1.112870018554388,
         "solidified_t001.112870", "saved at event"),
    ]

    initial, initial_meta = read_frame(run_dir/"p1h_marangoni_on_initial_t000.000000_ps.bin.gz",
                                       b"P1HSNAP1", 20)
    initial_ag = np.maximum(initial[2], 0.)
    output_rows = []
    for event, event_time, suffix, note in events:
        main_path = run_dir/f"p1h_marangoni_on_{suffix}_ps.bin.gz"
        pressure_path = run_dir/f"p1h_pressure_marangoni_on_{suffix}_ps.bin.gz"
        fields, meta = read_frame(main_path, b"P1HSNAP1", 20)
        pressure, pmeta = read_frame(pressure_path, b"P1HPRS01", 4)
        if meta["shape"] != pmeta["shape"] or meta["time_ps"] != pmeta["time_ps"]:
            raise ValueError(f"main/pressure field time mismatch for {event}")
        if meta["shape"] != initial_meta["shape"]:
            raise ValueError(f"initial/follow-up shape mismatch for {event}")

        dx_nm = meta["dx_m"]*1.e9
        cell_volume_nm3 = dx_nm**3
        (surface, radial_distance) = surface_metrics(initial[0], fields[0], meta)
        region_center = radial_distance <= 5.
        region_ring = (radial_distance > 5.) & (radial_distance <= 8.)
        region_outer = radial_distance > 8.
        ag = np.maximum(fields[2], 0.)
        liquid_ag = np.maximum(fields[4], 0.)
        skeleton = np.maximum(fields[3], 0.)
        def region_volume(field, region):
            return float(np.sum(field[region, :])*cell_volume_nm3)

        speed = np.sqrt(fields[7]**2 + fields[8]**2 + fields[9]**2)
        liquid_weight = float(np.sum(liquid_ag))
        liquid_weighted_speed = (float(np.sum(speed*liquid_ag)/liquid_weight)
                                 if liquid_weight else 0.)
        interface = (fields[0] > 1.e-8) & (fields[0] < 1.-1.e-8)
        curvature_valid = interface & (np.abs(fields[11]) < 1.e20)
        if curvature_valid.any():
            curvature_median = float(np.median(np.abs(fields[11][curvature_valid]))/1.e9)
            gamma_median = float(np.median(fields[13][curvature_valid]))
        else:
            curvature_median = gamma_median = float("nan")
        pmax = float(np.max(pressure[2]))
        active_pressure_cells = pressure[2] > 0.
        capillary = pressure[3] - pressure[2]
        stress_selection = ((pressure[1] > 0.) & (pressure[2] >= 0.1*pmax)
                            if pmax else np.zeros_like(pressure[1], dtype=bool))
        if np.any(stress_selection):
            loaded_pressure = pressure[2][stress_selection]/1.e6
            hydrodynamic_pressure = fields[10][stress_selection]/1.e6
            loaded_curvature = np.abs(fields[11][stress_selection])/1.e9
            loaded_weighted_gamma = fields[13][stress_selection]
            pressure_capillary_median = float(np.median(pressure[2][stress_selection]/
                                                    (np.abs(capillary[stress_selection])+1.e-30)))
            loaded_pressure_median = float(np.median(loaded_pressure))
            loaded_pressure_p90 = float(np.quantile(loaded_pressure,0.9))
            hydro_pressure_p10 = float(np.quantile(hydrodynamic_pressure,0.1))
            hydro_pressure_median = float(np.median(hydrodynamic_pressure))
            hydro_pressure_p90 = float(np.quantile(hydrodynamic_pressure,0.9))
            loaded_curvature_median = float(np.median(loaded_curvature))
            loaded_gamma_median = float(np.median(loaded_weighted_gamma))
            capillary_median = float(np.median(np.abs(capillary[stress_selection]))/1.e6)
            capillary_p90 = float(np.quantile(np.abs(capillary[stress_selection]),0.9)/1.e6)
        else:
            loaded_pressure_median = loaded_pressure_p90 = float("nan")
            hydro_pressure_p10 = hydro_pressure_median = hydro_pressure_p90 = float("nan")
            loaded_curvature_median = loaded_gamma_median = float("nan")
            pressure_capillary_median = capillary_median = capillary_p90 = float("nan")

        # Match event times to the nearest pressure and thermal diagnostics,
        # but preserve both timestamps so no missing snapshot is hidden.
        p_row = min(pressure_history, key=lambda r: abs(float(r["time_ps"])-meta["time_ps"]))
        t_row = min(thermal_history, key=lambda r: abs(float(r["time_ps"])-meta["time_ps"]))
        if event == "pressure_first_effective":
            event_pressure_row = first_active
        elif event == "pressure_fully_off":
            event_pressure_row = first_zero_after_active
        else:
            event_pressure_row = min(pressure_history,
                                     key=lambda r: abs(float(r["time_ps"])-event_time))
        output_rows.append({
            "event": event,
            "event_time_ps": event_time,
            "saved_field_time_ps": meta["time_ps"],
            "event_to_field_offset_ps": meta["time_ps"]-event_time,
            "field_alignment_note": note,
            "main_field_file": main_path.name,
            "pressure_field_file": pressure_path.name,
            "event_pressure_history_time_ps": float(event_pressure_row["time_ps"]),
            "event_pressure_history_max_effective_Pa": float(event_pressure_row["max_effective_interface_Pa"]),
            "event_pressure_history_force_z_N": float(event_pressure_row["integrated_pressure_force_z_N"]),
            "pressure_off_time_is_smooth_ramp_end": event == "pressure_fully_off",
            "pressure_effective_max_Pa": pmax,
            "pressure_effective_support_cells": int(np.count_nonzero(active_pressure_cells)),
            "loaded_cells_median_effective_pressure_MPa": loaded_pressure_median,
            "loaded_cells_p90_effective_pressure_MPa": loaded_pressure_p90,
            "loaded_cells_hydrodynamic_pressure_p10_MPa": hydro_pressure_p10,
            "loaded_cells_hydrodynamic_pressure_median_MPa": hydro_pressure_median,
            "loaded_cells_hydrodynamic_pressure_p90_MPa": hydro_pressure_p90,
            "loaded_cells_median_abs_curvature_per_nm": loaded_curvature_median,
            "loaded_cells_median_weighted_surface_tension_N_m": loaded_gamma_median,
            "integrated_force_z_nearest_logged_N": float(p_row["integrated_pressure_force_z_N"]),
            "center_mean_surface_change_nm": surface["center_mean_surface_change_nm"],
            "ring_mean_surface_change_nm": surface["ring_5_to_8nm_mean_surface_change_nm"],
            "center_max_downward_surface_change_nm": surface["center_max_downward_surface_change_nm"],
            "center_columns_below_initial_surface": surface["center_columns_below_initial_surface"],
            "projected_Ag_marker_volume_sum_nm3_not_conservative": float(np.sum(ag)*cell_volume_nm3),
            "projected_Ag_marker_L1_change_equivalent_nm3": float(0.5*np.sum(np.abs(ag-initial_ag))*cell_volume_nm3),
            "Ag_center_volume_nm3": region_volume(ag, region_center),
            "Ag_center_delta_from_initial_nm3": region_volume(ag-initial_ag, region_center),
            "Ag_ring_5_to_8nm_volume_nm3": region_volume(ag, region_ring),
            "Ag_ring_delta_from_initial_nm3": region_volume(ag-initial_ag, region_ring),
            "Ag_outer_volume_nm3": region_volume(ag, region_outer),
            "projected_liquid_Ag_volume_sum_nm3_not_conservative": float(np.sum(liquid_ag)*cell_volume_nm3),
            "liquid_Ag_center_volume_nm3": region_volume(liquid_ag, region_center),
            "liquid_Ag_ring_volume_nm3": region_volume(liquid_ag, region_ring),
            "skeleton_volume_nm3": float(np.sum(skeleton)*cell_volume_nm3),
            "max_speed_field_m_s": float(np.max(speed)),
            "liquid_Ag_weighted_mean_speed_m_s": liquid_weighted_speed,
            "liquid_Ag_weighted_drag_q99_s-1": float(np.quantile(fields[18][liquid_ag>1.e-6],0.99)) if np.any(liquid_ag>1.e-6) else float("nan"),
            "liquid_Ag_fraction_drag_ge_1e10": float(np.sum(liquid_ag[fields[18]>=1.e10])/liquid_weight) if liquid_weight else float("nan"),
            "median_abs_capillary_pressure_on_loaded_cells_MPa": capillary_median,
            "p90_abs_capillary_pressure_on_loaded_cells_MPa": capillary_p90,
            "median_pressure_to_capillary_ratio_on_loaded_cells": pressure_capillary_median,
            "interface_median_abs_curvature_per_nm": curvature_median,
            "interface_median_weighted_surface_tension_N_m": gamma_median,
            "reported_max_velocity_m_s": float(t_row["max_velocity_m_s"]),
            "reported_liquid_Ag_volume_nm3": float(t_row["liquid_Ag_volume_nm3"]),
            "reported_Ag_volume_nm3": float(t_row["Ag_volume_nm3"]),
            "reported_source_center_metric_nm": float(t_row["source_center_depression_nm"]),
        })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]),
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"wrote {args.output}")
    print("first active pressure:", first_active["time_ps"],
          "ps; max effective pressure", first_active["max_effective_interface_Pa"], "Pa; Fz",
          first_active["integrated_pressure_force_z_N"], "N")
    print("last active pressure row:", last_active["time_ps"],
          "ps; first recorded zero row:", first_zero_after_active["time_ps"], "ps")
    for r in output_rows:
        print(r["event"], f"event={r['event_time_ps']:.9f} ps",
              f"field={r['saved_field_time_ps']:.9f} ps",
              f"pmax={r['pressure_effective_max_Pa']:.6g} Pa",
              f"center_change={r['center_mean_surface_change_nm']:.6g} nm",
              f"reported_liquid_Ag={r['reported_liquid_Ag_volume_nm3']:.6g} nm^3")


if __name__ == "__main__":
    main()
