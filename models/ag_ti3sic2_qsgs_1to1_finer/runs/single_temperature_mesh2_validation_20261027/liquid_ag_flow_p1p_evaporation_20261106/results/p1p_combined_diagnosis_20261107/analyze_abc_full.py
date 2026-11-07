#!/usr/bin/env python3
"""Post-process completed full A/B/C heat-hold-cool runs without rerunning them."""
import csv
import hashlib
import math
import re
import struct
from pathlib import Path

import numpy as np
from scipy import ndimage

HERE = Path(__file__).resolve().parent
SHORT = HERE / "abc_short"
FULL = HERE / "abc_full"
CASES = ("A", "B", "C")
RHO_AG = 10490.0
DX_NM = 0.3018574846
ATOM_AG_KG = 107.8682e-3 / 6.02214076e23


def read_snapshot(path):
    raw = path.read_bytes()
    head = struct.unpack_from("<8s6I5d", raw)
    magic, version, nx, ny, nz, nf, flags, t_ps, dx_m, x0, y0, z0 = head
    if magic != b"P1HSNAP1":
        raise ValueError(f"{path}: unexpected snapshot magic {magic!r}")
    if len(raw) != 72 + nf * nx * ny * nz * 4:
        raise ValueError(f"{path}: incomplete or malformed snapshot")
    fields = np.frombuffer(raw, dtype="<f4", offset=72).reshape(nf, nx, ny, nz)
    return dict(fields=fields, nx=nx, ny=ny, nz=nz, nf=nf,
                time_ps=t_ps, dx_m=dx_m, x0=x0, y0=y0, z0=z0)


def rows(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def surface_height(snap):
    f = snap["fields"][0]
    z = snap["z0"] + (np.arange(snap["nz"]) + 0.5) * snap["dx_m"]
    h = np.full(f.shape[:2], np.nan)
    for k in range(snap["nz"]):
        present = f[:, :, k] > 0.5
        h[present] = z[k] + 0.5 * snap["dx_m"]
    return h


def atom_count_equiv(mass_kg):
    return mass_kg / ATOM_AG_KG


def label_periodic_xy(mask):
    """Six-face labels with the solver's periodic x/y and bounded z topology."""
    labels, count = ndimage.label(mask, structure=ndimage.generate_binary_structure(3, 1))
    if count == 0:
        return labels, 0
    parent = np.arange(count + 1, dtype=np.int64)

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(int(a)), find(int(b))
        if ra != rb:
            parent[rb] = ra

    for aa, bb in ((labels[0], labels[-1]), (labels[:, 0], labels[:, -1])):
        joined = (aa > 0) & (bb > 0)
        for a, b in zip(aa[joined], bb[joined]):
            union(a, b)
    roots = np.array([find(i) for i in range(count + 1)], dtype=np.int64)
    unique = np.unique(roots[1:])
    compact = np.zeros(count + 1, dtype=np.int32)
    for new, root in enumerate(unique, 1):
        compact[roots == root] = new
    return compact[labels], len(unique)


def main():
    init_path = SHORT / "A" / "p1h_marangoni_on_initial_t000.000000_ps.bin"
    initial = read_snapshot(init_path)
    f0 = initial["fields"][2].astype(np.float64)
    h0 = surface_height(initial)
    dx_nm = initial["dx_m"] * 1e9
    cell_nm3 = dx_nm ** 3
    nx, ny, nz = initial["nx"], initial["ny"], initial["nz"]
    x = initial["x0"] + (np.arange(nx) + 0.5) * initial["dx_m"]
    y = initial["y0"] + (np.arange(ny) + 0.5) * initial["dx_m"]
    X, Y = np.meshgrid(x, y, indexing="ij")
    r = np.hypot(X - 9.58832e-9, Y - 9.65943950771472e-9) * 1e9

    comparison, radial_rows, surface_rows, cluster_rows, peak_component_rows = [], [], [], [], []
    peak_cluster_rows, peak_radial_rows = [], []
    for case in CASES:
        d = FULL / case
        history = rows(d / "p1h_history.csv")
        final_row = history[-1]
        peak_pool = max(history, key=lambda q: float(q["liquid_Ag_volume_nm3"]))
        peak_temp = max(history, key=lambda q: float(q["Tmax_K"]))
        dimension_rows = rows(d / "p1h_peak_pool_dimensions.csv")
        largest_cluster = max(dimension_rows, key=lambda q: float(q["largest_liquid_volume_nm3"]))
        hold_paths = sorted(d.glob("p1h_marangoni_on_hold_start_*_t*.bin"))
        hold_snapshot = read_snapshot(hold_paths[0]) if len(hold_paths) == 1 else None
        snap_paths = sorted(d.glob("p1h_marangoni_on_solidified_t*.bin"))
        if len(snap_paths) > 1:
            raise RuntimeError(f"{case}: expected at most one verified solidified snapshot, found {len(snap_paths)}")
        final = read_snapshot(snap_paths[0]) if snap_paths else None
        if final is not None and abs(final["time_ps"] - float(final_row["time_ps"])) > 2e-6:
            raise RuntimeError(f"{case}: final history/snapshot times do not match")
        if final is not None:
            final_fields = final["fields"]
            agf = final_fields[2].astype(np.float64)
            liqf = final_fields[4].astype(np.float64)
            total_ag_nm3 = float(agf.sum(dtype=np.float64) * cell_nm3)
            total_liquid_nm3 = float((agf * liqf).sum(dtype=np.float64) * cell_nm3)
        else:
            agf = liqf = None
            total_ag_nm3 = total_liquid_nm3 = math.nan

        if hold_snapshot is not None:
            peak_sample = hold_snapshot["fields"][0].astype(np.float64)
            peak_agf = hold_snapshot["fields"][2].astype(np.float64)
            peakliq = hold_snapshot["fields"][4].astype(np.float64)
            peak_theta = peak_sample * peak_agf * peakliq
            for threshold in (0.1, 0.5, 0.9):
                mask = peak_theta >= threshold
                labels, count = label_periodic_xy(mask)
                sizes = np.bincount(labels.ravel(), weights=peak_theta.ravel(), minlength=count + 1)[1:] * cell_nm3
                peak_cluster_rows.append({
                    "case": case,
                    "snapshot_role": "peak_pool_hold_start_snapshot",
                    "snapshot_time_ps": hold_snapshot["time_ps"],
                    "threshold": threshold,
                    "components_6_face_periodic_xy": int(count),
                    "largest_component_nm3": float(sizes.max()) if sizes.size else 0.0,
                    "thresholded_liquid_volume_nm3": float(sizes.sum()),
                })
                if threshold == 0.5:
                    ranked = sorted(range(1,count+1),key=lambda lab:sizes[lab-1],reverse=True)
                    for rank, lab in enumerate(ranked,1):
                        ii,jj,kk=np.where(labels==lab)
                        peak_component_rows.append({
                            "case": case,
                            "snapshot_time_ps": hold_snapshot["time_ps"],
                            "threshold": threshold,
                            "connectivity": "6-face-periodic-xy",
                            "component_rank_by_volume": rank,
                            "component_label": lab,
                            "cell_count": int(ii.size),
                            "volume_nm3_weighted_by_sample_Ag_liquid_fraction": float(sizes[lab-1]),
                            "x_min_nm": float(x[ii.min()]*1e9), "x_max_nm": float((x[ii.max()]+initial["dx_m"])*1e9),
                            "y_min_nm": float(y[jj.min()]*1e9), "y_max_nm": float((y[jj.max()]+initial["dx_m"])*1e9),
                            "z_min_nm": float((initial["z0"]+kk.min()*initial["dx_m"])*1e9),
                            "z_max_nm": float((initial["z0"]+(kk.max()+1)*initial["dx_m"])*1e9),
                        })
            for lo, hi in ((0, 3), (3, 5), (5, 8), (8, math.inf)):
                m3 = np.broadcast_to(((r >= lo) & (r < hi))[:, :, None], peak_agf.shape)
                v0 = float(f0[m3].sum(dtype=np.float64) * cell_nm3)
                v1 = float(peak_agf[m3].sum(dtype=np.float64) * cell_nm3)
                peak_radial_rows.append({
                    "case": case,
                    "snapshot_time_ps": hold_snapshot["time_ps"],
                    "radial_zone_nm": f"{lo}-{hi if math.isfinite(hi) else 'inf'}",
                    "Ag_volume_initial_nm3": v0,
                    "Ag_volume_at_hold_start_nm3": v1,
                    "delta_Ag_volume_nm3": v1 - v0,
                })

        if final is not None:
            # Surface displacement is relative to the common initial topography;
            # positive values denote downward depression.
            hf = surface_height(final)
            valid2 = np.isfinite(h0) & np.isfinite(hf)
            dep = (h0 - hf) * 1e9
            zones = ((0, 3), (3, 5), (5, 8), (8, math.inf))
            for lo, hi in zones:
                m2 = valid2 & (r >= lo) & (r < hi)
                v = dep[m2]
                surface_rows.append({
                    "case": case,
                    "final_time_ps": final["time_ps"],
                    "radial_zone_nm": f"{lo}-{hi if math.isfinite(hi) else 'inf'}",
                    "columns": int(m2.sum()),
                    "mean_depression_nm_positive_down": float(np.mean(v)) if v.size else math.nan,
                    "max_depression_nm_positive_down": float(np.max(v)) if v.size else math.nan,
                    "min_depression_nm_positive_down": float(np.min(v)) if v.size else math.nan,
                    "grid_spacing_nm": dx_nm,
                })

                m3 = np.broadcast_to(((r >= lo) & (r < hi))[:, :, None], f0.shape)
                v0 = float(f0[m3].sum(dtype=np.float64) * cell_nm3)
                v1 = float(agf[m3].sum(dtype=np.float64) * cell_nm3)
                dv = v1 - v0
                radial_rows.append({
                    "case": case,
                    "radial_zone_nm": f"{lo}-{hi if math.isfinite(hi) else 'inf'}",
                    "Ag_volume_initial_nm3": v0,
                    "Ag_volume_final_nm3": v1,
                    "delta_Ag_volume_nm3": dv,
                    "delta_Ag_mass_kg_at_10490kg_m3": dv * 1e-27 * RHO_AG,
                })

            for threshold in (0.1, 0.5, 0.9):
                theta = final_fields[0].astype(np.float64) * agf * liqf
                mask = theta >= threshold
                labels, count = label_periodic_xy(mask)
                sizes = np.bincount(labels.ravel(), weights=theta.ravel(), minlength=count + 1)[1:] * cell_nm3
                cluster_rows.append({
                    "case": case,
                    "threshold": threshold,
                    "components_6_face_periodic_xy": int(count),
                    "largest_component_nm3": float(sizes.max()) if sizes.size else 0.0,
                    "thresholded_liquid_volume_nm3": float(sizes.sum()) if sizes.size else 0.0,
                })

        stderr = (d / "full.stderr").read_text(errors="replace")
        hold_start = re.search(r"HEAT_HOLD_STARTED_[^\n]*t_ps=([0-9.eE+-]+)", stderr)
        hold_end = re.search(r"HEAT_OFF_AFTER_BALANCED_HOLD[^\n]*t_ps=([0-9.eE+-]+)[^\n]*hold_age_ps=([0-9.eE+-]+)[^\n]*hold_complete=(\d+)[^\n]*balance_command_mW=([0-9.eE+-]+)", stderr)
        guard_end = re.search(r"HEAT_OFF_PRESTEP_PROPERTY_GUARD[^\n]*t_ps=([0-9.eE+-]+)[^\n]*Tmax_skeleton_K=([0-9.eE+-]+)[^\n]*margin_K=([0-9.eE+-]+)[^\n]*ceiling_K=([0-9.eE+-]+)[^\n]*hold_started=(\d+)[^\n]*hold_age_ps=([0-9.eE+-]+)", stderr)
        if not guard_end:
            guard_end = re.search(r"HEAT_OFF_PROPERTY_GUARD[^\n]*t_ps=([0-9.eE+-]+)[^\n]*Tmax_skeleton_K=([0-9.eE+-]+)[^\n]*predicted_after_ramp_K=([0-9.eE+-]+)[^\n]*Tmax_Ag_K=([0-9.eE+-]+)[^\n]*liquid_Ag_nm3=([0-9.eE+-]+)[^\n]*ceiling_K=([0-9.eE+-]+)[^\n]*hold_started=(\d+)[^\n]*hold_age_ps=([0-9.eE+-]+)", stderr)
        guard_type = "PRESTEP_PROPERTY_GUARD" if "HEAT_OFF_PRESTEP_PROPERTY_GUARD" in stderr else "PROPERTY_GUARD" if "HEAT_OFF_PROPERTY_GUARD" in stderr else ""
        p1p_stop = re.search(r"STOP P1P property/composition failure t_ps=([0-9.eE+-]+) invalid_interface_cells=(\d+) table_range_K=([^;]+); ([^\n]+)", stderr)
        balance_rows = re.findall(
            r"P1M_BALANCED_HOLD t_ps=([0-9.eE+-]+) max_power_mW=([0-9.eE+-]+) commanded_mW=([0-9.eE+-]+) measured_input_mW=([0-9.eE+-]+) measured_sink_mW=([0-9.eE+-]+) net_mW=([0-9.eE+-]+)",
            stderr)
        last_balance = balance_rows[-1] if balance_rows else None
        arc_rows = rows(d / "p1h_arc_pressure_history.csv") if (d / "p1h_arc_pressure_history.csv").exists() else []
        max_arc = max((float(q.get("max_effective_interface_Pa", 0)) for q in arc_rows), default=0.0)
        max_force = max((abs(float(q.get("integrated_pressure_force_z_N", 0))) for q in arc_rows), default=0.0)
        evaporated = 0.0
        evap_path = d / "p1p_evaporation_flux.csv"
        if evap_path.exists():
            evap_rows = rows(evap_path)
            if evap_rows:
                evaporated = float(evap_rows[-1].get("escaped_mass_cumulative_kg", 0))
        launcher_status = (d / "full.exit_status").read_text(errors="replace").strip().replace("\n", "; ") if (d / "full.exit_status").exists() else ""
        stop_lines = [line for line in stderr.splitlines() if line.startswith("STOP ")]
        if stop_lines:
            outcome = "INCOMPLETE_STOP"
        elif "DONE_SOLIDIFIED" in stderr:
            outcome = "DONE_SOLIDIFIED"
        else:
            outcome = "IN_PROGRESS_OR_INCOMPLETE"

        comparison.append({
            "case": case,
            "start_snapshot_sha256": sha256(init_path),
            "continuation_checkpoint_sha256": sha256(d / "start_from_3ps.dump") if (d / "start_from_3ps.dump").exists() else "",
            "arc_configured": int(case in ("B", "C")),
            "evaporation_enabled": int(case == "C"),
            "marangoni_multiplier": 3,
            "peak_pool_time_ps": float(peak_pool["time_ps"]),
            "peak_liquid_Ag_volume_nm3": float(peak_pool["liquid_Ag_volume_nm3"]),
            "maximum_largest_cluster_time_ps": float(largest_cluster["time_ps"]),
            "maximum_largest_cluster_nm3": float(largest_cluster["largest_liquid_volume_nm3"]),
            "component_count_at_largest_cluster_time": int(largest_cluster["components"]),
            "largest_cluster_xyz_nm": "x={:.4g};y={:.4g};z={:.4g}".format(
                float(largest_cluster["x_size_nm"]), float(largest_cluster["y_size_nm"]),
                float(largest_cluster["z_size_nm"])),
            "hold_start_snapshot_time_ps": hold_snapshot["time_ps"] if hold_snapshot else math.nan,
            "hold_start_snapshot": hold_paths[0].name if hold_snapshot else "",
            "maximum_temperature_time_ps": float(peak_temp["time_ps"]),
            "maximum_temperature_K": float(peak_temp["Tmax_K"]),
            "maximum_skeleton_temperature_K": float(max(float(q["Tmax_skeleton_K"]) for q in history)),
            "hold_started_ps": float(hold_start.group(1)) if hold_start else math.nan,
            "hold_end_ps": (float(hold_end.group(1)) if hold_end else
                            float(guard_end.group(1)) if guard_end else
                            float(p1p_stop.group(1)) if p1p_stop else math.nan),
            "actual_hold_duration_ps": (float(hold_end.group(2)) if hold_end else
                                         float(guard_end.group(6)) if guard_end and guard_type == "PRESTEP_PROPERTY_GUARD" else
                                         float(guard_end.group(8)) if guard_end else
                                         float(p1p_stop.group(1)) - float(hold_start.group(1)) if p1p_stop and hold_start else math.nan),
            "hold_complete": int(hold_end.group(3)) if hold_end else 0,
            "hold_termination_reason": ("completed_5ps" if hold_end else guard_type if guard_type else
                                         "P1P_PROPERTY_COMPOSITION_STOP" if p1p_stop else "unknown_or_not_reached"),
            "guard_reported_skeleton_temperature_K": (float(guard_end.group(2)) if guard_end and guard_type == "PRESTEP_PROPERTY_GUARD" else
                                                       float(guard_end.group(2)) if guard_end else math.nan),
            "guard_property_ceiling_K": (float(guard_end.group(4)) if guard_end and guard_type == "PRESTEP_PROPERTY_GUARD" else
                                          float(guard_end.group(6)) if guard_end else math.nan),
            "balanced_source_power_mW": float(hold_end.group(4)) if hold_end else math.nan,
            "last_balanced_hold_time_ps": float(last_balance[0]) if last_balance else math.nan,
            "last_balanced_command_mW": float(last_balance[2]) if last_balance else math.nan,
            "last_balanced_input_mW": float(last_balance[3]) if last_balance else math.nan,
            "last_balanced_sink_mW": float(last_balance[4]) if last_balance else math.nan,
            "last_balanced_net_mW": float(last_balance[5]) if last_balance else math.nan,
            "last_valid_history_time_ps": float(final_row["time_ps"]),
            "last_valid_Tmax_K": float(final_row["Tmax_K"]),
            "last_valid_Tmax_Ag_K": float(final_row["Tmax_Ag_K"]),
            "last_valid_Tmax_skeleton_K": float(final_row["Tmax_skeleton_K"]),
            "last_valid_liquid_Ag_volume_nm3_history": float(final_row["liquid_Ag_volume_nm3"]),
            "final_liquid_Ag_volume_nm3_snapshot": total_liquid_nm3,
            "final_Ag_volume_nm3": total_ag_nm3,
            "final_Ag_volume_drift_pct": float(final_row["Ag_volume_drift_pct"]) if final is not None else math.nan,
            "final_energy_residual_J": float(final_row["energy_residual_J"]) if final is not None else math.nan,
            "final_input_energy_J": float(final_row["input_energy_J"]) if final is not None else math.nan,
            "final_bottom_sink_energy_J": float(final_row["bottom_sink_energy_J"]) if final is not None else math.nan,
            "run_outcome": outcome,
            "launcher_status_record": launcher_status,
            "solver_stop_reason": stop_lines[-1] if stop_lines else "",
            "max_effective_arc_pressure_Pa": max_arc,
            "max_abs_integrated_arc_force_z_N": max_force,
            "evaporated_mass_kg": evaporated,
            "evaporated_Ag_atom_mass_equivalents": atom_count_equiv(evaporated),
            "solidified_snapshot": snap_paths[0].name if snap_paths else "",
        })

    if len({r["start_snapshot_sha256"] for r in comparison}) != 1:
        raise RuntimeError("A/B/C were not started from matching snapshots")
    for filename, data in (
        ("P1P_ABC_FULL_FINAL_COMPARISON.csv", comparison),
        ("P1P_ABC_FULL_RADIAL_MASS.csv", radial_rows),
        ("P1P_ABC_FULL_SURFACE_SHAPE.csv", surface_rows),
        ("P1P_ABC_FULL_LIQUID_CLUSTERS.csv", cluster_rows),
        ("P1P_ABC_FULL_PEAK_CLUSTERS.csv", peak_cluster_rows),
        ("P1P_ABC_FULL_PEAK_RADIAL_MASS.csv", peak_radial_rows),
        ("P1P_ABC_FULL_PEAK_COMPONENTS.csv", peak_component_rows),
    ):
        with (HERE / filename).open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(data[0]))
            w.writeheader()
            w.writerows(data)
    for row in comparison:
        print(row)


if __name__ == "__main__":
    main()
