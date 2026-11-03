#!/usr/bin/env python3
"""Extract the largest tracked liquid-Ag pool and periodic local windows.

Read-only inputs are the completed 2026-11-02 regrown-QSGS Case-B checkpoint
and its paired geometry. This script writes derived artifacts only into this
directory; it never changes the source geometry or checkpoint.
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import distance_transform_edt, label as ndimage_label


HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
SOURCE_RUN = PARENT / "ag70vol_1x_regrown_64x64x50_caseB_4x_composite_sc2300K_20261102"
GEOMETRY_PATH = SOURCE_RUN / "geometry/qsgs_morphology_ag70vol_regrown_64x64x50.npz"
CHECKPOINT_PATH = SOURCE_RUN / "checkpoints/heating_complete_2300K.npz"
CLUSTERS_PATH = SOURCE_RUN / "liquid_clusters.csv"
PREFLIGHT_PATH = SOURCE_RUN / "caseB_preflight.json"
EXPECTED_GEOMETRY_SHA256 = "5e509f16716e85464f604fcda2a53eaefd4f15e2d18d6c0bf24e0aa06c6656ce"
EXPECTED_CHECKPOINT_SHA256 = "8b15a8b4766ad1d2fc8e5a6ec3d927e1fbe3728563ddc3f65546bd2a446406c6"
TIME_EXPECTED_PS = 0.8495860511825916
CLASS_TOL = 1.0e-12
WINDOW_PADS = (4, 8, 12)
THRESHOLDS = (0.1, 0.5, 0.9, 0.99)

sys.path.insert(0, str(SOURCE_RUN / "runtime"))
from qsgs_liquid_connectivity import liquid_cluster_labels  # noqa: E402


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_array(array: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def periodic_interval(indices: np.ndarray, count: int) -> tuple[int, int, bool]:
    """Return minimal cyclic interval start, cell width, and seam crossing."""
    unique = np.unique(indices.astype(np.int64))
    if unique.size == 0:
        raise ValueError("cannot bound an empty pool")
    if unique.size == count:
        return 0, count, True
    gaps = np.diff(np.concatenate((unique, unique[:1] + count)))
    gap_index = int(np.argmax(gaps))
    start = int(unique[(gap_index + 1) % unique.size])
    width = int(count - gaps[gap_index] + 1)
    raw_width = int(unique[-1] - unique[0] + 1)
    crosses_seam = raw_width != width
    return start, width, crosses_seam


def axis_bounds(pool: np.ndarray, spacing_a: np.ndarray) -> dict[str, object]:
    coords = np.where(pool)
    result: dict[str, object] = {}
    shape = pool.shape
    for axis, name in enumerate(("x", "y")):
        start, width, crosses = periodic_interval(coords[axis], shape[axis])
        unwrapped = (coords[axis] - start) % shape[axis]
        low = int(unwrapped.min())
        high = int(unwrapped.max())
        result[name] = {
            "start_index": start,
            "width_cells": width,
            "crosses_periodic_seam": crosses,
            "unwrapped_low_cell": low,
            "unwrapped_high_cell_inclusive": high,
            "size_nm": float(width * spacing_a[axis] / 10.0),
            "unwrapped_min_coordinate_from_interval_origin_A": float(low * spacing_a[axis]),
            "unwrapped_max_coordinate_from_interval_origin_A_exclusive": float((high + 1) * spacing_a[axis]),
            "source_grid_min_cell": int(coords[axis].min()),
            "source_grid_max_cell_inclusive": int(coords[axis].max()),
            "source_grid_min_coordinate_A": float(coords[axis].min() * spacing_a[axis]),
            "source_grid_max_coordinate_A_exclusive": float((coords[axis].max() + 1) * spacing_a[axis]),
        }
    z_low = int(coords[2].min())
    z_high = int(coords[2].max())
    result["z"] = {
        "start_index": z_low,
        "width_cells": z_high - z_low + 1,
        "touches_low_domain_face": z_low == 0,
        "touches_high_domain_face": z_high == shape[2] - 1,
        "size_nm": float((z_high - z_low + 1) * spacing_a[2] / 10.0),
        "source_grid_min_cell": z_low,
        "source_grid_max_cell_inclusive": z_high,
        "source_grid_min_coordinate_A": float(z_low * spacing_a[2]),
        "source_grid_max_coordinate_A_exclusive": float((z_high + 1) * spacing_a[2]),
    }
    return result


def shift_neighbor(array: np.ndarray, axis: int, shift: int, periodic: bool) -> np.ndarray:
    shifted = np.roll(array, shift=shift, axis=axis)
    if not periodic:
        edge = [slice(None)] * array.ndim
        edge[axis] = 0 if shift == 1 else -1
        shifted[tuple(edge)] = False
    return shifted


def interface_masks(skeleton: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    ag = ~skeleton
    ag_interface = np.zeros_like(ag)
    scaffold_interface = np.zeros_like(skeleton)
    for axis in range(3):
        periodic = axis in (0, 1)
        for shift in (-1, 1):
            neighbor_scaffold = shift_neighbor(skeleton, axis, shift, periodic)
            neighbor_ag = shift_neighbor(ag, axis, shift, periodic)
            ag_interface |= ag & neighbor_scaffold
            scaffold_interface |= skeleton & neighbor_ag
    return ag_interface, scaffold_interface


def face_areas_and_pool_contact(pool: np.ndarray, skeleton: np.ndarray, spacing_a: np.ndarray) -> dict[str, object]:
    per_normal: dict[str, dict[str, float | int]] = {}
    normals = ("x", "y", "z")
    for axis, normal in enumerate(normals):
        contacts = 0
        for shift in (-1, 1):
            neighbor_skeleton = shift_neighbor(skeleton, axis, shift, axis in (0, 1))
            contacts += int(np.count_nonzero(pool & neighbor_skeleton))
        face_area_a2 = float(np.prod(np.delete(spacing_a, axis)))
        per_normal[normal] = {
            "pool_scaffold_face_count": contacts,
            "voxel_face_area_A2": face_area_a2,
            "estimated_interface_area_nm2": contacts * face_area_a2 / 100.0,
        }
    return {
        "by_face_normal": per_normal,
        "estimated_total_interface_area_nm2": float(
            sum(float(item["estimated_interface_area_nm2"]) for item in per_normal.values())
        ),
        "method": "six-face anisotropic voxel-face area; x/y periodic seams included; z is nonperiodic",
    }


def component_summary(mask: np.ndarray, spacing_a: np.ndarray) -> dict[str, object]:
    labels, count = liquid_cluster_labels(mask, periodic_axes=(0, 1))
    sizes = np.bincount(labels.ravel(), minlength=count + 1)
    if count == 0:
        return {"component_count": 0, "active_voxels": 0, "largest_voxels": 0}
    largest_label = int(np.argmax(sizes[1:]) + 1)
    largest = labels == largest_label
    bounds = axis_bounds(largest, spacing_a)
    return {
        "component_count": int(count),
        "active_voxels": int(mask.sum()),
        "largest_component_label_at_threshold": largest_label,
        "largest_voxels": int(sizes[largest_label]),
        "largest_volume_nm3": float(sizes[largest_label] * np.prod(spacing_a) / 1000.0),
        "largest_bbox_nm_xyz": [float(bounds[axis]["size_nm"]) for axis in ("x", "y", "z")],
        "largest_touches_z_top": bool(bounds["z"]["touches_high_domain_face"]),
    }


def widest_spanning_path_throats(pool: np.ndarray, spacing_a: np.ndarray) -> dict[str, object]:
    """Voxel-limited EDT throat estimates for paths across the pool bbox."""
    structure = np.zeros((3, 3, 3), dtype=bool)
    structure[1, 1, 1] = True
    for axis in range(3):
        lo = [1, 1, 1]
        hi = [1, 1, 1]
        lo[axis] = 0
        hi[axis] = 2
        structure[tuple(lo)] = True
        structure[tuple(hi)] = True
    coordinates = np.where(pool)
    clearance = distance_transform_edt(pool, sampling=spacing_a)
    results = {}
    min_spacing = float(np.min(spacing_a))
    for axis, name in enumerate(("x", "y", "z")):
        low_face = int(coordinates[axis].min())
        high_face = int(coordinates[axis].max())
        critical_radius = 0.0
        for radius in np.unique(clearance[pool])[::-1]:
            eroded_core = pool & (clearance >= radius)
            labels, _ = ndimage_label(eroded_core, structure=structure)
            low_labels = np.unique(np.take(labels, low_face, axis=axis))
            high_labels = np.unique(np.take(labels, high_face, axis=axis))
            common = np.intersect1d(low_labels[low_labels > 0], high_labels[high_labels > 0])
            if common.size:
                critical_radius = float(radius)
                break
        diameter_a = max(0.0, 2.0 * critical_radius - min_spacing)
        results[name] = {
            "critical_edt_clearance_radius_A": critical_radius,
            "voxel_corrected_widest_path_bottleneck_diameter_A": diameter_a,
            "voxel_corrected_widest_path_bottleneck_diameter_nm": diameter_a / 10.0,
            "estimated_voxels_across_using_axis_spacing": diameter_a / float(spacing_a[axis]),
            "bbox_faces_connected_after_edt_threshold": bool(critical_radius > 0.0),
        }
    return {
        "method": "largest EDT-clearance threshold retaining a 6-face path between opposite faces of pool bbox; voxel-corrected diameter = max(0, 2*r - min voxel spacing)",
        "estimates_by_axis": results,
        "resolution_warning": "These are one-voxel-scale geometric estimates, not resolved hydrodynamic channel widths.",
    }


def pool_core_sensitivity(pool: np.ndarray, spacing_a: np.ndarray) -> list[dict[str, object]]:
    """Count connected cores after geometric clearance thresholds.

    This is a voxel-resolution sensitivity diagnostic, not a physical erosion
    law or hydrodynamic channel measurement.
    """
    structure = np.zeros((3, 3, 3), dtype=bool)
    structure[1, 1, 1] = True
    for axis in range(3):
        lo = [1, 1, 1]
        hi = [1, 1, 1]
        lo[axis] = 0
        hi[axis] = 2
        structure[tuple(lo)] = True
        structure[tuple(hi)] = True
    coords = np.where(pool)
    clearance = distance_transform_edt(pool, sampling=spacing_a / 10.0)
    rows = []
    for radius_nm in (0.0, 0.15, 0.30, 0.45, 0.60, 0.90, 1.20, 1.50):
        core = pool.copy() if radius_nm == 0.0 else pool & (clearance >= radius_nm)
        labels, count = ndimage_label(core, structure=structure)
        sizes = np.bincount(labels.ravel(), minlength=count + 1)
        largest_label = int(np.argmax(sizes[1:]) + 1) if count else 0
        spans = {}
        for axis, name in enumerate(("x", "y", "z")):
            if count == 0:
                spans[name] = False
                continue
            low_face, high_face = int(coords[axis].min()), int(coords[axis].max())
            low_labels = np.unique(np.take(labels, low_face, axis=axis))
            high_labels = np.unique(np.take(labels, high_face, axis=axis))
            spans[name] = bool(np.intersect1d(low_labels[low_labels > 0], high_labels[high_labels > 0]).size)
        rows.append({
            "edt_clearance_threshold_nm": radius_nm,
            "remaining_pool_voxels": int(core.sum()),
            "remaining_fraction_of_threshold_0p5_pool": float(core.sum() / pool.sum()),
            "six_face_connected_components": int(count),
            "largest_core_voxels": int(sizes[largest_label]) if largest_label else 0,
            "largest_core_spans_pool_bbox_faces_xyz": spans,
        })
    return rows


def recorded_cluster_history_audit(path: Path, heat_off_time_ps: float) -> dict[str, object]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    by_time: dict[float, list[dict[str, str]]] = {}
    for row in rows:
        # CSV floating-point serialization can make one frame appear as e.g.
        # 0.3 and 0.30000000000000004. Group at sub-femtosecond precision.
        by_time.setdefault(round(float(row["time_ps"]), 9), []).append(row)
    observed_max_by_time = []
    for time_ps, frame_rows in sorted(by_time.items()):
        largest = max(frame_rows, key=lambda row: int(row["voxel_count"]))
        observed_max_by_time.append({
            "time_ps": time_ps,
            "largest_cluster_voxels": int(largest["voxel_count"]),
            "largest_cluster_track_id": int(largest["track_id"]),
            "largest_cluster_component_id": int(largest["component_id"]),
        })
    heating = [row for row in observed_max_by_time if row["time_ps"] <= heat_off_time_ps + 1e-12]
    if not heating:
        raise RuntimeError("no cluster observations found during the heating interval")
    max_heating = max(heating, key=lambda row: row["largest_cluster_voxels"])
    max_all = max(observed_max_by_time, key=lambda row: row["largest_cluster_voxels"])
    return {
        "cluster_tracking_observation_count": len(observed_max_by_time),
        "cluster_tracking_observations_during_heating": len(heating),
        "heating_sample_times_ps": [row["time_ps"] for row in heating],
        "maximum_recorded_heating_cluster": max_heating,
        "maximum_recorded_cluster_across_heating_and_cooling_records": max_all,
        "heat_off_event_time_ps": heat_off_time_ps,
        "reported_sample_interval_ps": 0.1,
        "sampling_basis": "The run script updates LiquidClusterTracker when writing 0.1 ps series rows and at the scaffold-temperature event; it does not label every internal integration step.",
        "maximum_is_global_over_saved_tracker_observations": bool(
            round(float(max_heating["time_ps"]), 9) == round(heat_off_time_ps, 9)
        ),
        "continuous_time_global_maximum_certified": False,
        "limit": "This audit proves the maximum among saved tracker observations, not an unsampled transient peak between those observations.",
        "all_observed_frame_maxima": observed_max_by_time,
    }


def orthogonal_slice_data(
    pool: np.ndarray,
    phase_code: np.ndarray,
    spacing_a: np.ndarray,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Describe centroid-crossing slices in xy, xz, and yz planes."""
    centroid_idx = np.asarray(np.where(pool), dtype=float).mean(axis=1)
    center_indices = np.rint(centroid_idx).astype(int)
    plan = (
        ("xy", 2, int(center_indices[2]), (0, 1)),
        ("xz", 1, int(center_indices[1]), (0, 2)),
        ("yz", 0, int(center_indices[0]), (1, 2)),
    )
    summaries = {
        "pool_centroid_index_xyz": centroid_idx.tolist(),
        "slice_center_index_xyz": center_indices.tolist(),
        "planes": {},
    }
    count_rows = []
    for name, fixed_axis, index, display_axes in plan:
        pool_slice = np.take(pool, index, axis=fixed_axis)
        phase_slice = np.take(phase_code, index, axis=fixed_axis)
        counts = {
            "Ti3SiC2_scaffold": int(np.count_nonzero(phase_slice == 0)),
            "solid_Ag": int(np.count_nonzero(phase_slice == 1)),
            "partially_melted_Ag": int(np.count_nonzero(phase_slice == 2)),
            "fully_liquid_Ag": int(np.count_nonzero(phase_slice == 3)),
            "maximum_pool_f_l_ge_0p5": int(pool_slice.sum()),
            "slice_total_cells": int(pool_slice.size),
        }
        axis_names = ("x", "y", "z")
        summaries["planes"][name] = {
            "fixed_axis": axis_names[fixed_axis],
            "fixed_index": index,
            "fixed_coordinate_A_cell_center": float((index + 0.5) * spacing_a[fixed_axis]),
            "display_axes": [axis_names[axis] for axis in display_axes],
            "phase_and_pool_counts": counts,
        }
        count_rows.append({
            "plane": name,
            "fixed_axis": axis_names[fixed_axis],
            "fixed_index": index,
            "fixed_coordinate_A_cell_center": float((index + 0.5) * spacing_a[fixed_axis]),
            **counts,
        })
    return summaries, count_rows


def write_slice_figure(
    arrays: dict[str, np.ndarray],
    source_x: np.ndarray,
    source_y: np.ndarray,
    source_z: np.ndarray,
    center_indices: list[int],
    spacing_a: np.ndarray,
    time_ps: float,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    prefix = "pad8_"
    phase = arrays[prefix + "phase_code_0_scaffold_1_solid_ag_2_partial_ag_3_fully_liquid_ag"]
    pool = arrays[prefix + "max_pool_mask"]
    local_indices = {
        axis: int(np.flatnonzero(source == center_indices[axis])[0])
        for axis, source in enumerate((source_x, source_y, source_z))
    }
    centers_nm = {
        "x": arrays[prefix + "x_unwrapped_center_A"] / 10.0,
        "y": arrays[prefix + "y_unwrapped_center_A"] / 10.0,
        "z": arrays[prefix + "z_center_A"] / 10.0,
    }
    plans = (
        ("XY", "z", local_indices[2], ("x", "y"), phase[:, :, local_indices[2]], pool[:, :, local_indices[2]]),
        ("XZ", "y", local_indices[1], ("x", "z"), phase[:, local_indices[1], :], pool[:, local_indices[1], :]),
        ("YZ", "x", local_indices[0], ("y", "z"), phase[local_indices[0], :, :], pool[local_indices[0], :, :]),
    )
    cmap = ListedColormap(("#17232b", "#84b8d8", "#f4df00", "#f47721"))
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.7), facecolor="white", constrained_layout=True)
    for ax, (plane, fixed_axis, local_fixed, displayed_axes, field, pool_field) in zip(axes, plans):
        horizontal, vertical = displayed_axes
        xcenters, ycenters = centers_nm[horizontal], centers_nm[vertical]
        dx = float(spacing_a["xyz".index(horizontal)] / 10.0)
        dy = float(spacing_a["xyz".index(vertical)] / 10.0)
        extent = (xcenters[0] - dx / 2, xcenters[-1] + dx / 2,
                  ycenters[0] - dy / 2, ycenters[-1] + dy / 2)
        ax.imshow(field.T, origin="lower", interpolation="nearest", aspect="equal",
                  extent=extent, cmap=cmap, vmin=-0.5, vmax=3.5)
        if pool_field.any():
            ax.contour(xcenters, ycenters, pool_field.T.astype(float), levels=(0.5,),
                       colors=("#c00028",), linewidths=1.2)
        original_index = center_indices["xyz".index(fixed_axis)]
        fixed_a = (original_index + 0.5) * spacing_a["xyz".index(fixed_axis)]
        ax.set_title(f"{plane} 截面｜{fixed_axis} 索引 {original_index} ({fixed_a/10:.3f} nm)")
        ax.set_xlabel(f"{horizontal} (nm)")
        ax.set_ylabel(f"{vertical} (nm)")
    legend = [
        Patch(facecolor="#17232b", label="Ti3SiC2 骨架"),
        Patch(facecolor="#84b8d8", label="固态 Ag"),
        Patch(facecolor="#f4df00", label="固液共存 Ag"),
        Patch(facecolor="#f47721", label="完全液态 Ag"),
        Patch(facecolor="none", edgecolor="#c00028", label="f_l >= 0.5 最大池边界"),
    ]
    fig.legend(handles=legend, loc="upper center", ncol=5, frameon=False, bbox_to_anchor=(0.5, 1.07))
    fig.suptitle(f"Ag70 最大连通液相簇中心正交截面｜t={time_ps:.6f} ps｜网格体素视图", y=1.13)
    fig.savefig(HERE / "maximum_pool_orthogonal_slices.png", dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def read_source() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, object], dict[str, object]]:
    geometry_sha = sha256_file(GEOMETRY_PATH)
    checkpoint_sha = sha256_file(CHECKPOINT_PATH)
    if geometry_sha != EXPECTED_GEOMETRY_SHA256:
        raise RuntimeError(f"geometry SHA mismatch: {geometry_sha}")
    if checkpoint_sha != EXPECTED_CHECKPOINT_SHA256:
        raise RuntimeError(f"checkpoint SHA mismatch: {checkpoint_sha}")

    with np.load(GEOMETRY_PATH, allow_pickle=False) as archive:
        skeleton = np.asarray(archive["skeleton"], dtype=bool)
        domain_size_a = np.asarray(archive["size_A"], dtype=float)
        seed = int(archive["seed"])
        qsgs_seed = int(archive["qsgs_growth_seed"])

    with np.load(CHECKPOINT_PATH, allow_pickle=False) as archive:
        temperature = np.asarray(archive["temperature_K"], dtype=float)
        liquid_fraction = np.asarray(archive["ag_liquid_fraction"], dtype=float)
        enthalpy_density = np.asarray(archive["enthalpy_density_J_m3"], dtype=float)
        tracker_labels = np.asarray(archive["tracker_track_labels"], dtype=np.int32)
        metadata = json.loads(str(archive["metadata_json"].item()))
    preflight = json.loads(PREFLIGHT_PATH.read_text(encoding="utf-8"))

    if not all(a.shape == skeleton.shape for a in (temperature, liquid_fraction, enthalpy_density, tracker_labels)):
        raise RuntimeError("checkpoint arrays and source geometry shapes differ")
    if abs(float(metadata["time_ps"]) - TIME_EXPECTED_PS) > 1e-12:
        raise RuntimeError(f"unexpected checkpoint time {metadata['time_ps']}")
    if metadata["identity"]["geometry_sha256"] != geometry_sha:
        raise RuntimeError("checkpoint geometry identity does not match the geometry file")

    source = {
        "geometry_sha256": geometry_sha,
        "checkpoint_sha256": checkpoint_sha,
        "checkpoint_time_ps": float(metadata["time_ps"]),
        "grid_shape_xyz": list(skeleton.shape),
        "domain_size_A_xyz": domain_size_a.tolist(),
        "spacing_A_xyz": (domain_size_a / np.asarray(skeleton.shape)).tolist(),
        "qsgs_seed_field": seed,
        "qsgs_growth_seed": qsgs_seed,
        "source_run_branch": preflight["configuration"].get("git_branch"),
        "source_run_commit": preflight["configuration"].get("git_commit"),
        "source_run_configuration": metadata["identity"],
        "thermal_model": "single-temperature enthalpy Case B; this run used a composite (Ag + Ti3SiC2) 4x source",
        "source_run_geometry_note": "fresh 64x64x50 QSGS realization; not a crop or the earlier 100x100x50 geometry",
    }
    return skeleton, temperature, liquid_fraction, enthalpy_density, tracker_labels, metadata, source


def extract() -> None:
    skeleton, temperature, fl, enthalpy, tracker_labels, metadata, source = read_source()
    ag = ~skeleton
    shape = np.asarray(ag.shape, dtype=int)
    spacing_a = np.asarray(source["spacing_A_xyz"], dtype=float)
    time_ps = float(metadata["time_ps"])
    history = metadata["tracker_metadata"]["history"]
    at_time = [
        row for row in history
        if abs(float(row["time_ps"]) - time_ps) < 1e-12 and bool(row["is_main_cluster"])
    ]
    if len(at_time) != 1:
        raise RuntimeError(f"expected one main cluster record at checkpoint time, found {len(at_time)}")
    main_record = at_time[0]
    pool_track_id = int(main_record["track_id"])
    pool = tracker_labels == pool_track_id
    threshold_05 = ag & (fl >= float(metadata["tracker_metadata"]["liquid_threshold"]))
    labels_05, _ = liquid_cluster_labels(threshold_05, periodic_axes=(0, 1))
    threshold_largest = labels_05 == int(np.argmax(np.bincount(labels_05.ravel())[1:]) + 1)
    if not np.array_equal(pool, threshold_largest):
        raise RuntimeError("persistent main track does not equal the largest current 0.5-threshold component")
    if int(pool.sum()) != int(main_record["voxel_count"]):
        raise RuntimeError("persistent track count and checkpoint history disagree")

    with CLUSTERS_PATH.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    csv_rows = [
        row for row in rows
        if abs(float(row["time_ps"]) - time_ps) < 1e-12
        and int(row["track_id"]) == pool_track_id
        and row["is_main_cluster"].lower() == "true"
    ]
    if len(csv_rows) != 1 or int(csv_rows[0]["voxel_count"]) != int(pool.sum()):
        raise RuntimeError("pool identity/count does not match the saved cluster-tracking CSV")

    pool_bounds = axis_bounds(pool, spacing_a)
    full_liquid = ag & (fl >= 1.0 - CLASS_TOL)
    partial = ag & (fl > CLASS_TOL) & (fl < 1.0 - CLASS_TOL)
    solid_ag = ag & (fl <= CLASS_TOL)
    liquid_dominant = ag & (fl >= 0.5)
    phase_code = np.zeros(ag.shape, dtype=np.uint8)
    phase_code[solid_ag] = 1
    phase_code[partial] = 2
    phase_code[full_liquid] = 3
    ag_interface, scaffold_interface = interface_masks(skeleton)

    threshold_results = []
    for threshold in THRESHOLDS:
        result = component_summary(ag & (fl >= threshold), spacing_a)
        result["liquid_fraction_threshold_inclusive"] = threshold
        threshold_results.append(result)

    history_audit = recorded_cluster_history_audit(CLUSTERS_PATH, time_ps)
    core_sensitivity = pool_core_sensitivity(pool, spacing_a)
    orthogonal_slices, orthogonal_slice_counts = orthogonal_slice_data(
        pool, phase_code, spacing_a
    )

    spacing_nm = spacing_a / 10.0
    pool_idx = np.where(pool)
    pool_bbox_size_cells = [int(pool_bounds[a]["width_cells"]) for a in ("x", "y", "z")]
    bbox_x_start = int(pool_bounds["x"]["start_index"])
    bbox_y_start = int(pool_bounds["y"]["start_index"])
    bbox_z_start = int(pool_bounds["z"]["start_index"])
    bbox_x_width, bbox_y_width, bbox_z_width = pool_bbox_size_cells
    x_bbox_idx = (bbox_x_start + np.arange(bbox_x_width)) % shape[0]
    y_bbox_idx = (bbox_y_start + np.arange(bbox_y_width)) % shape[1]
    z_bbox_idx = bbox_z_start + np.arange(bbox_z_width)
    bbox_tuple = np.ix_(x_bbox_idx, y_bbox_idx, z_bbox_idx)
    bbox_scaffold_count = int(skeleton[bbox_tuple].sum())
    bbox_ag_count = int(ag[bbox_tuple].sum())
    bbox_pool_count = int(pool[bbox_tuple].sum())

    window_arrays: dict[str, np.ndarray] = {}
    window_summaries = []
    for pad in WINDOW_PADS:
        x_start_raw = bbox_x_start - pad
        y_start_raw = bbox_y_start - pad
        x_len = min(int(shape[0]), bbox_x_width + 2 * pad)
        y_len = min(int(shape[1]), bbox_y_width + 2 * pad)
        x_index = (x_start_raw + np.arange(x_len)) % shape[0]
        y_index = (y_start_raw + np.arange(y_len)) % shape[1]
        z_low = max(0, bbox_z_start - pad)
        z_high_exclusive = min(int(shape[2]), bbox_z_start + bbox_z_width + pad)
        z_index = np.arange(z_low, z_high_exclusive)
        index = np.ix_(x_index, y_index, z_index)
        prefix = f"pad{pad}_"
        arrays = {
            "skeleton_mask": skeleton,
            "ag_mask": ag,
            "max_pool_mask": pool,
            "liquid_dominant_ag_mask_f_ge_0p5": liquid_dominant,
            "fully_liquid_ag_mask": full_liquid,
            "partially_melted_ag_mask": partial,
            "solid_ag_mask": solid_ag,
            "ag_interface_mask": ag_interface,
            "scaffold_interface_mask": scaffold_interface,
            "phase_code_0_scaffold_1_solid_ag_2_partial_ag_3_fully_liquid_ag": phase_code,
            "temperature_K": temperature,
            "enthalpy_density_J_m3": enthalpy,
            "ag_liquid_fraction": fl,
        }
        for name, array in arrays.items():
            window_arrays[prefix + name] = np.asarray(array[index])
        window_arrays[prefix + "source_x_index"] = x_index.astype(np.int32)
        window_arrays[prefix + "source_y_index"] = y_index.astype(np.int32)
        window_arrays[prefix + "source_z_index"] = z_index.astype(np.int32)
        window_arrays[prefix + "x_unwrapped_center_A"] = (
            (x_start_raw + np.arange(x_len) + 0.5) * spacing_a[0]
        )
        window_arrays[prefix + "y_unwrapped_center_A"] = (
            (y_start_raw + np.arange(y_len) + 0.5) * spacing_a[1]
        )
        window_arrays[prefix + "z_center_A"] = (z_index + 0.5) * spacing_a[2]

        local_pool = window_arrays[prefix + "max_pool_mask"]
        global_x_seam = bool(np.any(pool[0, :, :] & pool[-1, :, :]))
        global_y_seam = bool(np.any(pool[:, 0, :] & pool[:, -1, :]))
        touches = {
            "x_local_low_or_high": bool(local_pool[0, :, :].any() or local_pool[-1, :, :].any()),
            "y_local_low_or_high": bool(local_pool[:, 0, :].any() or local_pool[:, -1, :].any()),
            "z_local_low_or_high": bool(local_pool[:, :, 0].any() or local_pool[:, :, -1].any()),
            "pool_crosses_global_x_periodic_seam": global_x_seam,
            "pool_crosses_global_y_periodic_seam": global_y_seam,
            "pool_touches_global_z_top_physical_face": bool(pool[:, :, -1].any()),
            "pool_touches_artificial_lower_z_crop_face": bool(z_low > 0 and local_pool[:, :, 0].any()),
        }
        window_summaries.append({
            "requested_padding_voxels_each_lateral_side_and_z_side": pad,
            "window_shape_xyz": [int(x_len), int(y_len), int(len(z_index))],
            "source_x_index_sequence": x_index.tolist(),
            "source_y_index_sequence": y_index.tolist(),
            "source_z_index_range_inclusive": [int(z_index[0]), int(z_index[-1])],
            "effective_z_lower_padding_voxels": bbox_z_start - z_low,
            "effective_z_upper_padding_voxels": z_high_exclusive - (bbox_z_start + bbox_z_width),
            "z_upper_padding_missing_due_to_original_nonperiodic_domain_face": int(
                max(0, pad - (z_high_exclusive - (bbox_z_start + bbox_z_width)))
            ),
            "phase_cell_counts": {
                "Ti3SiC2_scaffold": int(window_arrays[prefix + "skeleton_mask"].sum()),
                "fully_solid_Ag": int(window_arrays[prefix + "solid_ag_mask"].sum()),
                "partially_melted_Ag": int(window_arrays[prefix + "partially_melted_ag_mask"].sum()),
                "fully_liquid_Ag": int(window_arrays[prefix + "fully_liquid_ag_mask"].sum()),
                "liquid_dominant_Ag_f_ge_0p5": int(window_arrays[prefix + "liquid_dominant_ag_mask_f_ge_0p5"].sum()),
                "maximum_pool_track": int(local_pool.sum()),
            },
            "other_liquid_dominant_cells_in_window": int(
                window_arrays[prefix + "liquid_dominant_ag_mask_f_ge_0p5"].sum() - local_pool.sum()
            ),
            "pool_boundary_contacts": touches,
        })

    # Full-domain derived masks are a compact convenience; source fields remain untouched.
    np.savez_compressed(
        HERE / "maximum_pool_full_domain_masks.npz",
        skeleton_mask=skeleton,
        ag_mask=ag,
        maximum_pool_track_mask=pool,
        liquid_dominant_ag_mask_f_ge_0p5=liquid_dominant,
        fully_liquid_ag_mask=full_liquid,
        partially_melted_ag_mask=partial,
        solid_ag_mask=solid_ag,
        ag_interface_mask=ag_interface,
        scaffold_interface_mask=scaffold_interface,
        phase_code=phase_code,
    )
    window_arrays["manifest_json"] = np.asarray(json.dumps({
        "source": source,
        "pool_track_id": pool_track_id,
        "pool_bbox": pool_bounds,
        "windows": window_summaries,
        "array_sha256": {
            "skeleton": sha256_array(skeleton),
            "temperature_K": sha256_array(temperature),
            "ag_liquid_fraction": sha256_array(fl),
            "enthalpy_density_J_m3": sha256_array(enthalpy),
            "tracker_track_labels": sha256_array(tracker_labels),
        },
    }))
    np.savez_compressed(HERE / "maximum_pool_local_windows.npz", **window_arrays)

    metrics = {
        "status": "READ_ONLY_MAXIMUM_LIQUID_AG_POOL_EXTRACTION_COMPLETE",
        "source": source,
        "checkpoint_metadata_time_ps": time_ps,
        "pool_tracking": {
            "main_track_id": pool_track_id,
            "component_id_at_checkpoint": int(main_record["component_id"]),
            "centroid_A_from_tracker": main_record["centroid_A"],
            "parent_track_ids_at_last_merge": main_record["parent_track_ids"],
            "continuous_tracked_duration_ps": float(main_record["continuous_duration_ps"]),
            "liquid_threshold_inclusive": float(metadata["tracker_metadata"]["liquid_threshold"]),
            "pool_voxels": int(pool.sum()),
            "pool_volume_nm3": float(pool.sum() * np.prod(spacing_nm)),
            "pool_mean_liquid_fraction": float(fl[pool].mean()),
            "pool_min_liquid_fraction": float(fl[pool].min()),
            "pool_max_liquid_fraction": float(fl[pool].max()),
            "pool_temperature_min_K": float(temperature[pool].min()),
            "pool_temperature_mean_K": float(temperature[pool].mean()),
            "pool_temperature_max_K": float(temperature[pool].max()),
            "pool_fully_liquid_voxels": int(np.count_nonzero(pool & full_liquid)),
            "pool_partially_melted_voxels": int(np.count_nonzero(pool & partial)),
            "pool_bbox": pool_bounds,
            "widest_spanning_path_throat_estimates": widest_spanning_path_throats(pool, spacing_a),
            "clearance_core_sensitivity": core_sensitivity,
            "centroid_crossing_orthogonal_slices": orthogonal_slices,
            "saved_cluster_history_audit": history_audit,
            "threshold_sensitivity": threshold_results,
            "all_threshold_0p5_component_count": int(component_summary(liquid_dominant, spacing_a)["component_count"]),
            "all_threshold_0p5_liquid_dominant_voxels": int(liquid_dominant.sum()),
            "other_threshold_0p5_voxels_outside_main_pool": int(liquid_dominant.sum() - pool.sum()),
            "full_domain_phase_cell_counts": {
                "Ti3SiC2_scaffold": int(skeleton.sum()),
                "fully_solid_Ag": int(solid_ag.sum()),
                "partially_melted_Ag": int(partial.sum()),
                "fully_liquid_Ag": int(full_liquid.sum()),
                "Ag_total": int(ag.sum()),
                "liquid_dominant_Ag_f_ge_0p5": int(liquid_dominant.sum()),
            },
            "pool_bbox_envelope_voxel_count": int(np.prod(pool_bbox_size_cells)),
            "pool_bbox_envelope_material_counts": {
                "skeleton_voxels": bbox_scaffold_count,
                "Ag_voxels": bbox_ag_count,
                "maximum_pool_voxels": bbox_pool_count,
                "pool_share_of_bbox_envelope": bbox_pool_count / float(np.prod(pool_bbox_size_cells)),
            },
            "pool_scaffold_interface_face_estimate": face_areas_and_pool_contact(pool, skeleton, spacing_a),
            "global_periodic_seam_contact": {
                "x": bool(np.any(pool[0, :, :] & pool[-1, :, :])),
                "y": bool(np.any(pool[:, 0, :] & pool[:, -1, :])),
            },
            "touches_nonperiodic_z_top_boundary": bool(pool[:, :, -1].any()),
            "touches_nonperiodic_z_bottom_boundary": bool(pool[:, :, 0].any()),
        },
        "phase_classification": {
            "fully_solid_Ag": f"f_l <= {CLASS_TOL:g}",
            "partially_melted_Ag": f"{CLASS_TOL:g} < f_l < {1.0-CLASS_TOL:g}",
            "fully_liquid_Ag": f"f_l >= {1.0-CLASS_TOL:g}",
            "connected_pool_tracking_rule": "Ag AND f_l >= 0.5, six-face connectivity, periodic x/y, nonperiodic z",
        },
        "buffer_window_sensitivity": window_summaries,
        "limitations": [
            "The exact 12.29 x 11.77 x 4.12 nm pool belongs to a 64x64x50 fresh QSGS geometry, not the stated earlier 100x100x50 geometry.",
            "The pool touches the original nonperiodic z-top face; no physical solid buffer exists above it. The extracted windows preserve this original face rather than inventing material beyond the sample.",
            "The checkpoint is a single-temperature enthalpy Case-B result with a 4x composite source, not an atomistic structure and not a flow calculation.",
            "The source run extrapolated Ti3SiC2 thermophysical correlations to a 2300 K stop event and assumed an intact scaffold; decomposition is not represented.",
            "The box envelope contains solid Ag and Ti3SiC2; it is not itself the liquid pool. Pool identity uses the recorded f_l >= 0.5 threshold.",
        ],
    }
    (HERE / "maximum_pool_metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")

    with (HERE / "threshold_sensitivity.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, lineterminator="\n", fieldnames=[
            "liquid_fraction_threshold_inclusive", "component_count", "active_voxels",
            "largest_component_label_at_threshold", "largest_voxels", "largest_volume_nm3",
            "largest_bbox_nm_xyz", "largest_touches_z_top",
        ])
        writer.writeheader()
        for row in threshold_results:
            writer.writerow({key: json.dumps(value) if isinstance(value, list) else value for key, value in row.items()})

    with (HERE / "buffer_sensitivity.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, lineterminator="\n", fieldnames=[
            "requested_padding_voxels_each_lateral_side_and_z_side", "window_shape_xyz",
            "window_size_nm_xyz", "effective_z_lower_padding_voxels",
            "effective_z_upper_padding_voxels", "z_upper_padding_missing_due_to_original_nonperiodic_domain_face",
            "scaffold_voxels", "solid_ag_voxels", "partial_ag_voxels", "fully_liquid_ag_voxels",
            "liquid_dominant_ag_voxels", "maximum_pool_voxels", "other_liquid_dominant_cells_in_window",
            "pool_touches_x_or_y_artificial_cut", "pool_touches_z_original_top_face",
        ])
        for row in window_summaries:
            xsize = row["window_shape_xyz"][0] * spacing_nm[0]
            ysize = row["window_shape_xyz"][1] * spacing_nm[1]
            zsize = row["window_shape_xyz"][2] * spacing_nm[2]
            counts = row["phase_cell_counts"]
            touch = row["pool_boundary_contacts"]
            writer.writerow({
                "requested_padding_voxels_each_lateral_side_and_z_side": row["requested_padding_voxels_each_lateral_side_and_z_side"],
                "window_shape_xyz": json.dumps(row["window_shape_xyz"]),
                "window_size_nm_xyz": json.dumps([float(xsize), float(ysize), float(zsize)]),
                "effective_z_lower_padding_voxels": row["effective_z_lower_padding_voxels"],
                "effective_z_upper_padding_voxels": row["effective_z_upper_padding_voxels"],
                "z_upper_padding_missing_due_to_original_nonperiodic_domain_face": row["z_upper_padding_missing_due_to_original_nonperiodic_domain_face"],
                "scaffold_voxels": counts["Ti3SiC2_scaffold"],
                "solid_ag_voxels": counts["fully_solid_Ag"],
                "partial_ag_voxels": counts["partially_melted_Ag"],
                "fully_liquid_ag_voxels": counts["fully_liquid_Ag"],
                "liquid_dominant_ag_voxels": counts["liquid_dominant_Ag_f_ge_0p5"],
                "maximum_pool_voxels": counts["maximum_pool_track"],
                "other_liquid_dominant_cells_in_window": row["other_liquid_dominant_cells_in_window"],
                "pool_touches_x_or_y_artificial_cut": bool(touch["x_local_low_or_high"] or touch["y_local_low_or_high"]),
                "pool_touches_z_original_top_face": touch["pool_touches_global_z_top_physical_face"],
            })

    with (HERE / "pool_core_sensitivity.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, lineterminator="\n", fieldnames=[
            "edt_clearance_threshold_nm", "remaining_pool_voxels",
            "remaining_fraction_of_threshold_0p5_pool", "six_face_connected_components",
            "largest_core_voxels", "largest_core_spans_pool_bbox_faces_xyz",
        ])
        writer.writeheader()
        for row in core_sensitivity:
            writer.writerow({
                **row,
                "largest_core_spans_pool_bbox_faces_xyz": json.dumps(row["largest_core_spans_pool_bbox_faces_xyz"], sort_keys=True),
            })

    with (HERE / "orthogonal_slice_counts.csv").open("w", newline="", encoding="utf-8") as stream:
        fields = [
            "plane", "fixed_axis", "fixed_index", "fixed_coordinate_A_cell_center",
            "Ti3SiC2_scaffold", "solid_Ag", "partially_melted_Ag", "fully_liquid_Ag",
            "maximum_pool_f_l_ge_0p5", "slice_total_cells",
        ]
        writer = csv.DictWriter(stream, lineterminator="\n", fieldnames=fields)
        writer.writeheader()
        writer.writerows(orthogonal_slice_counts)

    create_figure(window_arrays, window_summaries[1], pool_track_id, time_ps)
    write_slice_figure(
        window_arrays,
        window_arrays["pad8_source_x_index"],
        window_arrays["pad8_source_y_index"],
        window_arrays["pad8_source_z_index"],
        orthogonal_slices["slice_center_index_xyz"],
        spacing_a,
        time_ps,
    )
    print(f"pool_track={pool_track_id} time_ps={time_ps:.12f} count={pool.sum()}")
    print(f"bbox_nm_xyz={[pool_bounds[a]['size_nm'] for a in ('x','y','z')]}")
    print(f"source_grid_xyz={list(shape)} geometry_sha256={source['geometry_sha256']}")


def create_figure(arrays: dict[str, np.ndarray], window: dict[str, object], track_id: int, time_ps: float) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    from matplotlib.patches import Patch
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    font_path = Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
    if font_path.exists():
        font_manager.fontManager.addfont(str(font_path))
        plt.rcParams["font.sans-serif"] = [font_manager.FontProperties(fname=str(font_path)).get_name(), "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    prefix = "pad8_"
    skeleton = arrays[prefix + "skeleton_mask"]
    solid = arrays[prefix + "solid_ag_mask"]
    partial = arrays[prefix + "partially_melted_ag_mask"]
    full = arrays[prefix + "fully_liquid_ag_mask"]
    pool = arrays[prefix + "max_pool_mask"]
    x = arrays[prefix + "x_unwrapped_center_A"] / 10.0
    y = arrays[prefix + "y_unwrapped_center_A"] / 10.0
    z = arrays[prefix + "z_center_A"] / 10.0
    coords = np.meshgrid(x, y, z, indexing="ij")
    point = np.stack(coords, axis=-1)
    fig = plt.figure(figsize=(13.8, 10.2), facecolor="white")
    ax = fig.add_subplot(111, projection="3d")
    ax.set_proj_type("ortho")
    for mask, color, alpha, size in (
        (solid, "#84b8d8", 0.22, 7),
        (skeleton, "#17232b", 0.72, 8),
        (partial, "#f4df00", 0.95, 15),
        (full, "#f47721", 0.55, 16),
    ):
        ijk = np.argwhere(mask)
        if ijk.size:
            pts = point[mask]
            ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], s=size, c=color, alpha=alpha,
                      linewidths=0, depthshade=False, rasterized=True)
    # Draw only the maximum-pool envelope; it is explicitly labeled as an envelope, not liquid fill.
    pool_points_nm = point[pool]
    dx_nm = (x[1] - x[0]) if x.size > 1 else 0.0
    dy_nm = (y[1] - y[0]) if y.size > 1 else 0.0
    dz_nm = (z[1] - z[0]) if z.size > 1 else 0.0
    xlo, xhi = float(pool_points_nm[:, 0].min() - dx_nm / 2.0), float(pool_points_nm[:, 0].max() + dx_nm / 2.0)
    ylo, yhi = float(pool_points_nm[:, 1].min() - dy_nm / 2.0), float(pool_points_nm[:, 1].max() + dy_nm / 2.0)
    zlo, zhi = float(pool_points_nm[:, 2].min() - dz_nm / 2.0), float(pool_points_nm[:, 2].max() + dz_nm / 2.0)
    corners = np.array([[xlo,ylo,zlo],[xhi,ylo,zlo],[xhi,yhi,zlo],[xlo,yhi,zlo],
                        [xlo,ylo,zhi],[xhi,ylo,zhi],[xhi,yhi,zhi],[xlo,yhi,zhi]])
    edges = ((0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7))
    for a, b in edges:
        ax.plot(*corners[[a,b]].T, color="#c00028", linewidth=1.1, linestyle="--", alpha=0.85)
    ax.set_xlim(float(x.min()), float(x.max()))
    ax.set_ylim(float(y.min()), float(y.max()))
    ax.set_zlim(float(z.min()), float(z.max()))
    ax.set_box_aspect((float(x.max()-x.min()), float(y.max()-y.min()), float(z.max()-z.min())))
    ax.view_init(elev=27, azim=-54)
    ax.set_axis_off()
    legend = [
        Patch(facecolor="#17232b", label="Ti3SiC2 骨架"),
        Patch(facecolor="#84b8d8", alpha=0.6, label="固态 Ag"),
        Patch(facecolor="#f4df00", label="固液共存 Ag"),
        Patch(facecolor="#f47721", alpha=0.55, label="完全液态 Ag"),
        Patch(facecolor="none", edgecolor="#c00028", linestyle="--", label="最大池包围盒（不是液体填充区）"),
    ]
    ax.legend(handles=legend, loc="upper left", frameon=False, fontsize=10)
    ax.set_title(
        f"Ag70 最大连通液相团簇｜t={time_ps:.6f} ps｜track {track_id}｜阈值 f_l>=0.5\n"
        "局部窗口：侧向缓冲 8 体素；z 保留原域顶面",
        pad=15, fontsize=14,
    )
    fig.text(0.04, 0.035,
             "源网格为 64×64×50（2026-11-02 新 QSGS），非早期 100×100×50 网格；球/点表示连续体网格体素，不代表原子。",
             ha="left", fontsize=9)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.88, bottom=0.08)
    fig.savefig(HERE / "maximum_pool_local_3d.png", dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    extract()
