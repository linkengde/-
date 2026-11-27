#!/usr/bin/env python3
"""对三套 150³ 几何执行局部厚度、弦长、界面及连通性审计。"""

from __future__ import annotations

import csv
import json
from pathlib import Path
import resource
import time

import numpy as np
from scipy import ndimage


ROOT = Path(__file__).resolve().parents[1]
GENERATION = ROOT / "geometry_generation"
OUT = ROOT / "geometry_audit"
OUT.mkdir(parents=True, exist_ok=True)
SIZE_BY_MODEL_M = {"ultrafine_500nm": 10e-6, "micro_10um": 200e-6}
SHAPE = (150, 150, 150)
CONNECTIVITY = ndimage.generate_binary_structure(3, 1)


def component_audit(mask: np.ndarray) -> dict[str, object]:
    labels, count = ndimage.label(mask, structure=CONNECTIVITY)
    sizes = np.bincount(labels.ravel(), minlength=count + 1)
    sizes[0] = 0
    largest = int(sizes.max(initial=0))
    largest_id = int(np.argmax(sizes)) if count else 0
    spans = []
    for axis in range(3):
        first = np.take(labels, 0, axis=axis)
        last = np.take(labels, -1, axis=axis)
        spans.append(bool(largest_id and np.any(first == largest_id) and np.any(last == largest_id)))
    return {
        "component_count_6face": int(count),
        "largest_component_fraction": largest / int(mask.sum()) if mask.any() else 0.0,
        "isolated_component_count": max(0, int(count) - (1 if count else 0)),
        "isolated_volume_fraction": (int(mask.sum()) - largest) / int(mask.sum()) if mask.any() else 0.0,
        "largest_component_spans_xyz": spans,
        "largest_component_mask": labels == largest_id if largest_id else np.zeros(mask.shape, dtype=bool),
    }


def maximal_inscribed_sphere_diameters(
    mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, dict[str, float | int]]:
    """提取 3×3×3 EDT 局部极大点对应的最大内接球直径。"""
    started = time.perf_counter()
    distance = ndimage.distance_transform_edt(mask)
    local_maximum = ndimage.maximum_filter(distance, size=3, mode="constant")
    centers = mask & (distance > 0.0) & (distance >= local_maximum - 1.0e-12)
    center_indices = np.flatnonzero(centers)
    diameters = 2.0 * distance.ravel()[center_indices]
    if not diameters.size:
        raise RuntimeError("最大内接球局部极大值集合为空")
    return diameters, center_indices, {
        "edt_seconds": time.perf_counter() - started,
        "maximal_inscribed_sphere_center_count": int(centers.sum()),
        "process_peak_rss_mib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
    }


def chord_lengths(mask: np.ndarray, axis: int, periodic: bool) -> np.ndarray:
    lines = np.moveaxis(mask, axis, -1).reshape(-1, mask.shape[axis])
    if periodic:
        full = np.all(lines, axis=1)
        first_false = np.argmax(~lines, axis=1)
        columns = (first_false[:, None] + np.arange(lines.shape[1])[None, :]) % lines.shape[1]
        ordered = np.take_along_axis(lines, columns, axis=1)
        ordered[full] = False
    else:
        full = np.zeros(lines.shape[0], dtype=bool)
        ordered = lines
    padded = np.pad(ordered.astype(np.int8), ((0, 0), (1, 1)))
    transitions = np.diff(padded, axis=1)
    starts = np.where(transitions == 1)
    ends = np.where(transitions == -1)
    lengths = ends[1] - starts[1]
    if np.any(full):
        lengths = np.concatenate((lengths, np.full(int(full.sum()), mask.shape[axis], dtype=np.int64)))
    return lengths.astype(np.int32, copy=False)


def percentile_summary(values: np.ndarray) -> dict[str, float]:
    quantiles = np.percentile(values, [1, 10, 50, 90, 99])
    return {
        "minimum": float(values.min()),
        "p01": float(quantiles[0]),
        "p10": float(quantiles[1]),
        "p50": float(quantiles[2]),
        "p90": float(quantiles[3]),
        "p99": float(quantiles[4]),
        "maximum": float(values.max()),
        "mean": float(values.mean()),
    }


def interface_area_density(mask: np.ndarray, spacing_m: float) -> float:
    dx = spacing_m
    transverse_area = sum(
        float(np.count_nonzero(np.diff(mask.astype(np.int8), axis=axis))) * dx * dx
        for axis in range(3)
    )
    volume = float(mask.size) * dx**3
    return transverse_area / volume


def surface_exposure(mask: np.ndarray) -> dict[str, float]:
    axes = (
        np.take(mask, 0, axis=0), np.take(mask, -1, axis=0),
        np.take(mask, 0, axis=1), np.take(mask, -1, axis=1),
        np.take(mask, 0, axis=2), np.take(mask, -1, axis=2),
    )
    face_fractions = [float(face.mean()) for face in axes]
    return {
        "skeleton_fraction_on_six_outer_faces": float(np.mean(face_fractions)),
        "ag_fraction_on_six_outer_faces": float(1.0 - np.mean(face_fractions)),
        "skeleton_face_fractions_xyz_minus_plus": face_fractions,
    }


def save_histogram(path: Path, group_rows: list[dict[str, object]]) -> None:
    keys = sorted({key for row in group_rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(group_rows)


def main() -> None:
    comparisons: list[dict[str, object]] = []
    thickness_histograms: list[dict[str, object]] = []
    chord_histograms: list[dict[str, object]] = []
    timing: list[dict[str, object]] = []
    for seed in (1105, 1136, 1182):
        archive = np.load(GENERATION / f"qsgs_seed_{seed}.npz")
        skeleton = np.asarray(archive["skeleton"], dtype=bool)
        ag = np.asarray(archive["ag"], dtype=bool)
        if skeleton.shape != SHAPE or not np.array_equal(ag, ~skeleton):
            raise ValueError(f"seed {seed}: phase geometry shape or complement mismatch")
        qsgs_started = time.perf_counter()
        skeleton_thickness, skeleton_centers, skeleton_timing = maximal_inscribed_sphere_diameters(skeleton)
        ag_thickness, ag_centers, ag_timing = maximal_inscribed_sphere_diameters(ag)
        thickness_seconds = time.perf_counter() - qsgs_started
        components_s = component_audit(skeleton)
        components_a = component_audit(ag)
        largest_s = components_s.pop("largest_component_mask")
        largest_a = components_a.pop("largest_component_mask")
        axes_data: dict[str, dict[str, np.ndarray]] = {}
        for phase_name, phase in (("skeleton", skeleton), ("ag", ag)):
            axes_data[phase_name] = {}
            for axis, axis_name in enumerate("xyz"):
                axes_data[phase_name][axis_name] = chord_lengths(phase, axis, periodic=False)

        per_seed: dict[str, object] = {
            "seed": seed,
            "shape_xyz": SHAPE,
            "skeleton_fraction": float(skeleton.mean()),
            "ag_fraction": float(ag.mean()),
            "qsgs_generations": int(archive["qsgs_generations"]),
            "six_face_connectivity": {
                "skeleton": components_s,
                "ag": components_a,
            },
            "maximal_inscribed_sphere_diameter_voxels": {
                "skeleton": percentile_summary(skeleton_thickness),
                "ag_channel": percentile_summary(ag_thickness),
            },
            "thickness_method": "3D Euclidean distance transform; identify all 26-neighbour local maxima in each phase as maximal-inscribed-sphere centers; report twice the center-to-background distance. Percentiles weight sphere centers equally, not phase voxels.",
            "thickness_method_diagnostics": {
                "skeleton": skeleton_timing,
                "ag": ag_timing,
                "combined_seconds": thickness_seconds,
            },
            "six_cut_plane_exposure": surface_exposure(skeleton),
            "chord_method": "binary run lengths along straight grid lines; all axes treated as open at the RVE boundary",
            "chord_voxels": {
                phase: {axis: percentile_summary(values) for axis, values in data.items()}
                for phase, data in axes_data.items()
            },
        }

        for model_name, size_m in SIZE_BY_MODEL_M.items():
            spacing_m = size_m / SHAPE[0]
            model = {
                "model": model_name,
                "size_m": [size_m] * 3,
                "spacing_m": [spacing_m] * 3,
                "spacing_nm": spacing_m * 1e9,
                "interface_area_density_m2_per_m3": interface_area_density(skeleton, spacing_m),
            }
            for phase_name, phase, thick, thick_centers, comp, main_comp in (
                ("skeleton", skeleton, skeleton_thickness, skeleton_centers, components_s, largest_s),
                ("ag", ag, ag_thickness, ag_centers, components_a, largest_a),
            ):
                exposed = surface_exposure(skeleton)
                model[f"{phase_name}_fraction"] = float(phase.mean())
                model[f"{phase_name}_largest_component_fraction"] = float(comp["largest_component_fraction"])
                model[f"{phase_name}_spans_xyz"] = comp["largest_component_spans_xyz"]
                unit_scale = 1e6
                thickness_summary = percentile_summary(thick)
                model[f"{phase_name}_thickness_um"] = {
                    key: value * spacing_m * 1e6 for key, value in thickness_summary.items()
                }
                main_sphere_diameters = thick[main_comp.ravel()[thick_centers]]
                model[f"{phase_name}_p01_thickness_voxels_main_component"] = float(
                    np.percentile(main_sphere_diameters, 1)
                )
                model[f"{phase_name}_minimum_sphere_diameter_voxels_main_component"] = float(
                    main_sphere_diameters.min()
                )
                if phase_name == "skeleton":
                    model["skeleton_exposed_fraction"] = exposed["skeleton_fraction_on_six_outer_faces"]
                for axis, values in axes_data[phase_name].items():
                    chord_um = values.astype(float) * spacing_m * unit_scale
                    chord_summary = percentile_summary(chord_um)
                    model.setdefault(f"{phase_name}_chord_um", {})[axis] = chord_summary
                    for p in (10, 50, 90, 99):
                        chord_histograms.append(
                            {
                                "seed": seed,
                                "model": model_name,
                                "phase": phase_name,
                                "axis": axis,
                                "quantile_percentile": p,
                                "quantile_um": chord_summary[f"p{p:02d}"],
                                "chord_count": int(values.size),
                            }
                        )
                hist_values, hist_edges = np.histogram(thick * spacing_m * 1e6, bins=60)
                for left, right, count in zip(hist_edges[:-1], hist_edges[1:], hist_values):
                    thickness_histograms.append(
                        {
                            "seed": seed,
                            "model": model_name,
                            "phase": phase_name,
                            "bin_left_um": left,
                            "bin_right_um": right,
                            "sphere_center_count": int(count),
                        }
                    )
            flat = {"seed": seed, "model": model_name}
            for key, value in model.items():
                if key.endswith("_thickness_um"):
                    phase_name = key.removesuffix("_thickness_um")
                    flat.update({f"{phase_name}_sphere_diameter_{stat}_um": val for stat, val in value.items()})
                elif key.endswith("_chord_um"):
                    phase_name = key.removesuffix("_chord_um")
                    for axis, stats in value.items():
                        flat.update({f"{phase_name}_chord_{axis}_{stat}_um": val for stat, val in stats.items()})
                elif isinstance(value, (dict, list)):
                    flat[key] = json.dumps(value, separators=(",", ":"))
                else:
                    flat[key] = value
            comparisons.append(flat)

        per_seed["interface_area_density_scale_pair_m2_per_m3"] = {
            model_name: interface_area_density(skeleton, size_m / SHAPE[0])
            for model_name, size_m in SIZE_BY_MODEL_M.items()
        }
        per_seed["interface_area_scale_ratio"] = (
            per_seed["interface_area_density_scale_pair_m2_per_m3"]["ultrafine_500nm"]
            / per_seed["interface_area_density_scale_pair_m2_per_m3"]["micro_10um"]
        )
        (OUT / f"seed_{seed}_audit.json").write_text(
            json.dumps(per_seed, indent=2) + "\n", encoding="utf-8"
        )
        timing.append(
            {
                "seed": seed,
                "thickness_seconds_both_phases": thickness_seconds,
                "skeleton_edT_seconds": skeleton_timing["edt_seconds"],
                "ag_edT_seconds": ag_timing["edt_seconds"],
                "peak_rss_mib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
            }
        )
        print(json.dumps(timing[-1]), flush=True)

    save_histogram(OUT / "geometry_comparison.csv", comparisons)
    save_histogram(OUT / "local_thickness_histogram.csv", thickness_histograms)
    save_histogram(OUT / "chord_distribution_summary.csv", chord_histograms)
    (OUT / "audit_timing.json").write_text(json.dumps(timing, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
