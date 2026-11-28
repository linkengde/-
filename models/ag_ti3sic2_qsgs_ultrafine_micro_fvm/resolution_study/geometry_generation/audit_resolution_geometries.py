#!/usr/bin/env python3
"""审计200³/300³几何并与已提交150³基准比较。"""

from __future__ import annotations

import csv
import json
from pathlib import Path
import resource
import sys
import time

import numpy as np
from scipy import ndimage


HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
PROJECT = STUDY.parent
sys.path.insert(0, str(PROJECT / "geometry_audit"))
from audit_geometries import chord_lengths, interface_area_density, percentile_summary, surface_exposure  # noqa: E402


MODEL_SIZES_M = {"ultrafine_500nm": 10.0e-6, "micro_10um": 200.0e-6}
STRUCTURE = ndimage.generate_binary_structure(3, 1)
PHASES = (("skeleton", True), ("ag", False))


def component_metrics(mask: np.ndarray) -> tuple[dict[str, object], np.ndarray]:
    labels, count = ndimage.label(mask, structure=STRUCTURE)
    sizes = np.bincount(labels.ravel(), minlength=count + 1)
    sizes[0] = 0
    largest_id = int(np.argmax(sizes)) if count else 0
    largest_size = int(sizes[largest_id]) if count else 0
    largest_mask = labels == largest_id if count else np.zeros(mask.shape, dtype=bool)
    spans = []
    for axis in range(3):
        spans.append(bool(largest_id and np.any(np.take(largest_mask, 0, axis=axis))
                          and np.any(np.take(largest_mask, -1, axis=axis))))
    result = {
        "component_count_6face": int(count),
        "largest_component_fraction": largest_size / int(mask.sum()) if mask.any() else 0.0,
        "isolated_component_count": max(0, int(count) - 1),
        "isolated_volume_fraction": (int(mask.sum()) - largest_size) / int(mask.sum()) if mask.any() else 0.0,
        "largest_component_spans_xyz": spans,
    }
    del labels, sizes
    return result, largest_mask


def phase_metrics(mask: np.ndarray, largest_mask: np.ndarray) -> dict[str, object]:
    distance = ndimage.distance_transform_edt(mask)
    neighborhood = ndimage.maximum_filter(distance, size=3, mode="constant")
    centers = mask & (distance > 0.0) & (distance >= neighborhood - 1.0e-12)
    center_indices = np.flatnonzero(centers)
    diameters_vox = 2.0 * distance.ravel()[center_indices]
    main = largest_mask.ravel()[center_indices]
    del neighborhood, centers
    return {
        "sphere_diameter_voxels": percentile_summary(diameters_vox),
        "main_p01_sphere_diameter_voxels": float(np.percentile(diameters_vox[main], 1)),
        "main_min_sphere_diameter_voxels": float(diameters_vox[main].min()),
        "sphere_center_count": int(diameters_vox.size),
    }


def audit_one(path: Path) -> dict[str, object]:
    started = time.perf_counter()
    with np.load(path, allow_pickle=False) as f:
        skeleton = np.asarray(f["skeleton"], dtype=bool)
        n = int(f["shape_xyz"][0])
        seed = int(f["seed"])
        p0 = float(f["p0"])
        growth = np.asarray(f["growth_xyz"], dtype=float)
        sigma_vox = float(f["sigma_voxels"])
        sigma_fraction = float(f["sigma_fraction_of_domain"])
        generations = int(f["growth_generations"])
    if skeleton.shape != (n, n, n) or int(skeleton.sum()) != round(n**3 * 0.30):
        raise ValueError(f"Unexpected geometry shape/fraction in {path}")
    ag = ~skeleton
    components: dict[str, dict[str, object]] = {}
    largest_masks: dict[str, np.ndarray] = {}
    thickness: dict[str, dict[str, object]] = {}
    for name, mask in (("skeleton", skeleton), ("ag", ag)):
        components[name], largest_masks[name] = component_metrics(mask)
        thickness[name] = phase_metrics(mask, largest_masks[name])

    chords: dict[str, dict[str, dict[str, float]]] = {}
    for phase_name, phase in (("skeleton", skeleton), ("ag", ag)):
        chords[phase_name] = {}
        for axis, axis_name in enumerate("xyz"):
            chords[phase_name][axis_name] = percentile_summary(chord_lengths(phase, axis, periodic=False).astype(float))

    rows = []
    for model_name, size_m in MODEL_SIZES_M.items():
        spacing_m = size_m / n
        row: dict[str, object] = {
            "n": n,
            "seed": seed,
            "model": model_name,
            "size_m": json.dumps([size_m] * 3, separators=(",", ":")),
            "spacing_m": json.dumps([spacing_m] * 3, separators=(",", ":")),
            "p0": p0,
            "growth_xyz": json.dumps(growth.tolist(), separators=(",", ":")),
            "sigma_voxels": sigma_vox,
            "sigma_fraction_of_domain": sigma_fraction,
            "sigma_m": size_m * sigma_fraction,
            "qsgs_growth_generations": generations,
            "skeleton_fraction": float(skeleton.mean()),
            "ag_fraction": float(ag.mean()),
            "interface_area_density_m2_per_m3": interface_area_density(skeleton, spacing_m),
            "skeleton_exposed_fraction": surface_exposure(skeleton)["skeleton_fraction_on_six_outer_faces"],
        }
        for phase in ("skeleton", "ag"):
            comp = components[phase]
            thick = thickness[phase]
            row[f"{phase}_component_count_6face"] = comp["component_count_6face"]
            row[f"{phase}_largest_component_fraction"] = comp["largest_component_fraction"]
            row[f"{phase}_isolated_component_count"] = comp["isolated_component_count"]
            row[f"{phase}_isolated_volume_fraction"] = comp["isolated_volume_fraction"]
            row[f"{phase}_spans_xyz"] = json.dumps(comp["largest_component_spans_xyz"], separators=(",", ":"))
            row[f"{phase}_sphere_center_count"] = thick["sphere_center_count"]
            row[f"{phase}_minimum_sphere_diameter_voxels_main_component"] = thick["main_min_sphere_diameter_voxels"]
            row[f"{phase}_p01_sphere_diameter_voxels_main_component"] = thick["main_p01_sphere_diameter_voxels"]
            for stat, value in thick["sphere_diameter_voxels"].items():
                row[f"{phase}_sphere_diameter_{stat}_um"] = value * spacing_m * 1.0e6
            for axis, stats in chords[phase].items():
                for stat, value in stats.items():
                    row[f"{phase}_chord_{axis}_{stat}_um"] = value * spacing_m * 1.0e6
        rows.append(row)

    elapsed = time.perf_counter() - started
    return {
        "seed": seed,
        "n": n,
        "source_file": str(path.relative_to(STUDY)),
        "parameters": {"p0": p0, "growth_xyz": growth.tolist(), "sigma_voxels": sigma_vox,
                       "sigma_fraction_of_domain": sigma_fraction, "growth_generations": generations},
        "skeleton_fraction": float(skeleton.mean()),
        "phases_6face_connectivity": components,
        "thickness_method": "3D EDT local maxima in 26-neighborhood; maximum-inscribed-sphere diameter; center-weighted, plus main-component minimum and P01",
        "thickness_voxels": {phase: thickness[phase]["sphere_diameter_voxels"] for phase, _ in PHASES},
        "interface_area_density_m2_per_m3_by_model": {
            model: next(r["interface_area_density_m2_per_m3"] for r in rows if r["model"] == model)
            for model in MODEL_SIZES_M
        },
        "audit_wall_seconds": elapsed,
        "process_peak_rss_mib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0,
        "comparison_rows": rows,
    }


def baseline_rows() -> list[dict[str, object]]:
    source = PROJECT / "GEOMETRY_COMPARISON.csv"
    with source.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    result = []
    for row in rows:
        if row["model"] not in MODEL_SIZES_M:
            continue
        row["n"] = 150
        row["p0"] = 0.05
        row["growth_xyz"] = "[0.11,0.11,0.11]"
        row["sigma_voxels"] = 4.0
        row["sigma_fraction_of_domain"] = 4 / 150
        row["sigma_m"] = MODEL_SIZES_M[row["model"]] * (4 / 150)
        result.append(row)
    return result


def main() -> None:
    audit_started = time.perf_counter()
    details = []
    rows: list[dict[str, object]] = baseline_rows()
    for path in sorted((STUDY / "geometries").glob("n*/qsgs_seed_*.npz")):
        result = audit_one(path)
        rows.extend(result.pop("comparison_rows"))
        details.append(result)
        print(json.dumps({"seed": result["seed"], "n": result["n"],
                          "audit_wall_seconds": result["audit_wall_seconds"],
                          "process_peak_rss_mib": result["process_peak_rss_mib"]}), flush=True)

    output = STUDY / "geometry_audit"
    output.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row})
    with (output / "resolution_comparison.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "status": "RESOLUTION_GEOMETRY_AUDIT_COMPLETED",
        "grid_sizes": [150, 200, 300],
        "seeds": [1105, 1136, 1182],
        "physical_model_sizes_m": MODEL_SIZES_M,
        "geometry_audits": details,
        "total_audit_wall_seconds": time.perf_counter() - audit_started,
        "no_heat_or_flow_run": True,
    }
    (output / "resolution_audit.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
