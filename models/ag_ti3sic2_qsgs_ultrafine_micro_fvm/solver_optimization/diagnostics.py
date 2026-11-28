"""Ag-only liquid metrics, threshold sensitivity, and conservative-map checks."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy import ndimage


def phase_metrics_from_arrays(ag, liquid, spacing_m, liquid_threshold=0.5, ag_threshold=0.5,
                              stored_energy_change_J=0.0, source_energy_J=0.0,
                              boundary_energy_out_J=0.0):
    ag = np.asarray(ag)
    liquid = np.asarray(liquid)
    ag_present = ag > 1.0e-14
    if np.any(ag_present):
        maximum_ag_liquid = float(np.max(liquid[ag_present]))
    else:
        maximum_ag_liquid = 0.0
    weighted_fraction = ag * liquid
    voxel_volume = float(np.prod(spacing_m))
    ag_volume = float(np.sum(ag) * voxel_volume)
    weighted_liquid_volume = float(np.sum(weighted_fraction) * voxel_volume)
    stored_energy_change = float(stored_energy_change_J)
    energy_residual = (stored_energy_change - source_energy_J + boundary_energy_out_J)
    active = (ag >= ag_threshold) & (liquid >= liquid_threshold)
    labels, count = ndimage.label(active, structure=ndimage.generate_binary_structure(3, 3))
    if count:
        counts = np.bincount(labels.ravel())[1:]
        max_label = int(np.argmax(counts)) + 1
        max_cells = int(counts.max())
        max_weighted_volume = float(np.sum(weighted_fraction[labels == max_label]) * voxel_volume)
    else:
        max_cells = 0
        max_weighted_volume = 0.0
    return {
        "Ag_max_liquid_fraction": maximum_ag_liquid,
        "Ag_bulk_liquid_fraction": weighted_liquid_volume / max(ag_volume, 1.0e-30),
        "Ag_liquid_weighted_volume_m3": weighted_liquid_volume,
        "stored_energy_change_J": stored_energy_change,
        "source_energy_J": float(source_energy_J),
        "boundary_energy_out_J": float(boundary_energy_out_J),
        "energy_residual_J": energy_residual,
        "liquid_threshold": float(liquid_threshold),
        "Ag_fraction_threshold": float(ag_threshold),
        "thresholded_geometric_volume_m3": float(np.count_nonzero(active) * voxel_volume),
        "thresholded_connected_component_count_26n": int(count),
        "thresholded_max_component_cells_26n": max_cells,
        "thresholded_max_component_geometric_volume_m3_26n": max_cells * voxel_volume,
        "thresholded_max_component_Ag_liquid_weighted_volume_m3_26n": max_weighted_volume,
        "mixed_cell_fraction": float(np.mean((ag > 1.0e-12) & (ag < 1.0 - 1.0e-12))),
        "connectivity_interpretation": "thresholded mixed control volumes; not proof of a physical Ag channel",
    }


def phase_metrics(sim, liquid_threshold=0.5, ag_threshold=0.5):
    stored = float(np.sum(sim.enthalpy_density_J_m3) * np.prod(sim.spacing_m)
                   - sim.initial_stored_enthalpy_J)
    return phase_metrics_from_arrays(
        sim.ag_fraction, sim.materials.liquid_fraction(sim.temperature_K), sim.spacing_m,
        liquid_threshold, ag_threshold, stored, sim.source_energy_J, sim.boundary_energy_out_J,
    )


def threshold_sensitivity(sim, ag_thresholds=(0.1, 0.25, 0.5, 0.75, 0.9),
                          liquid_thresholds=(0.25, 0.5, 0.75)):
    rows = []
    for ag_threshold in ag_thresholds:
        for liquid_threshold in liquid_thresholds:
            rows.append(phase_metrics(sim, liquid_threshold, ag_threshold))
    return rows


def threshold_sensitivity_from_arrays(ag, liquid, spacing_m,
                                      ag_thresholds=(0.1, 0.25, 0.5, 0.75, 0.9),
                                      liquid_thresholds=(0.25, 0.5, 0.75)):
    rows = []
    for ag_threshold in ag_thresholds:
        for liquid_threshold in liquid_thresholds:
            rows.append(phase_metrics_from_arrays(
                ag, liquid, spacing_m, liquid_threshold, ag_threshold,
            ))
    return rows


def write_csv(path, rows):
    rows = list(rows)
    if not rows:
        return
    fieldnames = list(dict.fromkeys(key for row in rows for key in row))
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def mapping_diagnostics(seed_ag, n_target, edge_m, source_rate, radius_fraction):
    """Report mixed-cell share, phase connectivity, and source-integral mapping error."""
    from run_3d_phase_change import conservative_remap_3d

    seed_ag = np.asarray(seed_ag, dtype=np.float64)
    mapped = conservative_remap_3d(seed_ag, int(n_target))
    edge = float(edge_m)
    dx = edge / n_target
    native_dx = edge / seed_ag.shape[0]
    def source_mask(n, spacing):
        coordinates = (np.arange(n, dtype=np.float64) + 0.5 - 0.5 * n) * spacing
        radius = radius_fraction * edge
        return (coordinates[:, None, None] ** 2 + coordinates[None, :, None] ** 2
                + coordinates[None, None, :] ** 2) <= radius * radius
    native_integral = float(source_rate * np.sum(seed_ag * source_mask(seed_ag.shape[0], native_dx)) * native_dx**3)
    mapped_integral = float(source_rate * np.sum(mapped * source_mask(n_target, dx)) * dx**3)
    rows = []
    for phase, volume_fraction in (("Ag", mapped), ("Ti3SiC2", 1.0 - mapped)):
        binary = volume_fraction >= 0.5
        for connectivity_rank, name in ((1, "6-neighbour"), (2, "18-neighbour"), (3, "26-neighbour")):
            labels, count = ndimage.label(binary, structure=ndimage.generate_binary_structure(3, connectivity_rank))
            spans = {}
            for axis, axis_name in enumerate("xyz"):
                low_selector = [slice(None)] * 3
                high_selector = [slice(None)] * 3
                low_selector[axis] = 0
                high_selector[axis] = -1
                low_ids = np.unique(labels[tuple(low_selector)])
                high_ids = np.unique(labels[tuple(high_selector)])
                spans[axis_name] = bool(np.intersect1d(low_ids[low_ids > 0], high_ids[high_ids > 0]).size)
            rows.append({
                "resolution": int(n_target), "phase": phase, "connectivity": name,
                "volume_fraction_mean": float(volume_fraction.mean()),
                "mixed_cell_fraction": float(np.mean((mapped > 1e-12) & (mapped < 1.0 - 1e-12))),
                "component_count": int(count),
                "spans_x": spans["x"], "spans_y": spans["y"], "spans_z": spans["z"],
                "native_source_rate_integral_W": native_integral,
                "mapped_source_rate_integral_W": mapped_integral,
                "source_integral_relative_error": abs(mapped_integral - native_integral) / max(abs(native_integral), 1e-30),
                "mapping_note": "volume-overlap remap; binary connectivity uses phase fraction >=0.5",
            })
    return rows
