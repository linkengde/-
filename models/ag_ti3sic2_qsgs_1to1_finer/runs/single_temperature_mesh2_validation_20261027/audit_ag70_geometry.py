#!/usr/bin/env python3
"""Compute grid-explicit morphology metrics for the existing Ag70 QSGS mask."""

from __future__ import annotations

import json
import hashlib
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
GEOMETRY = HERE / "qsgs_morphology_ag70vol.npz"


def periodic_chords(mask: np.ndarray, axis: int, periodic: bool, spacing: float) -> np.ndarray:
    moved = np.moveaxis(mask, axis, -1).reshape(-1, mask.shape[axis])
    lengths: list[int] = []
    for row in moved:
        if not row.any():
            continue
        if row.all():
            lengths.append(row.size)
            continue
        if periodic:
            start = int(np.flatnonzero(row != np.roll(row, 1))[0])
            ordered = np.roll(row, -start)
        else:
            ordered = row
        padded = np.concatenate(([False], ordered, [False]))
        transitions = np.flatnonzero(padded[1:] != padded[:-1])
        run_lengths = transitions[1::2] - transitions[::2]
        lengths.extend(int(x) for x in run_lengths)
    return np.asarray(lengths, dtype=float) * spacing


def summary(values: np.ndarray) -> dict[str, float | int | None]:
    if not values.size:
        return {"count": 0, "p10_A": None, "median_A": None, "p90_A": None, "max_A": None}
    return {
        "count": int(values.size),
        "p10_A": float(np.percentile(values, 10)),
        "median_A": float(np.median(values)),
        "p90_A": float(np.percentile(values, 90)),
        "max_A": float(values.max()),
    }


def main() -> None:
    with np.load(GEOMETRY, allow_pickle=False) as data:
        skeleton = np.asarray(data["skeleton"], dtype=bool)
        size_A = np.asarray(data["size_A"], dtype=float)
    spacing = size_A / np.asarray(skeleton.shape)
    ag = ~skeleton

    mismatch_counts: dict[str, int] = {}
    interface_area_A2: dict[str, float] = {}
    for name, axis in zip(("x", "y", "z"), range(3), strict=True):
        neighbour = np.roll(skeleton, -1, axis=axis) if axis < 2 else skeleton.take(range(1, skeleton.shape[2]), axis=2)
        phase = skeleton if axis < 2 else skeleton.take(range(0, skeleton.shape[2] - 1), axis=2)
        mismatch = int(np.count_nonzero(phase != neighbour))
        mismatch_counts[name] = mismatch
        other_axes = [i for i in range(3) if i != axis]
        face_area = float(np.prod(spacing[other_axes]))
        interface_area_A2[name] = mismatch * face_area
    interface_area_A2_total = sum(interface_area_A2.values())
    volume_A3 = float(np.prod(size_A))
    phase_volumes = {"Ag": float(ag.mean()), "Ti3SiC2": float(skeleton.mean())}
    phase_feature_scale = {
        phase: (2.0 * fraction / (interface_area_A2_total / volume_A3))
        for phase, fraction in phase_volumes.items()
    }

    chords: dict[str, dict[str, dict[str, float | int | None]]] = {}
    for phase_name, mask in (("Ag", ag), ("Ti3SiC2", skeleton)):
        chords[phase_name] = {}
        for name, axis in zip(("x", "y", "z"), range(3), strict=True):
            lengths = periodic_chords(mask, axis, periodic=(axis < 2), spacing=float(spacing[axis]))
            chords[phase_name][name] = summary(lengths)

    result = {
        "status": "RASTER_RESOLVED_MORPHOLOGY_AUDIT",
        "geometry": GEOMETRY.name,
        "geometry_sha256": hashlib.sha256(GEOMETRY.read_bytes()).hexdigest(),
        "shape_xyz": list(skeleton.shape),
        "voxel_spacing_A": spacing.tolist(),
        "domain_size_A": size_A.tolist(),
        "phase_volume_fraction": phase_volumes,
        "connectivity_boundary_convention": "periodic x/y seams, open z",
        "interfacial_faces_by_normal": mismatch_counts,
        "interface_area_by_normal_A2": interface_area_A2,
        "total_Ag_Ti3SiC2_interface_area_A2": interface_area_A2_total,
        "interface_area_density_A_minus1": interface_area_A2_total / volume_A3,
        "phase_scale_2phi_over_Sv_A": phase_feature_scale,
        "phase_lineal_chord_lengths_A": chords,
        "interpretation_limits": [
            "Face counts and chord lengths are exact for this binary raster and stated periodic/open boundaries.",
            "Lineal chord distributions are directional raster metrics, not minimum hydraulic neck thicknesses.",
            "2*phase_volume_fraction/interface_area_density is an interface-based characteristic scale, not a channel-width measurement.",
            "Features below the approximately 3 A native voxel size and independent-realization variability are unresolved.",
        ],
    }
    (HERE / "ag70_geometry_metrics.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
