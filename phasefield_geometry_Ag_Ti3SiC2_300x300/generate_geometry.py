#!/usr/bin/env python3
"""Generate reproducible, geometry-only branched two-phase candidates.

This is an implicit-distance/branching construction, not a calibrated
Cahn-Hilliard calculation and not an atomistic structure generator.
"""
from __future__ import annotations

import json
import platform
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import scipy
from scipy import ndimage as ndi


OUT = Path(__file__).resolve().parent
LENGTHS_A = np.array([150.1, 300.0, 300.0], dtype=np.float64)
GRID = (76, 150, 150)
SPACING_A = LENGTHS_A / np.asarray(GRID)
RHO_TIC_A = 4.53
RHO_AG = 10.4976
WT_TIC = 0.4
TARGET_TIC_VOL = (WT_TIC / RHO_TIC_A) / (
    (WT_TIC / RHO_TIC_A) + ((1.0 - WT_TIC) / RHO_AG)
)
SEEDS = [20261008, 20261009, 20261010]
MAX_BRANCH_DEPTH = 5
CURVE_AMPLITUDE_A = 3.0
ROUGHNESS_RMS_A = 0.18


def clip_endpoint(start: np.ndarray, end: np.ndarray) -> np.ndarray:
    """Clip a segment endpoint at the first box boundary it crosses."""
    delta = end - start
    t_exit = 1.0
    for axis, length in enumerate(LENGTHS_A):
        if delta[axis] > 0 and end[axis] > length:
            t_exit = min(t_exit, (length - start[axis]) / delta[axis])
        elif delta[axis] < 0 and end[axis] < 0:
            t_exit = min(t_exit, (0.0 - start[axis]) / delta[axis])
    return start + max(0.0, t_exit) * delta


def branch_pair(direction: np.ndarray, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    direction = direction / np.linalg.norm(direction)
    side = rng.normal(size=3)
    side -= np.dot(side, direction) * direction
    norm = np.linalg.norm(side)
    if norm < 1e-10:
        side = np.cross(direction, np.array([1.0, 0.0, 0.0]))
        if np.linalg.norm(side) < 1e-10:
            side = np.cross(direction, np.array([0.0, 1.0, 0.0]))
        norm = np.linalg.norm(side)
    side /= norm
    angle = np.deg2rad(rng.uniform(30.0, 52.0))
    return (
        np.cos(angle) * direction + np.sin(angle) * side,
        np.cos(angle) * direction - np.sin(angle) * side,
    )


def curved_points(start: np.ndarray, end: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    vector = end - start
    length = np.linalg.norm(vector)
    if length < 1e-9:
        return np.asarray([start])
    direction = vector / length
    wiggle = rng.normal(size=3)
    wiggle -= np.dot(wiggle, direction) * direction
    norm = np.linalg.norm(wiggle)
    if norm > 1e-10:
        wiggle /= norm
    else:
        wiggle[:] = 0.0
    n = max(2, int(np.ceil(length / 0.8)) + 1)
    t = np.linspace(0.0, 1.0, n)
    amp = CURVE_AMPLITUDE_A * min(1.0, length / 30.0)
    bow = amp * np.sin(np.pi * t)
    points = start[None, :] + t[:, None] * vector[None, :] + bow[:, None] * wiggle[None, :]
    return np.clip(points, 0.0, LENGTHS_A[None, :])


def build_branch_network(seed: int) -> tuple[np.ndarray, dict]:
    rng = np.random.default_rng(seed)
    center = LENGTHS_A / 2.0
    centerline = np.zeros(GRID, dtype=bool)
    segments: list[tuple[np.ndarray, np.ndarray, int]] = []
    branchpoints = 0
    total_length = 0.0

    def draw_segment(start: np.ndarray, end: np.ndarray, depth: int) -> float:
        nonlocal total_length
        actual_end = clip_endpoint(start, end)
        actual_length = float(np.linalg.norm(actual_end - start))
        if actual_length < 1e-6:
            return 0.0
        segments.append((start.copy(), actual_end.copy(), depth))
        total_length += actual_length
        points = curved_points(start, actual_end, rng)
        indices = np.rint(points / SPACING_A[None, :] - 0.5).astype(np.int32)
        indices = np.clip(indices, 0, np.asarray(GRID) - 1)
        centerline[indices[:, 0], indices[:, 1], indices[:, 2]] = True
        return actual_length

    def add_tree(start: np.ndarray, direction: np.ndarray, length: float, depth: int) -> None:
        nonlocal branchpoints
        direction = direction / np.linalg.norm(direction)
        nominal_end = start + direction * length
        end = clip_endpoint(start, nominal_end)
        actual_length = draw_segment(start, end, depth)
        if actual_length < 4.0:
            return
        if depth >= MAX_BRANCH_DEPTH:
            return
        child_directions = branch_pair(direction, rng)
        branchpoints += 1
        for child_direction in child_directions:
            child_length = actual_length * rng.uniform(0.68, 0.80)
            add_tree(end, child_direction, child_length, depth + 1)

    # Six outward trunks meet at the center and span every axis. Branch trees
    # attach at several interior points, so they grow into the box rather than
    # immediately being clipped when rooted at a face.
    for axis in range(3):
        for sign in (-1.0, 1.0):
            direction = np.zeros(3)
            direction[axis] = sign
            arm_length = LENGTHS_A[axis] / 2.0
            endpoint = center + direction * arm_length
            draw_segment(center, endpoint, -1)
            for fraction in (0.22, 0.42, 0.62, 0.80):
                anchor = center + direction * arm_length * fraction
                branchpoints += 1
                branch_directions = branch_pair(direction, rng)
                remaining = arm_length * (1.0 - fraction)
                for child_direction in branch_directions:
                    first_length = remaining * rng.uniform(0.58, 0.82)
                    add_tree(anchor, child_direction, first_length, 0)

    meta = {
        "seed": seed,
        "algorithm": "six spanning trunks with interior binary branch trees, then signed-distance threshold",
        "branch_depth": MAX_BRANCH_DEPTH,
        "trunk_count": 6,
        "branch_angle_deg_range": [30.0, 52.0],
        "trunk_anchor_fractions": [0.22, 0.42, 0.62, 0.80],
        "recursive_child_length_ratio_range": [0.68, 0.80],
        "initial_branch_length_fraction_of_remaining_trunk_range": [0.58, 0.82],
        "centerline_sampling_max_step_A": 0.8,
        "generated_segment_count": len(segments),
        "generated_branchpoint_count": branchpoints,
        "centerline_length_A": total_length,
        "curve_amplitude_A": CURVE_AMPLITUDE_A,
        "boundary_handling": "nonperiodic box; branch segments clipped at box faces",
    }
    return centerline, meta


def exact_volume_field(centerline: np.ndarray, seed: int) -> tuple[np.ndarray, np.ndarray, float]:
    # EDT gives physical distance in Angstroms using the anisotropic voxel size.
    distance = ndi.distance_transform_edt(~centerline, sampling=SPACING_A).astype(np.float32)
    rng = np.random.default_rng(seed ^ 0x5A17C3)
    roughness = rng.standard_normal(GRID, dtype=np.float32)
    roughness = ndi.gaussian_filter(roughness, sigma=(1.2, 1.6, 1.6), mode="reflect")
    roughness /= max(float(roughness.std()), 1e-8)
    score = distance - np.float32(ROUGHNESS_RMS_A) * roughness

    n_total = score.size
    n_tic = int(round(TARGET_TIC_VOL * n_total))
    flat = score.ravel()
    order = np.argpartition(flat, n_tic - 1)
    chosen = order[:n_tic]
    threshold = float(flat[chosen].max())
    mask_flat = np.zeros(n_total, dtype=bool)
    mask_flat[chosen] = True
    field = (np.float32(threshold) - score).astype(np.float32, copy=False)

    # Make any threshold ties explicit so that field sign and exact target count agree.
    ties = np.flatnonzero(flat == threshold)
    selected_ties = np.intersect1d(chosen, ties, assume_unique=False)
    unselected_ties = np.setdiff1d(ties, selected_ties, assume_unique=False)
    epsilon = np.float32(1e-5)
    field.ravel()[selected_ties] = epsilon
    field.ravel()[unselected_ties] = -epsilon
    mask = mask_flat.reshape(GRID)
    assert int(mask.sum()) == n_tic
    assert np.array_equal(mask, field > 0)
    return field, mask, threshold


def component_stats(mask: np.ndarray, connectivity: int) -> dict:
    structure = ndi.generate_binary_structure(3, connectivity)
    labels, count = ndi.label(mask, structure=structure)
    if count == 0:
        return {"component_count": 0, "largest_component_voxels": 0, "largest_component_fraction": 0.0,
                "largest_component_spans": {"X": False, "Y": False, "Z": False}}
    sizes = np.bincount(labels.ravel())[1:]
    largest_id = int(np.argmax(sizes)) + 1
    largest_size = int(sizes[largest_id - 1])
    box_spans = {}
    for axis, name in enumerate(("X", "Y", "Z")):
        low = np.unique(labels.take(0, axis=axis))
        high = np.unique(labels.take(labels.shape[axis] - 1, axis=axis))
        box_spans[name] = bool(np.any(np.intersect1d(low[low != 0], high[high != 0]) == largest_id))
    return {
        "component_count": int(count),
        "largest_component_voxels": largest_size,
        "largest_component_fraction": largest_size / int(mask.sum()),
        "largest_component_spans": box_spans,
    }


def thickness_stats(mask: np.ndarray) -> dict:
    padded = np.pad(mask, 1, mode="constant", constant_values=False)
    padded_radii = ndi.distance_transform_edt(padded, sampling=SPACING_A)
    radii = padded_radii[1:-1, 1:-1, 1:-1]
    maxima = ndi.maximum_filter(radii, size=3, mode="constant")
    ridge = mask & (radii >= maxima - 1e-9) & (radii > 0)
    values = 2.0 * radii[ridge]
    if values.size == 0:
        return {"ridge_samples": 0}
    p = np.percentile(values, [10, 50, 90, 95])
    return {
        "definition": "2x interior EDT at 3x3x3 local-radius maxima; geometric thickness proxy, Angstrom",
        "ridge_samples": int(values.size),
        "p10_A": float(p[0]),
        "p50_A": float(p[1]),
        "p90_A": float(p[2]),
        "p95_A": float(p[3]),
    }


def render_candidate(mask: np.ndarray, candidate_dir: Path, candidate_no: int, seed: int) -> None:
    colors = ListedColormap(["#1d7d91", "#ef8a3b"])  # Ag, Ti3SiC2
    fig = plt.figure(figsize=(15, 10), constrained_layout=True)
    ax3 = fig.add_subplot(2, 2, 1, projection="3d")
    surface = mask & ~ndi.binary_erosion(
        mask, structure=ndi.generate_binary_structure(3, 1), border_value=0
    )
    step = 2
    sampled = surface[::step, ::step, ::step]
    ix, iy, iz = np.nonzero(sampled)
    ax3.scatter((ix + 0.5) * SPACING_A[0] * step,
                (iy + 0.5) * SPACING_A[1] * step,
                (iz + 0.5) * SPACING_A[2] * step,
                c="#e88436", s=1.8, alpha=0.7, linewidths=0, depthshade=False)
    ax3.set_xlim(0, LENGTHS_A[0])
    ax3.set_ylim(0, LENGTHS_A[1])
    ax3.set_zlim(0, LENGTHS_A[2])
    ax3.set_box_aspect(LENGTHS_A.copy())
    ax3.set_xlabel("X (Å)")
    ax3.set_ylabel("Y (Å)")
    ax3.set_zlabel("Z (Å)")
    ax3.view_init(elev=22, azim=35)
    ax3.set_title("3D Ti₃SiC₂ branched network")

    mx, my, mz = np.asarray(GRID) // 2
    panels = [
        (fig.add_subplot(2, 2, 2), mask[mx, :, :], (LENGTHS_A[2], LENGTHS_A[1]), "Y (Å)", "Z (Å)", f"YZ slice, X={((mx + .5) * SPACING_A[0]):.1f} Å"),
        (fig.add_subplot(2, 2, 3), mask[:, my, :], (LENGTHS_A[2], LENGTHS_A[0]), "X (Å)", "Z (Å)", f"XZ slice, Y={((my + .5) * SPACING_A[1]):.1f} Å"),
        (fig.add_subplot(2, 2, 4), mask[:, :, mz], (LENGTHS_A[1], LENGTHS_A[0]), "X (Å)", "Y (Å)", f"XY slice, Z={((mz + .5) * SPACING_A[2]):.1f} Å"),
    ]
    for ax, image, extent_sizes, xlabel, ylabel, title in panels:
        ax.imshow(image.T.astype(np.uint8), origin="lower", interpolation="nearest",
                  cmap=colors, vmin=0, vmax=1, aspect="equal",
                  extent=(0, extent_sizes[0], 0, extent_sizes[1]))
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.ticklabel_format(style="plain", useOffset=False)
    fig.suptitle(f"Candidate {candidate_no:02d} | seed={seed} | orange: Ti₃SiC₂, blue: Ag")
    fig.savefig(candidate_dir / "geometry_preview.png", dpi=180)
    plt.close(fig)


def generate_candidate(candidate_no: int, seed: int) -> dict:
    candidate_dir = OUT / f"candidate_{candidate_no:02d}"
    candidate_dir.mkdir(parents=True, exist_ok=True)
    centerline, network_meta = build_branch_network(seed)
    field, tic_mask, threshold = exact_volume_field(centerline, seed)
    ag_mask = ~tic_mask
    n_voxels = tic_mask.size
    actual_phi = float(tic_mask.mean())
    actual_wt = (actual_phi * RHO_TIC_A) / (
        actual_phi * RHO_TIC_A + (1.0 - actual_phi) * RHO_AG
    )
    stats = {
        "candidate": candidate_no,
        **network_meta,
        "box_A": {"X": float(LENGTHS_A[0]), "Y": float(LENGTHS_A[1]), "Z": float(LENGTHS_A[2])},
        "grid_shape_XYZ": list(GRID),
        "voxel_spacing_A_XYZ": [float(x) for x in SPACING_A],
        "voxel_count": int(n_voxels),
        "density_inputs_g_cm3": {"Ti3SiC2": RHO_TIC_A, "Ag": RHO_AG},
        "mass_fraction_target_Ti3SiC2": WT_TIC,
        "volume_fraction_target_Ti3SiC2": TARGET_TIC_VOL,
        "volume_fraction_actual_Ti3SiC2": actual_phi,
        "volume_fraction_actual_Ag": float(ag_mask.mean()),
        "mass_fraction_actual_Ti3SiC2": actual_wt,
        "mass_fraction_actual_Ag": 1.0 - actual_wt,
        "porosity": 0.0,
        "field_threshold_A": threshold,
        "boundary_roughness_rms_A": ROUGHNESS_RMS_A,
        "software_versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "Ti3SiC2_connectivity_6": component_stats(tic_mask, 1),
        "Ti3SiC2_connectivity_26": component_stats(tic_mask, 3),
        "Ag_connectivity_6": component_stats(ag_mask, 1),
        "Ag_connectivity_26": component_stats(ag_mask, 3),
        "Ti3SiC2_branch_thickness_proxy": thickness_stats(tic_mask),
        "Ag_channel_thickness_proxy": thickness_stats(ag_mask),
        "interpretation": "exploratory geometry only; not a calibrated phase-field or atomistic model",
    }
    np.savez_compressed(
        candidate_dir / "phase_masks_and_field.npz",
        signed_geometry_field_A=field,
        Ti3SiC2_mask=tic_mask.astype(np.uint8),
        Ag_mask=ag_mask.astype(np.uint8),
        centerline_mask=centerline.astype(np.uint8),
        voxel_spacing_A=SPACING_A,
        box_lengths_A=LENGTHS_A,
    )
    (candidate_dir / "geometry_stats.json").write_text(
        json.dumps(stats, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    render_candidate(tic_mask, candidate_dir, candidate_no, seed)
    return stats


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    stats = [generate_candidate(i + 1, seed) for i, seed in enumerate(SEEDS)]
    summary = {
        "title": "Ag-Ti3SiC2 branched geometry candidates",
        "classification": "exploratory geometry design; uncalibrated; not a physical spinodal prediction",
        "seeds": SEEDS,
        "target_ti3sic2_vol_fraction": TARGET_TIC_VOL,
        "candidates": stats,
    }
    (OUT / "candidate_comparison.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "output": str(OUT),
        "candidates": [
            {
                "candidate": s["candidate"],
                "seed": s["seed"],
                "Ti3SiC2_vol_fraction": s["volume_fraction_actual_Ti3SiC2"],
                "Ti3SiC2_mass_fraction": s["mass_fraction_actual_Ti3SiC2"],
                "Ti3SiC2_components_26": s["Ti3SiC2_connectivity_26"]["component_count"],
                "Ag_components_26": s["Ag_connectivity_26"]["component_count"],
                "Ti3SiC2_spans": s["Ti3SiC2_connectivity_26"]["largest_component_spans"],
                "Ag_spans": s["Ag_connectivity_26"]["largest_component_spans"],
                "thickness_proxy_A": s["Ti3SiC2_branch_thickness_proxy"],
            } for s in stats
        ]
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
