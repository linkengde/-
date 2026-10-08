#!/usr/bin/env python3
"""Generate reproducible, geometry-only two-phase morphology candidates.

The candidates come from thresholded correlated random fields. This is a
geometric design method, not a calibrated Cahn-Hilliard calculation or an
atomistic structure generator.
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
FIELD_SIGMA_XYZ_VOXELS = (3.0, 4.0, 4.0)


def exact_volume_field(seed: int) -> tuple[np.ndarray, np.ndarray, float]:
    """Threshold a smooth, seeded 3D random field to the target Ti fraction."""
    rng = np.random.default_rng(seed)
    noise = rng.standard_normal(GRID, dtype=np.float32)
    smooth = ndi.gaussian_filter(noise, sigma=FIELD_SIGMA_XYZ_VOXELS, mode="reflect")
    smooth /= max(float(smooth.std()), 1e-8)
    score = -smooth  # Lowest scores are assigned to the Ti3SiC2 phase.

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
    colors = ListedColormap(["#d2d4d8", "#2563eb"])  # Ag, Ti3SiC2
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
                c="#2563eb", s=1.8, alpha=0.7, linewidths=0, depthshade=False)
    ax3.set_xlim(0, LENGTHS_A[0])
    ax3.set_ylim(0, LENGTHS_A[1])
    ax3.set_zlim(0, LENGTHS_A[2])
    ax3.set_box_aspect(LENGTHS_A.copy())
    ax3.set_xlabel("X (Å)")
    ax3.set_ylabel("Y (Å)")
    ax3.set_zlabel("Z (Å)")
    ax3.view_init(elev=22, azim=35)
    ax3.set_title("3D Ti₃SiC₂ connected morphology")

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
    fig.suptitle(f"Candidate {candidate_no:02d} | seed={seed} | blue: Ti₃SiC₂, gray: Ag")
    fig.savefig(candidate_dir / "geometry_preview.png", dpi=180)
    plt.close(fig)

    # A clean square view for direct comparison with the user's reference.
    # The YZ slice is 150 x 150 voxels (2 Å per pixel); nearest-neighbor
    # enlargement gives a 300 x 300 image without smoothing phase boundaries.
    mx = np.asarray(GRID)[0] // 2
    slice_yz = np.repeat(np.repeat(mask[mx, :, :].T, 2, axis=0), 2, axis=1)
    rgb = np.empty((*slice_yz.shape, 3), dtype=np.uint8)
    rgb[~slice_yz] = (210, 212, 216)  # Ag, light gray
    rgb[slice_yz] = (37, 99, 235)     # Ti3SiC2, blue
    plt.imsave(candidate_dir / "morphology_reference_view_2d.png", rgb)


def generate_candidate(candidate_no: int, seed: int) -> dict:
    candidate_dir = OUT / f"candidate_{candidate_no:02d}"
    candidate_dir.mkdir(parents=True, exist_ok=True)
    field, tic_mask, threshold = exact_volume_field(seed)
    ag_mask = ~tic_mask
    n_voxels = tic_mask.size
    actual_phi = float(tic_mask.mean())
    actual_wt = (actual_phi * RHO_TIC_A) / (
        actual_phi * RHO_TIC_A + (1.0 - actual_phi) * RHO_AG
    )
    stats = {
        "candidate": candidate_no,
        "seed": seed,
        "algorithm": "thresholded 3D correlated random field; geometry design only",
        "field_correlation_sigma_voxels_XYZ": list(FIELD_SIGMA_XYZ_VOXELS),
        "boundary_handling": "nonperiodic box; Gaussian smoothing uses reflection at box faces",
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
        "geometry_score_threshold": threshold,
        "preview_color_mapping": {"Ti3SiC2": "blue #2563eb", "Ag": "light gray #d2d4d8"},
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
        signed_geometry_field=field,
        Ti3SiC2_mask=tic_mask.astype(np.uint8),
        Ag_mask=ag_mask.astype(np.uint8),
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
        "title": "Ag-Ti3SiC2 morphology geometry candidates",
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
