#!/usr/bin/env python3
"""Generate the user-selected, geometry-only Ag-Ti3SiC2 candidate 02.

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
from matplotlib.lines import Line2D
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
SEEDS = [20261009]
FIELD_SIGMA_XYZ_VOXELS = (1.7, 2.7, 2.7)


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
    step = 2
    surface_structure = ndi.generate_binary_structure(3, 1)
    tic_surface = mask & ~ndi.binary_erosion(mask, structure=surface_structure, border_value=0)
    ag_mask = ~mask
    ag_surface = ag_mask & ~ndi.binary_erosion(ag_mask, structure=surface_structure, border_value=0)
    tic_ix, tic_iy, tic_iz = np.nonzero(tic_surface[::step, ::step, ::step])
    ag_ix, ag_iy, ag_iz = np.nonzero(ag_surface[::step, ::step, ::step])

    def scatter_phase_surface(ax, ix, iy, iz, color, size, alpha):
        ax.scatter((ix + 0.5) * SPACING_A[0] * step,
                   (iy + 0.5) * SPACING_A[1] * step,
                   (iz + 0.5) * SPACING_A[2] * step,
                   c=color, s=size, alpha=alpha, linewidths=0, depthshade=False)

    # Draw Ti first and the more opaque gray Ag points last so both
    # complementary phase boundaries remain legible in the static view.
    scatter_phase_surface(ax3, tic_ix, tic_iy, tic_iz, "#2563eb", 1.8, 0.68)
    scatter_phase_surface(ax3, ag_ix, ag_iy, ag_iz, "#8f97a3", 1.8, 0.62)
    ax3.set_xlim(0, LENGTHS_A[0])
    ax3.set_ylim(0, LENGTHS_A[1])
    ax3.set_zlim(0, LENGTHS_A[2])
    ax3.set_box_aspect(LENGTHS_A.copy())
    ax3.set_xlabel("X (Å)")
    ax3.set_ylabel("Y (Å)")
    ax3.set_zlabel("Z (Å)")
    ax3.view_init(elev=22, azim=35)
    ax3.set_title("3D two-phase connected morphology")

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

    # Standalone 3D view for easier inspection than the four-panel summary.
    fig3 = plt.figure(figsize=(10, 8), constrained_layout=True)
    ax3 = fig3.add_subplot(1, 1, 1, projection="3d")
    scatter_phase_surface(ax3, tic_ix, tic_iy, tic_iz, "#2563eb", 2.0, 0.68)
    scatter_phase_surface(ax3, ag_ix, ag_iy, ag_iz, "#8f97a3", 2.0, 0.62)
    ax3.set_xlim(0, LENGTHS_A[0])
    ax3.set_ylim(0, LENGTHS_A[1])
    ax3.set_zlim(0, LENGTHS_A[2])
    ax3.set_box_aspect(LENGTHS_A.copy())
    ax3.set_xlabel("X (Å)")
    ax3.set_ylabel("Y (Å)")
    ax3.set_zlabel("Z (Å)")
    ax3.view_init(elev=22, azim=35)
    ax3.set_title("Connected two-phase surface points")
    ax3.legend(handles=[
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#2563eb", markersize=7, label="Ti₃SiC₂"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#8f97a3", markersize=7, label="Ag"),
    ], loc="upper left")
    fig3.suptitle(f"Candidate {candidate_no:02d} | seed={seed} | blue: Ti₃SiC₂, gray: Ag")
    fig3.savefig(candidate_dir / "geometry_3d_view.png", dpi=180)
    plt.close(fig3)

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
    stats = generate_candidate(2, SEEDS[0])
    selected = {
        "selected_candidate": 2,
        "selection_status": "selected by user",
        "seed": SEEDS[0],
        "interpretation": "exploratory geometry only; not a calibrated phase-field prediction or atomistic model",
        "files": {
            "3d_view": "candidate_02/geometry_3d_view.png",
            "2d_reference_view": "candidate_02/morphology_reference_view_2d.png",
            "multi_panel_preview": "candidate_02/geometry_preview.png",
            "phase_masks_and_field": "candidate_02/phase_masks_and_field.npz",
            "statistics": "candidate_02/geometry_stats.json",
        },
        "volume_fraction_Ti3SiC2": stats["volume_fraction_actual_Ti3SiC2"],
        "mass_fraction_Ti3SiC2": stats["mass_fraction_actual_Ti3SiC2"],
        "Ti3SiC2_thickness_proxy_median_A": stats["Ti3SiC2_branch_thickness_proxy"]["p50_A"],
    }
    (OUT / "selected_candidate.json").write_text(
        json.dumps(selected, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "output": str(OUT),
        "candidate": 2,
        "seed": SEEDS[0],
        "Ti3SiC2_volume_fraction": stats["volume_fraction_actual_Ti3SiC2"],
        "Ti3SiC2_mass_fraction": stats["mass_fraction_actual_Ti3SiC2"],
        "Ti3SiC2_connectivity_26": stats["Ti3SiC2_connectivity_26"],
        "Ag_connectivity_26": stats["Ag_connectivity_26"],
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
