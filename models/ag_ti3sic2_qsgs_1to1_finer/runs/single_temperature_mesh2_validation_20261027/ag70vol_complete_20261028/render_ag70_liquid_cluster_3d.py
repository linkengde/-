#!/usr/bin/env python3
"""Render the Ag70 1x melt cluster locally and inside the full 3D sample."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


RUN_DIR = Path(__file__).resolve().parent.parent
MODEL_DIR = RUN_DIR.parents[1]
GEOMETRY = RUN_DIR / "qsgs_morphology_ag70vol.npz"
CHECKPOINT = (
    RUN_DIR
    / "ag70vol_resume_1x_from1ps_20261028"
    / "checkpoints"
    / "caseA_ag70_1x_final_4p52ps.npz"
)
OUTPUT = Path(__file__).resolve().parent / "ag70_local_liquid_cluster_3d_4p52ps.png"
OUTPUT_WHOLE = Path(__file__).resolve().parent / "ag70_whole_sample_with_liquid_pool_3d_4p52ps.png"


def draw_box(ax, bounds: tuple[tuple[float, float], ...], color: str, linewidth: float) -> None:
    (x0, x1), (y0, y1), (z0, z1) = bounds
    corners = [
        (x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
        (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1),
    ]
    edges = (
        (0, 1), (1, 2), (2, 3), (3, 0),
        (4, 5), (5, 6), (6, 7), (7, 4),
        (0, 4), (1, 5), (2, 6), (3, 7),
    )
    for start, end in edges:
        a, b = corners[start], corners[end]
        ax.plot([a[0], b[0]], [a[1], b[1]], [a[2], b[2]], color=color, linewidth=linewidth)


def main() -> None:
    if not GEOMETRY.is_file() or not CHECKPOINT.is_file():
        raise FileNotFoundError(f"Required QSGS geometry or saved checkpoint is missing: {CHECKPOINT}")

    phase_dir = MODEL_DIR / "phase_field_qsgs"
    sys.path.insert(0, str(phase_dir))
    from qsgs_liquid_connectivity import liquid_cluster_labels

    with np.load(GEOMETRY, allow_pickle=False) as archive:
        scaffold = np.asarray(archive["skeleton"], dtype=bool)
        domain_size_A = np.asarray(archive["size_A"], dtype=float)
    with np.load(CHECKPOINT, allow_pickle=False) as archive:
        temperature = np.asarray(archive["temperature_K"], dtype=float)
        liquid_fraction = np.asarray(archive["ag_liquid_fraction"], dtype=float)
        metadata = json.loads(str(archive["metadata_json"].item()))

    ag = ~scaffold
    voxel_A = domain_size_A / np.asarray(ag.shape, dtype=float)
    liquid_mask = ag & (liquid_fraction >= 0.5)
    labels, count = liquid_cluster_labels(liquid_mask)
    sizes = np.bincount(labels.ravel())[1:]
    if count == 0 or sizes.size == 0:
        raise ValueError("The selected saved state contains no connected liquid Ag cluster.")
    largest_label = int(np.argmax(sizes) + 1)
    largest_mask = labels == largest_label
    coords = np.where(largest_mask)
    lower = np.maximum(np.asarray([axis.min() for axis in coords]) - 4, 0)
    upper = np.minimum(np.asarray([axis.max() + 5 for axis in coords]), np.asarray(ag.shape))
    slices = tuple(slice(int(lo), int(hi)) for lo, hi in zip(lower, upper, strict=True))

    roi_ag = ag[slices]
    roi_liquid = largest_mask[slices]
    roi_partial = ag[slices] & (liquid_fraction[slices] > 1.0e-12) & (liquid_fraction[slices] < 0.5)
    axes_edges = [
        np.arange(int(lower[i]), int(upper[i]) + 1, dtype=float) * voxel_A[i]
        for i in range(3)
    ]
    edge_grid = np.meshgrid(*axes_edges, indexing="ij")
    centers = [
        (np.arange(int(lower[i]), int(upper[i])) + 0.5) * voxel_A[i]
        for i in range(3)
    ]

    fig = plt.figure(figsize=(15.0, 7.6))
    fig.subplots_adjust(left=0.035, right=0.97, top=0.80, bottom=0.14, wspace=0.02)
    ax_local = fig.add_subplot(1, 2, 1, projection="3d")
    context = roi_ag & ~roi_liquid & ~roi_partial
    xyz = np.where(context)
    ax_local.scatter(
        centers[0][xyz[0]], centers[1][xyz[1]], centers[2][xyz[2]],
        s=7, c="#7792a5", alpha=0.18, depthshade=False, linewidths=0,
    )
    if np.any(roi_partial):
        ax_local.voxels(
            *edge_grid, roi_partial,
            facecolors=(1.0, 0.78, 0.26, 0.85),
            edgecolors=(0.48, 0.30, 0.08, 0.8),
            linewidth=0.35,
        )
    ax_local.voxels(
        *edge_grid, roi_liquid,
        facecolors="#f28e2b",
        edgecolors="#7f3b08",
        linewidth=0.42,
    )
    ax_local.set_title(
        f"Largest connected liquid Ag cluster\n"
        f"{int(sizes.max())} voxels · {sizes.max() * np.prod(voxel_A):.0f} Å³"
    )
    ax_local.set_xlabel("x (Å)")
    ax_local.set_ylabel("y (Å)")
    ax_local.set_zlabel("z (Å)")
    ax_local.set_box_aspect((upper - lower) * voxel_A)
    ax_local.view_init(elev=27, azim=-54)
    ax_local.set_proj_type("ortho")

    ax_global = fig.add_subplot(1, 2, 2, projection="3d")
    draw_box(
        ax_global,
        ((0.0, domain_size_A[0]), (0.0, domain_size_A[1]), (0.0, domain_size_A[2])),
        "#8796a5",
        0.75,
    )
    cluster_bounds = tuple(
        (float(lower[i] * voxel_A[i]), float(upper[i] * voxel_A[i]))
        for i in range(3)
    )
    draw_box(ax_global, cluster_bounds, "#f28e2b", 2.2)
    centroid_A = np.asarray([(coord.mean() + 0.5) * voxel_A[i] for i, coord in enumerate(coords)])
    ax_global.scatter(*centroid_A, color="#e4572e", s=45, marker="o", depthshade=False)
    source = np.asarray([domain_size_A[0] / 2, domain_size_A[1] / 2, domain_size_A[2]])
    ax_global.scatter(*source, color="#4c78a8", s=70, marker="*", depthshade=False)
    ax_global.text(*centroid_A, "  liquid cluster", color="#9a4b00", fontsize=8)
    ax_global.set_xlim(0, domain_size_A[0])
    ax_global.set_ylim(0, domain_size_A[1])
    ax_global.set_zlim(0, domain_size_A[2])
    ax_global.set_box_aspect(domain_size_A.copy())
    ax_global.set_proj_type("ortho")
    ax_global.view_init(elev=23, azim=-53)
    ax_global.set_xticks((0, domain_size_A[0] / 2, domain_size_A[0]))
    ax_global.set_yticks((0, domain_size_A[1] / 2, domain_size_A[1]))
    ax_global.set_zticks((0, domain_size_A[2] / 2, domain_size_A[2]))
    ax_global.set_xticklabels([])
    ax_global.set_yticklabels([])
    ax_global.set_zticklabels([])
    ax_global.tick_params(labelsize=7, pad=0)
    ax_global.set_xlabel("x (Å)", labelpad=3)
    ax_global.set_ylabel("y (Å)", labelpad=3)
    ax_global.set_zlabel("z (Å)", labelpad=2)
    ax_global.set_title("Location in the full 70:30 vol% sample")
    ax_global.legend(
        handles=[
            Line2D([0], [0], marker="o", color="w", markerfacecolor="#e4572e", label="liquid-cluster centroid"),
            Line2D([0], [0], marker="*", color="w", markerfacecolor="#4c78a8", label="Gaussian source at top"),
        ],
        loc="upper left",
        fontsize=8,
    )

    time_ps = float(metadata["time_s"]) * 1.0e12
    peak_K = float(temperature.max())
    liquid_volfrac = float(np.sum(ag * liquid_fraction) / np.sum(ag))
    fig.suptitle(
        f"Ag 70:30 vol% Case A · 1× native geometry · t={time_ps:.2f} ps · "
        f"peak T={peak_K:.1f} K · Ag liquid fraction={liquid_volfrac:.3e}\n"
        "Saved heating state; orange voxels have liquid fraction ≥ 0.5"
    )
    fig.text(
        0.5, 0.045,
        f"Voxel spacing: {voxel_A[0]:.2f} × {voxel_A[1]:.2f} × {voxel_A[2]:.2f} Å; "
        f"sample size: {domain_size_A[0]:.1f} × {domain_size_A[1]:.1f} × {domain_size_A[2]:.1f} Å; "
        "grey dots: nearby Ag voxels; orange: liquid fraction ≥ 0.5; yellow: partial liquid; "
        "right panel marks the local crop.",
        ha="center", fontsize=9,
    )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=220, facecolor="white", bbox_inches="tight")
    print(f"saved={OUTPUT}")

    whole = plt.figure(figsize=(16.0, 8.4))
    whole.subplots_adjust(left=0.035, right=0.98, top=0.83, bottom=0.13, wspace=0.06)
    ax_sample = whole.add_subplot(1, 2, 1, projection="3d")
    sample_stride = 2
    ag_sample = ag[::sample_stride, ::sample_stride, ::sample_stride]
    sc_sample = scaffold[::sample_stride, ::sample_stride, ::sample_stride]
    sample_indices = np.indices(ag_sample.shape)
    sample_centers = [
        (sample_indices[i] * sample_stride + 0.5) * voxel_A[i]
        for i in range(3)
    ]
    ax_sample.scatter(
        *[sample_centers[i][ag_sample] for i in range(3)],
        s=2.0,
        c="#74a9cf",
        alpha=0.10,
        marker="s",
        depthshade=False,
        linewidths=0,
        rasterized=True,
    )
    ax_sample.scatter(
        *[sample_centers[i][sc_sample] for i in range(3)],
        s=2.4,
        c="#334e62",
        alpha=0.20,
        marker="s",
        depthshade=False,
        linewidths=0,
        rasterized=True,
    )
    all_liquid = ag & (liquid_fraction > 1.0e-12) & ~largest_mask
    partial_coords = np.where(all_liquid)
    if partial_coords[0].size:
        ax_sample.scatter(
            *((partial_coords[i] + 0.5) * voxel_A[i] for i in range(3)),
            s=34,
            c="#ffd166",
            marker="s",
            edgecolors="#7b5a12",
            linewidths=0.35,
            depthshade=False,
        )
    ax_sample.scatter(
        *((coords[i] + 0.5) * voxel_A[i] for i in range(3)),
        s=34,
        c="#f28e2b",
        marker="s",
        edgecolors="#7f3b08",
        linewidths=0.35,
        depthshade=False,
    )
    draw_box(
        ax_sample,
        ((0.0, domain_size_A[0]), (0.0, domain_size_A[1]), (0.0, domain_size_A[2])),
        "#667788",
        0.7,
    )
    source = np.asarray([domain_size_A[0] / 2, domain_size_A[1] / 2, domain_size_A[2]])
    ax_sample.scatter(*source, color="#2468a2", s=95, marker="*", depthshade=False)
    ax_sample.set_xlim(0, domain_size_A[0])
    ax_sample.set_ylim(0, domain_size_A[1])
    ax_sample.set_zlim(0, domain_size_A[2])
    ax_sample.set_box_aspect(domain_size_A.copy())
    ax_sample.set_proj_type("ortho")
    ax_sample.view_init(elev=24, azim=-53)
    ax_sample.set_xlabel("x (Å)")
    ax_sample.set_ylabel("y (Å)")
    ax_sample.set_zlabel("z (Å)")
    ax_sample.set_title("Complete QSGS sample with the saved liquid region inside")
    ax_sample.legend(
        handles=[
            Line2D([0], [0], marker="s", color="w", markerfacecolor="#74a9cf", label="Ag matrix (every 2nd voxel)"),
            Line2D([0], [0], marker="s", color="w", markerfacecolor="#334e62", label="Ti₃SiC₂ scaffold (every 2nd voxel)"),
            Line2D([0], [0], marker="s", color="w", markerfacecolor="#f28e2b", label="connected liquid Ag, f ≥ 0.5"),
            Line2D([0], [0], marker="s", color="w", markerfacecolor="#ffd166", label="other partially liquid Ag"),
            Line2D([0], [0], marker="*", color="w", markerfacecolor="#2468a2", label="Gaussian source on top face"),
        ],
        loc="upper left",
        fontsize=7.3,
        framealpha=0.9,
    )

    ax_detail = whole.add_subplot(1, 2, 2, projection="3d")
    local_context = roi_ag & ~roi_liquid & ~roi_partial
    local_xyz = np.where(local_context)
    ax_detail.scatter(
        centers[0][local_xyz[0]], centers[1][local_xyz[1]], centers[2][local_xyz[2]],
        s=9, c="#7792a5", alpha=0.18, depthshade=False, linewidths=0,
    )
    if np.any(roi_partial):
        ax_detail.voxels(
            *edge_grid, roi_partial,
            facecolors=(1.0, 0.78, 0.26, 0.88),
            edgecolors=(0.48, 0.30, 0.08, 0.8),
            linewidth=0.35,
        )
    ax_detail.voxels(
        *edge_grid, roi_liquid,
        facecolors="#f28e2b",
        edgecolors="#7f3b08",
        linewidth=0.42,
    )
    ax_detail.set_xlabel("x (Å)")
    ax_detail.set_ylabel("y (Å)")
    ax_detail.set_zlabel("z (Å)")
    ax_detail.set_box_aspect((upper - lower) * voxel_A)
    ax_detail.view_init(elev=27, azim=-54)
    ax_detail.set_proj_type("ortho")
    ax_detail.set_title(f"Pool detail: {int(sizes.max())} connected voxels")

    whole.suptitle(
        f"Ag–Ti₃SiC₂, 70:30 vol% · 1× Case A · t={time_ps:.2f} ps · "
        f"peak T={peak_K:.1f} K · connected liquid cluster={int(sizes.max())} voxels "
        f"({sizes.max() * np.prod(voxel_A):.0f} Å³)\n"
        "The phase cloud is sampled every 2nd voxel for visibility; the highlighted melt voxels use the original grid."
    )
    whole.text(
        0.5, 0.035,
        f"Whole sample: {domain_size_A[0]:.1f} × {domain_size_A[1]:.1f} × {domain_size_A[2]:.1f} Å; "
        "orange marks the main connected liquid Ag cluster within the full dual-phase structure.",
        ha="center", fontsize=9,
    )
    whole.savefig(OUTPUT_WHOLE, dpi=220, facecolor="white", bbox_inches="tight")
    print(f"saved={OUTPUT_WHOLE}")
    print(f"time_ps={time_ps:.8f} peak_temperature_K={peak_K:.6f}")
    print(f"liquid_fraction={liquid_volfrac:.12g} cluster_voxels={int(sizes.max())}")
    print(f"cluster_volume_A3={sizes.max() * np.prod(voxel_A):.6f} centroid_A={centroid_A.tolist()}")


if __name__ == "__main__":
    main()
