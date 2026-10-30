#!/usr/bin/env python3
"""Render the final heated Ag70 Case B state as whole-sample and pool views."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch


HERE = Path(__file__).resolve().parent
RUN_DIR = HERE.parent.parent
PROJECT_DIR = RUN_DIR.parents[1]
PHASE_DIR = PROJECT_DIR / "phase_field_qsgs"
GEOMETRY = RUN_DIR / "qsgs_morphology_ag70vol.npz"
CHECKPOINT = HERE / "checkpoints" / "caseB_ag70_1x_heat_final_14p506ps.npz"
OUTPUT = HERE / "caseB_ag70_1x_melt_pool_3d_14p506ps.png"

COLORS = {
    "scaffold": "#203b4d",
    "solid_ag": "#86b9d5",
    "partial_ag": "#ffd43b",
    "liquid_ag": "#f28e2b",
}


def draw_box(ax, bounds: tuple[tuple[float, float], ...], color: str, width: float) -> None:
    (x0, x1), (y0, y1), (z0, z1) = bounds
    corners = np.array(
        [
            [x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],
            [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1],
        ]
    )
    edges = ((0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
             (0, 4), (1, 5), (2, 6), (3, 7))
    for i, j in edges:
        ax.plot(*corners[[i, j]].T, color=color, linewidth=width)


def render_phases(ax, phase: np.ndarray, edges: list[np.ndarray], visible: np.ndarray) -> None:
    x, y, z = np.meshgrid(*edges, indexing="ij")
    specs = (
        (0, COLORS["scaffold"], "Ti₃SiC₂ scaffold", 1.0),
        (1, COLORS["solid_ag"], "unmelted Ag", 1.0),
        (2, COLORS["partial_ag"], "partially liquid Ag", 1.0),
        (3, COLORS["liquid_ag"], "liquid Ag", 0.36),
    )
    for code, color, _label, alpha in specs:
        mask = (phase == code) & visible
        if not np.any(mask):
            continue
        rgba = matplotlib.colors.to_rgba(color, alpha)
        edge_rgba = matplotlib.colors.to_rgba(color, 0.12 if code == 3 else 0.18)
        ax.voxels(x, y, z, mask, facecolors=rgba, edgecolors=edge_rgba, linewidth=0.035, shade=True)


def main() -> None:
    if not GEOMETRY.is_file() or not CHECKPOINT.is_file():
        raise FileNotFoundError(f"Missing Ag70 geometry or Case B checkpoint: {CHECKPOINT}")
    sys.path.insert(0, str(PHASE_DIR))
    from qsgs_liquid_tracking import liquid_cluster_labels

    with np.load(GEOMETRY, allow_pickle=False) as data:
        scaffold = np.asarray(data["skeleton"], dtype=bool)
        size_nm = np.asarray(data["size_A"], dtype=float) / 10.0
    with np.load(CHECKPOINT, allow_pickle=False) as data:
        temperature = np.asarray(data["temperature_K"], dtype=float)
        liquid_fraction = np.asarray(data["ag_liquid_fraction"], dtype=float)
        metadata = json.loads(str(data["metadata_json"].item()))

    ag = ~scaffold
    liquid = ag & (liquid_fraction >= 0.5)
    partial = ag & (liquid_fraction > 1e-12) & (liquid_fraction < 0.5)
    labels, count = liquid_cluster_labels(liquid)
    sizes = np.bincount(labels.ravel(), minlength=count + 1)
    sizes[0] = 0
    if count == 0 or int(sizes.max()) == 0:
        raise ValueError("The selected checkpoint contains no connected liquid Ag cluster.")
    largest_id = int(np.argmax(sizes))
    pool = labels == largest_id
    pool_voxels = int(sizes[largest_id])
    spacing_nm = size_nm / np.asarray(scaffold.shape, dtype=float)
    pool_volume_nm3 = pool_voxels * float(np.prod(spacing_nm))
    coords = np.where(pool)
    lower = np.maximum(np.array([axis.min() for axis in coords]) - 5, 0)
    upper = np.minimum(np.array([axis.max() + 6 for axis in coords]), np.asarray(scaffold.shape))
    roi_slices = tuple(slice(int(lo), int(hi)) for lo, hi in zip(lower, upper, strict=True))
    roi_phase = np.zeros(tuple(upper - lower), dtype=np.uint8)
    roi_ag = ag[roi_slices]
    roi_phase[roi_ag & ~liquid[roi_slices] & ~partial[roi_slices]] = 1
    roi_phase[partial[roi_slices]] = 2
    roi_phase[liquid[roi_slices]] = 3

    # Overview: one representative native voxel per 4×4×4 block so that the
    # full-sample cutaway stays quick to render while retaining the topology.
    overview_edges_idx = [np.r_[np.arange(0, n, 4), n] for n in scaffold.shape]
    overview_centers = [(e[:-1] + e[1:] - 1) // 2 for e in overview_edges_idx]
    full_phase = np.zeros(scaffold.shape, dtype=np.uint8)
    full_phase[ag & ~liquid & ~partial] = 1
    full_phase[partial] = 2
    full_phase[liquid] = 3
    overview = full_phase[np.ix_(*overview_centers)]
    overview_edges = [e * spacing_nm[i] for i, e in enumerate(overview_edges_idx)]
    roi_edges = [
        np.arange(int(lower[i]), int(upper[i]) + 1, dtype=float) * spacing_nm[i]
        for i in range(3)
    ]
    domain_nm = size_nm
    pool_center_nm = (np.array([axis.mean() + 0.5 for axis in coords]) * spacing_nm)
    bbox_nm = (upper - lower) * spacing_nm
    time_ps = float(metadata["time_ps"])

    fig = plt.figure(figsize=(17, 9.5))
    ax_full = fig.add_subplot(1, 2, 1, projection="3d")
    full_visible = np.ones(overview.shape, dtype=bool)
    full_visible[overview.shape[0] // 2 :, overview.shape[1] // 2 :, overview.shape[2] // 2 :] = False
    render_phases(ax_full, overview, overview_edges, full_visible)
    draw_box(
        ax_full,
        tuple((float(lower[i] * spacing_nm[i]), float(upper[i] * spacing_nm[i])) for i in range(3)),
        "#e4572e",
        2.0,
    )
    ax_full.scatter(*pool_center_nm, color="#d94801", marker="o", s=45, depthshade=False)
    ax_full.set_xlim(0, domain_nm[0]); ax_full.set_ylim(0, domain_nm[1]); ax_full.set_zlim(0, domain_nm[2])
    ax_full.set_box_aspect(domain_nm)
    ax_full.view_init(elev=28, azim=-55)
    ax_full.set_proj_type("ortho")
    ax_full.set_title("Full sample overview (4×4×4 representative voxels)\nRed box marks the magnified pool region")
    ax_full.set_xlabel("x (nm)"); ax_full.set_ylabel("y (nm)"); ax_full.set_zlabel("z (nm)")

    ax_roi = fig.add_subplot(1, 2, 2, projection="3d")
    roi_visible = np.ones(roi_phase.shape, dtype=bool)
    # A small corner cutaway exposes the phases behind the translucent liquid.
    if min(roi_phase.shape) > 4:
        roi_visible[roi_phase.shape[0] // 2 :, roi_phase.shape[1] // 2 :, roi_phase.shape[2] // 2 :] = False
    render_phases(ax_roi, roi_phase, roi_edges, roi_visible)
    ax_roi.set_box_aspect(bbox_nm)
    ax_roi.view_init(elev=28, azim=-55)
    ax_roi.set_proj_type("ortho")
    ax_roi.set_title(
        f"Melt pool and surrounding phases (native-grid ROI)\n{pool_voxels:,} liquid voxels · {pool_volume_nm3:.3f} nm³"
    )
    ax_roi.set_xlabel("x (nm)"); ax_roi.set_ylabel("y (nm)"); ax_roi.set_zlabel("z (nm)")

    specs = (
        (COLORS["scaffold"], "Ti₃SiC₂ scaffold"),
        (COLORS["solid_ag"], "unmelted Ag"),
        (COLORS["partial_ag"], "partially liquid Ag"),
        (matplotlib.colors.to_rgba(COLORS["liquid_ag"], 0.36), "liquid Ag (translucent)"),
    )
    fig.legend(
        handles=[Patch(facecolor=c, edgecolor="#34495e", label=label) for c, label in specs],
        loc="upper center", bbox_to_anchor=(0.5, 0.93), ncol=4, framealpha=0.96,
    )
    fig.suptitle(
        f"Ag70 1× Case B 3D melt pool | t={time_ps:.3f} ps | Ag peak T={temperature.max():.1f} K",
        fontsize=16, y=0.99,
    )
    fig.text(
        0.5, 0.015,
        "Phases: Ti₃SiC₂ scaffold, unmelted Ag, partially liquid Ag, and liquid Ag. "
        "Liquid Ag is translucent; the local view uses a visual cutaway. "
        "The overview samples one voxel per 4×4×4 region; the local view uses native voxels.",
        ha="center", fontsize=9,
    )
    fig.subplots_adjust(left=0.02, right=0.98, top=0.83, bottom=0.06, wspace=0.04)
    fig.savefig(OUTPUT, dpi=190, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"saved={OUTPUT}")
    print(f"time_ps={time_ps:.6f} pool_voxels={pool_voxels} volume_nm3={pool_volume_nm3:.6f}")
    print(f"roi_shape={roi_phase.shape} roi_size_nm={bbox_nm.tolist()}")


if __name__ == "__main__":
    main()
