#!/usr/bin/env python3
"""Render all four phase classes as solid 3D voxel blocks with a cutaway."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from scipy import ndimage


BASE = Path(__file__).resolve().parent
RUN_DIR = BASE.parent
GEOMETRY = RUN_DIR / "qsgs_morphology_ag70vol.npz"
CHECKPOINT = BASE / "checkpoints" / "caseA_ag70_1x_heat_18p636ps.npz"
OUTPUT = BASE / "ag70_phase_map_3d_18p636ps.png"

COLORS = {
    "scaffold": "#263f51",
    "solid_ag": "#83b7d5",
    "partial_ag": "#ffd166",
    "liquid_ag": "#f28e2b",
}


def add_sample_outline(ax, size_nm: np.ndarray) -> None:
    corners = np.array([
        [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
        [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
    ], dtype=float) * size_nm
    edges = ((0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6),
             (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7))
    for i, j in edges:
        ax.plot(*corners[[i, j]].T, color="#536574", linewidth=0.75, alpha=0.78)
def main() -> None:
    if not GEOMETRY.exists() or not CHECKPOINT.exists():
        raise FileNotFoundError(f"Required geometry/checkpoint not found: {CHECKPOINT}")
    with np.load(GEOMETRY, allow_pickle=False) as data:
        scaffold = np.asarray(data["skeleton"], dtype=bool)
        size_nm = np.asarray(data["size_A"], dtype=float) / 10.0
    with np.load(CHECKPOINT, allow_pickle=False) as data:
        temperature = np.asarray(data["temperature_K"], dtype=float)
        liquid_fraction = np.asarray(data["ag_liquid_fraction"], dtype=float)
        metadata = json.loads(str(data["metadata_json"].item()))

    ag = ~scaffold
    spacing_nm = size_nm / np.asarray(scaffold.shape, dtype=float)
    liquid = ag & (liquid_fraction >= 0.5)
    partial = ag & (liquid_fraction > 1e-12) & (liquid_fraction < 0.5)
    solid_ag = ag & ~liquid & ~partial
    labels, _ = ndimage.label(liquid, ndimage.generate_binary_structure(3, 1))
    counts = np.bincount(labels.ravel())
    counts[0] = 0
    pool_id = int(np.argmax(counts))
    pool = labels == pool_id
    pool_voxels = int(counts[pool_id])
    pool_volume = pool_voxels * float(np.prod(spacing_nm))
    coords = np.where(pool)
    lo = np.array([c.min() for c in coords])
    hi = np.array([c.max() for c in coords])
    bbox = (hi - lo + 1) * spacing_nm
    time_ps = float(metadata["time_ps"])

    # Discrete native-grid phase values, with the same four labels as the source slice.
    phase = np.zeros(scaffold.shape, dtype=np.uint8)
    phase[solid_ag] = 1
    phase[partial] = 2
    phase[liquid] = 3

    # Render one representative native voxel per 4x4x4 block. This makes the
    # all-phase solid-voxel rendering light enough to run alongside the study.
    native_shape = np.asarray(phase.shape)
    native_edges = [np.r_[np.arange(0, n, 4), n] for n in native_shape]
    native_centres = [(edges[:-1] + edges[1:] - 1) // 2 for edges in native_edges]
    phase_display = phase[np.ix_(*native_centres)]
    masks_display = (
        phase_display == 0,
        phase_display == 1,
        phase_display == 2,
        phase_display == 3,
    )

    # Remove one corner of the rendering only. This exposes interior voxels on
    # three faces without changing the checkpoint or any reported metric.
    visible = np.ones(phase_display.shape, dtype=bool)
    nx, ny, nz = phase_display.shape
    visible[nx // 2 :, ny // 2 :, nz // 2 :] = False

    x, y, z = np.meshgrid(
        native_edges[0] * spacing_nm[0],
        native_edges[1] * spacing_nm[1],
        native_edges[2] * spacing_nm[2],
        indexing="ij",
    )

    fig = plt.figure(figsize=(15.0, 11.0))
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    ax.set_proj_type("ortho")
    phase_specs = (
        (masks_display[0], COLORS["scaffold"], "Ti₃SiC₂ scaffold"),
        (masks_display[1], COLORS["solid_ag"], "unmelted Ag"),
        (masks_display[2], COLORS["partial_ag"], "partially liquid Ag"),
        (masks_display[3], COLORS["liquid_ag"], "liquid Ag"),
    )
    for mask, color, _label in phase_specs:
        ax.voxels(
            x, y, z, mask & visible,
            facecolors=color,
            edgecolors=(0.16, 0.22, 0.27, 0.16),
            linewidth=0.045,
            shade=True,
        )

    add_sample_outline(ax, size_nm)
    ax.set_xlim(0, size_nm[0])
    ax.set_ylim(0, size_nm[1])
    ax.set_zlim(0, size_nm[2])
    ax.set_box_aspect(size_nm)
    ax.set_xlabel("x (nm)", labelpad=8)
    ax.set_ylabel("y (nm)", labelpad=8)
    ax.set_zlabel("z (nm)", labelpad=8)
    ax.view_init(elev=29, azim=-56)
    ax.grid(True, alpha=0.2)
    ax.legend(handles=[Patch(facecolor=color, edgecolor="#36434c", label=label)
                       for _mask, color, label in phase_specs],
              loc="upper left", bbox_to_anchor=(0.02, 0.98), framealpha=0.96, fontsize=10)
    ax.text2D(
        0.02, 0.03,
        f"One representative voxel from each 4 × 4 × 4 native region is rendered. "
        f"Largest liquid cluster = {pool_volume:.0f} nm³, "
        f"{pool_voxels:,} voxels.\n"
        f"Native-grid pool dimensions: {bbox[0]:.1f} × {bbox[1]:.1f} × {bbox[2]:.1f} nm. "
        "One corner is hidden to expose interior phases; simulation data are unchanged.",
        transform=ax.transAxes, fontsize=9, ha="left", va="bottom",
        bbox={"boxstyle": "round,pad=0.45", "facecolor": "white", "alpha": 0.88, "edgecolor": "#c8d0d6"},
    )
    fig.suptitle(
        f"3D phase voxels and Ag melt pool | Ag–Ti₃SiC₂ 70:30 vol% | t={time_ps:.3f} ps\n"
        f"Peak temperature {temperature.max():.1f} K | largest pool={pool_voxels:,} native voxels",
        fontsize=17, y=0.96,
    )
    fig.text(
        0.5, 0.012,
        "Visualization cutaway only; thermal enthalpy model, with no liquid flow or splashing physics.",
        ha="center", fontsize=9,
    )
    fig.subplots_adjust(left=0.02, right=0.98, top=0.88, bottom=0.06)
    fig.savefig(OUTPUT, dpi=200, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"saved={OUTPUT}")
    print(f"t_ps={time_ps:.6f} peak_K={temperature.max():.4f} pool_voxels={pool_voxels}")
    print(f"pool_volume_nm3={pool_volume:.6f} bbox_nm={bbox.tolist()}")


if __name__ == "__main__":
    main()
