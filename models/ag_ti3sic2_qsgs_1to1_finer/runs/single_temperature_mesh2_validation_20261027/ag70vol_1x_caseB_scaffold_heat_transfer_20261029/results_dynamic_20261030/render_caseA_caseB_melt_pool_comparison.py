#!/usr/bin/env python3
"""Render Case A/B melt pools at matched time and identical spatial scale."""

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
CHECKPOINTS = {
    "A": HERE / "checkpoints/caseA_ag70_1x_heat_match_14p51ps.npz",
    "B": HERE / "checkpoints/caseB_ag70_1x_heat_final_14p506ps.npz",
}
OUTPUT = HERE / "caseA_caseB_melt_pool_same_scale_14p51ps.png"
COLORS = ("#203b4d", "#86b9d5", "#ffd43b", "#f28e2b")


def main() -> None:
    sys.path.insert(0, str(PHASE_DIR))
    from qsgs_liquid_tracking import liquid_cluster_labels

    with np.load(GEOMETRY, allow_pickle=False) as data:
        scaffold = np.asarray(data["skeleton"], dtype=bool)
        size_nm = np.asarray(data["size_A"], dtype=float) / 10.0
    ag = ~scaffold
    spacing_nm = size_nm / np.asarray(scaffold.shape, dtype=float)
    states = {}
    boxes = []
    for case, checkpoint in CHECKPOINTS.items():
        with np.load(checkpoint, allow_pickle=False) as data:
            temperature = np.asarray(data["temperature_K"], dtype=float)
            fraction = np.asarray(data["ag_liquid_fraction"], dtype=float)
            meta = json.loads(str(data["metadata_json"].item()))
        liquid = ag & (fraction >= 0.5)
        partial = ag & (fraction > 1e-12) & (fraction < 0.5)
        labels, cluster_count = liquid_cluster_labels(liquid)
        counts = np.bincount(labels.ravel(), minlength=cluster_count + 1)
        counts[0] = 0
        largest = int(counts.max()) if counts.size > 1 else 0
        largest_nm3 = largest * float(np.prod(spacing_nm))
        coords = np.where(liquid)
        if coords[0].size == 0:
            raise ValueError(f"Case {case} has no liquid Ag")
        boxes.append((np.array([axis.min() for axis in coords]), np.array([axis.max() + 1 for axis in coords])))
        phase = np.zeros(scaffold.shape, dtype=np.uint8)
        phase[ag & ~liquid & ~partial] = 1
        phase[partial] = 2
        phase[liquid] = 3
        states[case] = (phase, temperature, meta, largest, largest_nm3)

    lower = np.maximum(np.min(np.stack([b[0] for b in boxes]), axis=0) - 5, 0)
    upper = np.minimum(np.max(np.stack([b[1] for b in boxes]), axis=0) + 5, np.asarray(scaffold.shape))
    edge_idx = [np.r_[np.arange(int(lower[i]), int(upper[i]), 2), int(upper[i])] for i in range(3)]
    center_idx = [(edge[:-1] + edge[1:] - 1) // 2 for edge in edge_idx]
    edges_nm = [edge_idx[i] * spacing_nm[i] for i in range(3)]
    roi_size_nm = (upper - lower) * spacing_nm

    fig = plt.figure(figsize=(15.5, 8.5))
    for panel, case in enumerate(("A", "B"), start=1):
        phase, temperature, meta, largest, volume = states[case]
        roi = phase[tuple(slice(int(lower[i]), int(upper[i])) for i in range(3))]
        local_centers = [center_idx[i] - lower[i] for i in range(3)]
        display = roi[np.ix_(*local_centers)]
        x, y, z = np.meshgrid(*edges_nm, indexing="ij")
        ax = fig.add_subplot(1, 2, panel, projection="3d")
        for code, color in enumerate(COLORS):
            mask = display == code
            if np.any(mask):
                rgba = matplotlib.colors.to_rgba(color, 0.38 if code == 3 else 1.0)
                edge = matplotlib.colors.to_rgba(color, 0.12 if code == 3 else 0.18)
                ax.voxels(x, y, z, mask, facecolors=rgba, edgecolors=edge, linewidth=0.035, shade=True)
        ax.set_xlim(edges_nm[0][0], edges_nm[0][-1])
        ax.set_ylim(edges_nm[1][0], edges_nm[1][-1])
        ax.set_zlim(edges_nm[2][0], edges_nm[2][-1])
        ax.set_box_aspect(roi_size_nm)
        ax.view_init(elev=28, azim=-55)
        ax.set_proj_type("ortho")
        ax.set_xlabel("x (nm)"); ax.set_ylabel("y (nm)"); ax.set_zlabel("z (nm)")
        boundary = "scaffold thermally inactive" if case == "A" else "scaffold stores and conducts heat"
        ax.set_title(
            f"Case {case}: {boundary}\n"
            f"t={meta['time_ps']:.3f} ps · Ag peak={temperature[ag].max():.1f} K\n"
            f"largest cluster={largest:,} voxels · {volume:.1f} nm³"
        )
    labels = ("Ti₃SiC₂ scaffold", "unmelted Ag", "partially liquid Ag", "liquid Ag (translucent)")
    fig.legend(
        handles=[Patch(facecolor=matplotlib.colors.to_rgba(c, 0.38 if i == 3 else 1.0), edgecolor="#36434c", label=labels[i]) for i, c in enumerate(COLORS)],
        loc="upper center", bbox_to_anchor=(0.5, 0.94), ncol=4, framealpha=0.96,
    )
    fig.suptitle("Ag70 Case A vs Case B: matched-time melt-pool comparison", fontsize=16, y=0.99)
    fig.text(
        0.5, 0.015,
        "Same geometry and coordinates; both panels use the same ROI limits. "
        "One representative native voxel per 2×2×2 region is displayed; listed cluster metrics use the full native grid.",
        ha="center", fontsize=9,
    )
    fig.subplots_adjust(left=0.02, right=0.98, top=0.80, bottom=0.07, wspace=0.02)
    fig.savefig(OUTPUT, dpi=190, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"saved={OUTPUT}")
    print(f"shared_roi_origin_nm={(lower * spacing_nm).tolist()} size_nm={roi_size_nm.tolist()}")


if __name__ == "__main__":
    main()
