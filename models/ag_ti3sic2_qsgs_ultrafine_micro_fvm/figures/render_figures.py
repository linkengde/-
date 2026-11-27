#!/usr/bin/env python3
"""用真实二值体素数据生成论文级切面和三维正交切片图。"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures"
SEEDS = (1105, 1136, 1182)
MODELS = {
    "ultrafine_500nm": {"label": "Ultrafine scaffold", "size_um": 10.0},
    "micro_10um": {"label": "Micron scaffold", "size_um": 200.0},
}
AG_RGB = np.array([0.72, 0.75, 0.79])
SC_RGB = np.array([0.08, 0.36, 0.66])
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.titlesize": 10, "axes.labelsize": 9,
    "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white",
})


def load_geometry(seed: int) -> np.ndarray:
    with np.load(ROOT / "geometry_generation" / f"qsgs_seed_{seed}.npz", allow_pickle=False) as f:
        return np.asarray(f["skeleton"], dtype=bool)


def add_cutaway(ax, skeleton: np.ndarray, size_um: float, mode: str = "both") -> None:
    """绘制三个真实中心切面；不伪装成体素外表面或SEM重建。"""
    n = skeleton.shape[0]
    step = size_um / n
    middle = n // 2
    plane_specs = [
        (2, middle, "xy"),
        (1, middle, "xz"),
        (0, middle, "yz"),
    ]
    for fixed_axis, fixed_index, plane_name in plane_specs:
        if plane_name == "xy":
            image = skeleton[:, :, fixed_index]
            rows, cols = np.indices(image.shape)
            origins = np.column_stack((rows.ravel(), cols.ravel(), np.full(image.size, fixed_index)))
            corners = ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0))
        elif plane_name == "xz":
            image = skeleton[:, fixed_index, :]
            rows, cols = np.indices(image.shape)
            origins = np.column_stack((rows.ravel(), np.full(image.size, fixed_index), cols.ravel()))
            corners = ((0, 0, 0), (1, 0, 0), (1, 0, 1), (0, 0, 1))
        else:
            image = skeleton[fixed_index, :, :]
            rows, cols = np.indices(image.shape)
            origins = np.column_stack((np.full(image.size, fixed_index), rows.ravel(), cols.ravel()))
            corners = ((0, 0, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1))

        phase = image.ravel()
        if mode == "skeleton":
            visible = phase
            colors = np.tile(np.r_[SC_RGB, 1.0], (int(visible.sum()), 1))
        elif mode == "ag":
            visible = ~phase
            colors = np.tile(np.r_[AG_RGB, 1.0], (int(visible.sum()), 1))
        else:
            visible = np.ones(phase.size, dtype=bool)
            rgba = np.empty((phase.size, 4), dtype=float)
            rgba[phase, :3] = SC_RGB
            rgba[~phase, :3] = AG_RGB
            rgba[:, 3] = 0.28 if mode == "translucent_ag" else 0.98
            colors = rgba

        bases = origins[visible].astype(float) * step
        if not len(bases):
            continue
        face_offsets = np.asarray(corners, dtype=float) * step
        faces = bases[:, None, :] + face_offsets[None, :, :]
        collection = Poly3DCollection(
            faces, facecolors=colors, edgecolors="none", linewidths=0,
            antialiased=False, zsort="average",
        )
        ax.add_collection3d(collection)

    ax.set_xlim(0, size_um)
    ax.set_ylim(0, size_um)
    ax.set_zlim(0, size_um)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=24, azim=-56, roll=0)
    ax.set_xlabel("x (μm)", labelpad=1)
    ax.set_ylabel("y (μm)", labelpad=1)
    ax.set_zlabel("z (μm)", labelpad=1)
    ax.tick_params(pad=0, labelsize=7)
    ax.grid(False)


def save_volume_panels(skeleton: np.ndarray, model_name: str, size_um: float) -> None:
    fig = plt.figure(figsize=(10.5, 8.3), constrained_layout=True)
    panels = [
        ("both", "Ag + Ti₃SiC₂: orthogonal center cutaway"),
        ("skeleton", "Ti₃SiC₂ scaffold only"),
        ("ag", "Ag phase only"),
        ("translucent_ag", "Ag translucent; scaffold visible"),
    ]
    for index, (mode, title) in enumerate(panels, start=1):
        ax = fig.add_subplot(2, 2, index, projection="3d")
        add_cutaway(ax, skeleton, size_um, mode)
        ax.set_title(title + f"\n{size_um:g} μm cube, 150³ cells", pad=2)
    fig.suptitle(f"{MODELS[model_name]['label']} · seed 1105 · true voxel cutaway", fontsize=13)
    fig.text(0.5, -0.015, "Three mutually perpendicular center planes from the generated binary volume; not an external-surface rendering.",
             ha="center", fontsize=8, color="#40454b")
    fig.savefig(OUT / f"{model_name}_volume_views.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_slices(skeleton: np.ndarray, model_name: str, size_um: float) -> None:
    n = skeleton.shape[0]
    mid = n // 2
    slices = [skeleton[:, :, mid], skeleton[:, mid, :], skeleton[mid, :, :]]
    labels = ["XY center plane", "XZ center plane", "YZ center plane"]
    extent = (0, size_um, 0, size_um)
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.8), constrained_layout=True)
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap([AG_RGB, SC_RGB])
    for ax, image, label in zip(axes, slices, labels):
        ax.imshow(image.T, origin="lower", interpolation="nearest", extent=extent, cmap=cmap, vmin=0, vmax=1)
        ax.set_title(label)
        ax.set_xlabel("μm")
        ax.set_ylabel("μm")
        ax.set_aspect("equal")
    fig.suptitle(f"{MODELS[model_name]['label']} · seed 1105 · 150×150 samples per plane")
    fig.savefig(OUT / f"{model_name}_center_slices.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_closeup(skeleton: np.ndarray, model_name: str, size_um: float) -> None:
    n = skeleton.shape[0]
    center = n // 2
    half = 15
    region = skeleton[center-half:center+half, center-half:center+half, center]
    cell_nm = size_um * 1000 / n
    width_um = region.shape[0] * size_um / n
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap([AG_RGB, SC_RGB])
    fig, ax = plt.subplots(figsize=(5, 4.5), constrained_layout=True)
    ax.imshow(region.T, origin="lower", interpolation="nearest", extent=(0, width_um, 0, width_um), cmap=cmap, vmin=0, vmax=1)
    ax.set_xlabel("μm")
    ax.set_ylabel("μm")
    ax.set_title(f"Center XY region · 30×30 cells · cell pitch {cell_nm:.3g} nm")
    ax.set_aspect("equal")
    fig.suptitle(f"{MODELS[model_name]['label']} · seed 1105 local detail")
    fig.savefig(OUT / f"{model_name}_local_detail.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_width_distributions() -> None:
    colors = {"skeleton": SC_RGB, "ag_channel": AG_RGB}
    titles = {"skeleton": "Ti₃SiC₂ maximal-inscribed-sphere diameter", "ag_channel": "Ag channel maximal-inscribed-sphere diameter"}
    for model_name, model in MODELS.items():
        fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
        for ax, phase in zip(axes, ("skeleton", "ag_channel")):
            for seed in SEEDS:
                # The JSON contains center-weighted quantiles; plot exact distributions from the source voxels.
                with np.load(ROOT / "geometry_generation" / f"qsgs_seed_{seed}.npz", allow_pickle=False) as f:
                    mask = np.asarray(f["skeleton"], dtype=bool)
                phase_mask = mask if phase == "skeleton" else ~mask
                from scipy import ndimage
                distance = ndimage.distance_transform_edt(phase_mask)
                maxima = phase_mask & (distance >= ndimage.maximum_filter(distance, size=3, mode="constant") - 1e-12)
                widths = 2 * distance[maxima] * model["size_um"] / 150
                ax.hist(widths, bins=45, density=True, histtype="step", linewidth=1.3,
                        color=colors[phase], alpha=0.6 + 0.12 * (seed - SEEDS[0]) / (SEEDS[-1] - SEEDS[0]),
                        label=f"seed {seed}; P50={np.median(widths):.3g} μm")
            ax.set_title(titles[phase])
            ax.set_xlabel("diameter (μm)")
            ax.set_ylabel("density per μm")
            ax.legend(fontsize=7, frameon=False)
            ax.grid(axis="y", alpha=0.2)
        fig.suptitle(f"{model['label']} · sphere-center weighted distributions, three seeds")
        fig.savefig(OUT / f"{model_name}_width_distributions.png", dpi=300, bbox_inches="tight")
        plt.close(fig)


def save_interface_comparison() -> None:
    means = []
    values = {name: [] for name in MODELS}
    for seed in SEEDS:
        data = json.loads((ROOT / "geometry_audit" / f"seed_{seed}_audit.json").read_text())
        for model_name in MODELS:
            values[model_name].append(data["interface_area_density_scale_pair_m2_per_m3"][model_name] / 1e6)
    fig, ax = plt.subplots(figsize=(6.5, 4), constrained_layout=True)
    positions = [0, 1]
    box = ax.boxplot([values[k] for k in MODELS], positions=positions, widths=0.48,
                     patch_artist=True, showmeans=True)
    for patch, color in zip(box["boxes"], ("#4489ba", "#8a9eae")):
        patch.set_facecolor(color)
        patch.set_alpha(0.65)
    for x, model_name in zip(positions, MODELS):
        for index, value in enumerate(values[model_name]):
            ax.scatter(x + (-0.07 + index * 0.07), value, s=24, color="#153a54", zorder=3)
    ax.set_xticks(positions, ["10 μm cube\n500 nm target", "200 μm cube\n10 μm target"])
    ax.set_ylabel("internal interface area density (10⁶ m²/m³)")
    ax.set_title("Same voxel geometry scaled by 20×: interface density scales by 1/20")
    ax.grid(axis="y", alpha=0.25)
    fig.savefig(OUT / "interface_area_density_comparison.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_scale_side_by_side() -> None:
    skeleton = load_geometry(1105)
    fig = plt.figure(figsize=(10, 4.7), constrained_layout=True)
    for i, (model_name, model) in enumerate(MODELS.items(), start=1):
        ax = fig.add_subplot(1, 2, i, projection="3d")
        add_cutaway(ax, skeleton, model["size_um"], "both")
        ax.set_title(f"{model['label']}\n{model['size_um']:g} μm cube · Δx={model['size_um']*1000/150:.3g} nm")
    fig.suptitle("Geometrically identical binary RVE scaled by 20×", fontsize=13)
    fig.text(0.5, -0.012, "Axes are labelled in each model's physical micrometre scale; apparent cube sizes are normalized for comparison.",
             ha="center", fontsize=8, color="#40454b")
    fig.savefig(OUT / "scale_pair_side_by_side.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    representative = load_geometry(1105)
    for model_name, model in MODELS.items():
        save_volume_panels(representative, model_name, model["size_um"])
        save_slices(representative, model_name, model["size_um"])
        save_closeup(representative, model_name, model["size_um"])
    save_width_distributions()
    save_interface_comparison()
    save_scale_side_by_side()
    print(json.dumps({"png_count": len(list(OUT.glob("*.png"))), "dpi": 300}))


if __name__ == "__main__":
    main()
