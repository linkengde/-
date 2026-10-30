#!/usr/bin/env python3
"""Render exact QSGS phase voxels and a connected-pool close-up from a checkpoint.

The script is read-only with respect to simulation inputs. Give it a stable
checkpoint with --checkpoint; it writes a PNG and a JSON audit beside it.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-ag70-qsgs")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/fontconfig-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgba
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import proj3d


HERE = Path(__file__).resolve().parent
OUTPUT_DEFAULT = HERE / "ag70_caseB_2x_full_sample_and_pool_t4p04ps.png"

PHASES = (
    (0, "Ti3SiC2 骨架", "#183b50", 1.0),
    (1, "未熔 Ag（f_l = 0）", "#8bbbd5", 0.98),
    (2, "部分液化 Ag（0 < f_l < 0.5）", "#ffdf00", 0.98),
    (3, "液相占优 Ag（f_l >= 0.5）", "#f28e2b", 0.40),
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def make_phase_ids(skeleton: np.ndarray, liquid_fraction: np.ndarray) -> np.ndarray:
    phase = np.empty(skeleton.shape, dtype=np.uint8)
    phase[skeleton] = 0
    ag = ~skeleton
    phase[ag & (liquid_fraction <= 1.0e-12)] = 1
    phase[ag & (liquid_fraction > 1.0e-12) & (liquid_fraction < 0.5)] = 2
    phase[ag & (liquid_fraction >= 0.5)] = 3
    return phase


def exposed_voxel_faces(
    phase_ids: np.ndarray,
    spacing_nm: np.ndarray,
    origin_nm: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Build axis-aligned quadrilateral faces at each categorical voxel boundary."""
    origin = np.zeros(3, dtype=float) if origin_nm is None else np.asarray(origin_nm, dtype=float)
    polygons: list[np.ndarray] = []
    color_ids: list[np.ndarray] = []
    size = phase_ids.shape
    for phase_id, _, _, _ in PHASES:
        occupied = phase_ids == phase_id
        for axis in range(3):
            other = [dimension for dimension in range(3) if dimension != axis]
            for side in (-1, 1):
                neighbor = np.roll(occupied, 1 if side < 0 else -1, axis=axis)
                boundary = [slice(None)] * 3
                boundary[axis] = 0 if side < 0 else -1
                neighbor[tuple(boundary)] = False
                face_indices = np.argwhere(occupied & ~neighbor)
                if not face_indices.size:
                    continue

                # Each voxel face is a flat quad. Its corners follow the grid
                # edges exactly, so the image is a blocky QSGS voxel surface,
                # not a smoothed mesh or a cloud of atom-like markers.
                quads = np.repeat(face_indices[:, None, :], 4, axis=1).astype(np.float32)
                quads += origin.astype(np.float32)[None, None, :] / spacing_nm.astype(np.float32)[None, None, :]
                quads[:, :, axis] += 1.0 if side > 0 else 0.0
                corners = ((0, 0), (1, 0), (1, 1), (0, 1))
                for corner_id, (u, v) in enumerate(corners):
                    quads[:, corner_id, other[0]] += u
                    quads[:, corner_id, other[1]] += v
                quads *= spacing_nm.astype(np.float32)[None, None, :]
                polygons.append(quads)
                color_ids.append(np.full(face_indices.shape[0], phase_id, dtype=np.uint8))

    return np.concatenate(polygons, axis=0), np.concatenate(color_ids, axis=0)


def add_voxel_surface(ax, phase_ids: np.ndarray, spacing_nm: np.ndarray, origin_nm: np.ndarray | None = None) -> None:
    polygons, color_ids = exposed_voxel_faces(phase_ids, spacing_nm, origin_nm)
    base_colors = {phase_id: to_rgba(color, alpha) for phase_id, _, color, alpha in PHASES}
    rgba = np.empty((color_ids.size, 4), dtype=np.float32)
    for phase_id, color in base_colors.items():
        rgba[color_ids == phase_id] = color
    mesh = Poly3DCollection(
        polygons,
        facecolors=rgba,
        edgecolors="none",
        linewidths=0.0,
        zsort="average",
        rasterized=True,
    )
    ax.add_collection3d(mesh)


def configure_axes(ax, lo_nm: np.ndarray, hi_nm: np.ndarray, aspect_nm: np.ndarray) -> None:
    ax.set_xlim(float(lo_nm[0]), float(hi_nm[0]))
    ax.set_ylim(float(lo_nm[1]), float(hi_nm[1]))
    ax.set_zlim(float(lo_nm[2]), float(hi_nm[2]))
    ax.set_box_aspect(aspect_nm)
    ax.set_proj_type("ortho")
    ax.set_xlabel("X (nm)", labelpad=5)
    ax.set_ylabel("Y (nm)", labelpad=5)
    ax.set_zlabel("Z (nm)", labelpad=4)
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.yaxis.set_major_locator(MultipleLocator(5))
    ax.zaxis.set_major_locator(MultipleLocator(2.5))
    ax.tick_params(labelsize=8, pad=0)
    ax.grid(True, color="#c8cdd2", linewidth=0.5, alpha=0.45)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor((0.97, 0.98, 0.99, 1.0))
        axis.pane.set_edgecolor("#aab3ba")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", type=Path, required=True, help="Case B solver checkpoint NPZ")
    parser.add_argument("--geometry", type=Path, required=True, help="native Ag70 QSGS morphology NPZ")
    parser.add_argument("--output", type=Path, default=OUTPUT_DEFAULT)
    parser.add_argument("--dpi", type=int, default=240)
    args = parser.parse_args()
    checkpoint = args.checkpoint.resolve()
    geometry_path = args.geometry.resolve()
    output = args.output.resolve()
    if not checkpoint.is_file() or not geometry_path.is_file():
        parser.error("checkpoint and native QSGS geometry must both exist")
    output.parent.mkdir(parents=True, exist_ok=True)

    with np.load(checkpoint, allow_pickle=False) as archive:
        metadata = json.loads(archive["metadata_json"].item())
        temperature_K = np.asarray(archive["temperature_K"], dtype=np.float32)
        liquid_fraction = np.asarray(archive["ag_liquid_fraction"], dtype=np.float32)
        track_labels = np.asarray(archive["tracker_track_labels"], dtype=np.int32)
    with np.load(geometry_path, allow_pickle=False) as archive:
        skeleton = np.asarray(archive["skeleton"], dtype=bool)
        domain_size_A = np.asarray(archive["size_A"], dtype=float)

    if skeleton.shape != liquid_fraction.shape:
        raise ValueError("QSGS geometry and checkpoint array shapes differ")
    identity = metadata["identity"]
    if identity.get("geometry_sha256") != sha256_file(geometry_path):
        raise ValueError("Checkpoint geometry hash does not match the supplied native QSGS geometry")

    spacing_A = np.asarray(identity["spacing_A"], dtype=float)
    spacing_nm = spacing_A / 10.0
    domain_size_nm = domain_size_A / 10.0
    ag_mask = ~skeleton
    phase_ids = make_phase_ids(skeleton, liquid_fraction)
    ag_liquid_fraction = float(np.sum(liquid_fraction[ag_mask]) / np.count_nonzero(ag_mask))

    ids, counts = np.unique(track_labels[track_labels > 0], return_counts=True)
    if not ids.size:
        raise ValueError("Checkpoint contains no connected liquid Ag cluster")
    largest_index = int(np.argmax(counts))
    largest_track_id = int(ids[largest_index])
    largest_voxel_count = int(counts[largest_index])
    largest_cluster_volume_A3 = float(largest_voxel_count * np.prod(spacing_A))
    largest_cluster_mask = track_labels == largest_track_id
    pool_coords = np.argwhere(largest_cluster_mask)
    pool_lo = pool_coords.min(axis=0)
    pool_hi = pool_coords.max(axis=0) + 1

    # Keep a fixed physical margin around the largest tracked cluster. Z is
    # clipped at the actual sample top, while X/Y stay local to the pool.
    margin_A = np.asarray([45.0, 45.0, 35.0])
    margin_voxels = np.ceil(margin_A / spacing_A).astype(int)
    crop_lo = np.maximum(pool_lo - margin_voxels, 0)
    crop_hi = np.minimum(pool_hi + margin_voxels, np.asarray(phase_ids.shape))
    crop = tuple(slice(int(a), int(b)) for a, b in zip(crop_lo, crop_hi))
    crop_ids = phase_ids[crop]
    crop_origin_nm = crop_lo * spacing_nm
    crop_size_nm = (crop_hi - crop_lo) * spacing_nm

    colors = {phase_id: to_rgba(color, alpha) for phase_id, _, color, alpha in PHASES}
    cjk_regular = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
    cjk_bold = "/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc"
    font_body = FontProperties(fname=cjk_regular, size=10)
    font_panel = FontProperties(fname=cjk_bold, size=14)
    font_heading = FontProperties(fname=cjk_bold, size=19)
    font_subtitle = FontProperties(fname=cjk_regular, size=12)
    font_legend = FontProperties(fname=cjk_regular, size=11)
    rc = {
        "axes.unicode_minus": False,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    }
    with plt.rc_context(rc):
        fig = plt.figure(figsize=(19.2, 11.6), dpi=args.dpi)
        fig.text(0.5, 0.965, f"Ag70 1× Case B | 2× 热源：真实 QSGS 三维体素图 | 检查点 t = {metadata['time_ps']:.2f} ps", ha="center", va="top", fontproperties=font_heading)
        fig.text(0.5, 0.925, f"Ag 液相体积分数 {100*ag_liquid_fraction:.4f}%  |  最大连通液态团簇 {largest_cluster_volume_A3:,.0f} Å³（{largest_voxel_count} voxels）  |  f_l ≥ 0.5 表示液相占优，不表示完全熔化", ha="center", va="top", fontproperties=font_subtitle)

        ax_full = fig.add_axes([0.015, 0.12, 0.53, 0.73], projection="3d")
        add_voxel_surface(ax_full, phase_ids, spacing_nm)
        ax_full.view_init(elev=24, azim=225)
        configure_axes(ax_full, np.asarray([-0.4, -0.4, -0.6]), domain_size_nm + np.asarray([0.4, 0.4, 0.9]), domain_size_nm + np.asarray([0.8, 0.8, 1.5]))
        ax_full.set_title("A  完整样品（100 × 100 × 50 体素）", pad=10, fontproperties=font_panel)

        # Mark the top Gaussian source footprint and direction. The physical
        # source is a smooth Gaussian, so the circle is only a visual radius.
        cx_A, cy_A = identity["source_center_xy_A"]
        cx_nm, cy_nm = cx_A / 10.0, cy_A / 10.0
        theta = np.linspace(0.0, 2.0 * np.pi, 120)
        radius_nm = float(identity["source_gaussian_radius_A"]) / 10.0
        source_ring_z = domain_size_nm[2] + 0.06
        source_ring_xyz = np.column_stack((cx_nm + radius_nm*np.cos(theta), cy_nm + radius_nm*np.sin(theta), np.full_like(theta, source_ring_z)))
        projected = proj3d.proj_transform(source_ring_xyz[:, 0], source_ring_xyz[:, 1], source_ring_xyz[:, 2], ax_full.get_proj())
        projected_display = ax_full.transData.transform(np.column_stack(projected[:2]))
        source_ring_axes = ax_full.transAxes.inverted().transform(projected_display)
        ax_full.add_artist(Line2D(source_ring_axes[:, 0], source_ring_axes[:, 1], transform=ax_full.transAxes, color="#c92323", linewidth=2.4, linestyle="--", zorder=1000))
        center_projected = proj3d.proj_transform(cx_nm, cy_nm, domain_size_nm[2], ax_full.get_proj())
        center_display = ax_full.transData.transform([[center_projected[0], center_projected[1]]])
        source_center_axes = ax_full.transAxes.inverted().transform(center_display)[0]
        ax_full.annotate("", xy=source_center_axes, xycoords="axes fraction", xytext=(0.36, 0.89), textcoords="axes fraction", arrowprops={"arrowstyle": "->", "color": "#c92323", "lw": 1.8}, zorder=1001)
        ax_full.text2D(0.02, 0.97, "顶部：2× 高斯热源（红圈为 R=3 nm 参考范围）", transform=ax_full.transAxes, fontproperties=font_body, color="#a31919", va="top")

        # Show downward sink arrows only at real bottom-face Ag contacts.
        bottom_contacts = np.argwhere(ag_mask[:, :, 0])
        targets = np.asarray([[0.25*phase_ids.shape[0], 0.25*phase_ids.shape[1]], [0.50*phase_ids.shape[0], 0.50*phase_ids.shape[1]], [0.75*phase_ids.shape[0], 0.75*phase_ids.shape[1]]])
        chosen = []
        for target in targets:
            nearest = int(np.argmin(np.sum((bottom_contacts - target) ** 2, axis=1)))
            chosen.append(bottom_contacts[nearest])
        chosen = np.unique(np.asarray(chosen), axis=0)
        for ix, iy in chosen:
            x_nm = (ix + 0.5) * spacing_nm[0]
            y_nm = (iy + 0.5) * spacing_nm[1]
            ax_full.quiver(x_nm, y_nm, 0.05, 0.0, 0.0, -0.55, color="#008c95", linewidth=1.8, arrow_length_ratio=0.5)
        # Outline the physical bottom plane at its visible lower edges; the
        # actual sink contacts remain limited to the native Ag mask at z=0.
        for x_edge in (0.0, domain_size_nm[0]):
            ax_full.plot([x_edge, x_edge], [0.0, domain_size_nm[1]], [0.0, 0.0], color="#008c95", linewidth=1.8, linestyle="--")
        for y_edge in (0.0, domain_size_nm[1]):
            ax_full.plot([0.0, domain_size_nm[0]], [y_edge, y_edge], [0.0, 0.0], color="#008c95", linewidth=1.8, linestyle="--")
        ax_full.text2D(0.02, 0.92, "底面 z=0：300 K 理想热沉，仅耦合实际接触底面的 Ag", transform=ax_full.transAxes, fontproperties=font_body, color="#006b72", va="top")

        ax_zoom = fig.add_axes([0.54, 0.12, 0.445, 0.73], projection="3d")
        add_voxel_surface(ax_zoom, crop_ids, spacing_nm, crop_origin_nm)
        ax_zoom.view_init(elev=26, azim=225)
        lo_nm = crop_origin_nm - np.asarray([0.25, 0.25, 0.25])
        hi_nm = crop_origin_nm + crop_size_nm + np.asarray([0.25, 0.25, 0.35])
        configure_axes(ax_zoom, lo_nm, hi_nm, crop_size_nm + np.asarray([0.5, 0.5, 0.6]))
        ax_zoom.set_title("B  最大连通熔池局部放大（同一检查点）", pad=10, fontproperties=font_panel)
        ax_zoom.text2D(0.02, 0.97, "橙色为 f_l >= 0.5 的 Ag；半透明显示其后的骨架", transform=ax_zoom.transAxes, fontproperties=font_body, color="#7a4015", va="top")

        legend_handles = [Patch(facecolor=color, edgecolor="#39434a", alpha=alpha, label=label) for _, label, color, alpha in PHASES]
        fig.legend(handles=legend_handles, loc="lower center", bbox_to_anchor=(0.5, 0.078), ncol=4, frameon=True, prop=font_legend, columnspacing=1.7, handlelength=1.8)
        fig.text(0.5, 0.039, f"样品尺寸 {domain_size_nm[0]:.2f} × {domain_size_nm[1]:.2f} × {domain_size_nm[2]:.2f} nm³  |  坐标轴单位：nm  |  Ag 底面接触率 {100*bottom_contacts.shape[0]/(phase_ids.shape[0]*phase_ids.shape[1]):.2f}%", ha="center", va="center", fontproperties=font_body)
        fig.savefig(output, dpi=args.dpi, bbox_inches="tight")
        plt.close(fig)

    audit = {
        "image": str(output),
        "checkpoint": str(checkpoint),
        "checkpoint_sha256": sha256_file(checkpoint),
        "checkpoint_time_ps": float(metadata["time_ps"]),
        "geometry": str(geometry_path),
        "geometry_sha256": sha256_file(geometry_path),
        "grid_shape_xyz": list(phase_ids.shape),
        "domain_size_nm": domain_size_nm.tolist(),
        "source_multiplier": 2.0,
        "source_center_nm": [cx_nm, cy_nm],
        "source_gaussian_radius_nm": radius_nm,
        "bottom_sink_K": 300.0,
        "ag_liquid_volume_fraction": ag_liquid_fraction,
        "largest_connected_track_id": largest_track_id,
        "largest_connected_liquid_voxels": largest_voxel_count,
        "largest_connected_liquid_volume_A3": largest_cluster_volume_A3,
        "phase_counts_voxels": {label: int(np.count_nonzero(phase_ids == phase_id)) for phase_id, label, _, _ in PHASES},
        "note": "f_l >= 0.5 is liquid-dominant Ag, not fully melted Ag. Original checkpoint and QSGS geometry are read-only inputs.",
    }
    audit_path = output.with_suffix(".json")
    audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
