"""Render two actual local zoom-ins below the existing Ag70 scaffold RVE.

The selected windows come from lower/upper quantiles of the voxel-resolved
local-inscribed-sphere diameter proxy. They are two regions of one geometry,
not separate micro- and submicron structures.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/ag70_mplconfig")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import ConnectionPatch
import matplotlib.patheffects as path_effects
import numpy as np
from mpl_toolkits.mplot3d import proj3d
import vtk
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray, vtk_to_numpy
from scipy.ndimage import distance_transform_edt, maximum_filter

from render_ag70_model_like_reference import mesh_for_mask, add_mesh


HERE = Path(__file__).resolve().parent
RUN_DIR = HERE.parent
GEOMETRY_PATH = RUN_DIR / "qsgs_morphology_ag70vol.npz"
METADATA_PATH = RUN_DIR / "composition_70vol_geometry.json"
OUTPUT_PATH = HERE / "ag70_scaffold_two_local_zoom_insets.png"
NOTES_PATH = HERE / "SCAFFOLD_ZOOM_INSETS_NOTES.md"

SCAFFOLD_COLOR = "#39AFC6"
BLUE = "#2878B5"
ORANGE = "#E58B32"
INK = "#1C2733"
HALF_WINDOW = (9, 9, 8)  # 19 x 19 x 17 voxels per local rendering


def thickness_peaks(mask: np.ndarray, spacing_nm: np.ndarray):
    """Return peak locations and local inscribed-sphere diameter proxies."""
    nx, ny, nz = mask.shape
    tiled = np.tile(mask, (3, 3, 1))
    edt = distance_transform_edt(tiled, sampling=spacing_nm)
    center = edt[nx:2 * nx, ny:2 * ny, 2:-2]
    inner = mask[:, :, 2:-2]
    maxima = inner & (center >= maximum_filter(center, size=3, mode="nearest") - 1e-12)
    coords = np.argwhere(maxima)
    diameters = 2 * center[maxima]
    coords[:, 2] += 2
    return coords, diameters


def select_region(coords, diameters, mask_shape, quantile):
    """Select an interior peak nearest the requested global thickness quantile."""
    nx, ny, nz = mask_shape
    x, y, z = coords.T
    # Keep the displayed crop away from non-periodic z faces and lateral box edges.
    hx, hy, hz = HALF_WINDOW
    interior = (
        (x >= hx) & (x < nx - hx)
        & (y >= hy) & (y < ny - hy)
        & (z >= hz) & (z < nz - hz)
    )
    target = float(np.percentile(diameters, quantile))
    ids = np.flatnonzero(interior)
    chosen = ids[np.argmin(np.abs(diameters[ids] - target))]
    return coords[chosen].astype(int), float(diameters[chosen]), target


def periodic_crop(mask, center):
    nx, ny, _ = mask.shape
    hx, hy, hz = HALF_WINDOW
    tiled = np.tile(mask, (3, 3, 1))
    x, y, z = (int(v) for v in center)
    x += nx
    y += ny
    return tiled[x - hx:x + hx + 1, y - hy:y + hy + 1,
                 max(0, z - hz):min(mask.shape[2], z + hz + 1)]


def draw_crop(ax, crop, spacing_nm, color, title, diameter):
    crop_lengths = np.asarray(crop.shape, dtype=float) * spacing_nm
    points, faces = mesh_for_mask(crop, spacing_nm, reduction=0.25)
    points = smooth_surface(points, faces)
    add_mesh(ax, points, faces, color, alpha=1.0)
    ax.set_box_aspect(tuple(crop_lengths))
    ax.view_init(elev=25, azim=-55)
    ax.set_xlim(0, crop_lengths[0])
    ax.set_ylim(0, crop_lengths[1])
    ax.set_zlim(0, crop_lengths[2])
    ax.set_xticks([0, round(crop_lengths[0] / 2, 1), round(crop_lengths[0], 1)])
    ax.set_yticks([0, round(crop_lengths[1] / 2, 1), round(crop_lengths[1], 1)])
    ax.set_zticks([0, round(crop_lengths[2] / 2, 1), round(crop_lengths[2], 1)])
    ax.tick_params(labelsize=8, pad=-1)
    ax.set_xlabel("局部 x (nm)", labelpad=0, fontsize=8)
    ax.set_ylabel("局部 y (nm)", labelpad=0, fontsize=8)
    ax.set_zlabel("局部 z (nm)", labelpad=0, fontsize=8)
    ax.set_title(title, loc="left", color=INK, fontsize=12, fontweight="bold", pad=4)
    ax.text2D(
        0.02, 0.94,
        f"局部厚度代理 ≈ {diameter:.2f} nm  ({diameter / np.mean(spacing_nm):.1f} 个体素)",
        transform=ax.transAxes, fontsize=9, color=INK,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.92, edgecolor="#D8DEE5"),
    )


def smooth_surface(points, faces):
    """Smooth only the display mesh; the underlying voxel mask is unchanged."""
    vtk_points = vtk.vtkPoints()
    vtk_points.SetData(numpy_to_vtk(points.astype(np.float64), deep=True, array_type=vtk.VTK_DOUBLE))
    packed_faces = np.column_stack((np.full(len(faces), 3, dtype=np.int64), faces.astype(np.int64))).ravel()
    cells = vtk.vtkCellArray()
    cells.SetCells(len(faces), numpy_to_vtkIdTypeArray(packed_faces, deep=True))
    poly = vtk.vtkPolyData()
    poly.SetPoints(vtk_points)
    poly.SetPolys(cells)

    smoother = vtk.vtkWindowedSincPolyDataFilter()
    smoother.SetInputData(poly)
    smoother.SetNumberOfIterations(12)
    smoother.SetPassBand(0.08)
    smoother.NormalizeCoordinatesOn()
    smoother.NonManifoldSmoothingOn()
    smoother.BoundarySmoothingOff()
    smoother.FeatureEdgeSmoothingOff()
    smoother.Update()
    return vtk_to_numpy(smoother.GetOutput().GetPoints().GetData()).copy()


def annotate_dimension(ax, center_nm, length_nm, label, color):
    """Draw a projected double arrow through the inscribed-sphere center."""
    center_nm = np.asarray(center_nm, dtype=float)
    candidates = []
    for axis in range(3):
        start = center_nm.copy()
        end = center_nm.copy()
        start[axis] -= length_nm / 2
        end[axis] += length_nm / 2
        p0 = proj3d.proj_transform(*start, ax.get_proj())[:2]
        p1 = proj3d.proj_transform(*end, ax.get_proj())[:2]
        a0 = ax.transAxes.inverted().transform(ax.transData.transform(p0))
        a1 = ax.transAxes.inverted().transform(ax.transData.transform(p1))
        candidates.append((np.linalg.norm(a1 - a0), a0, a1))
    _, a0, a1 = max(candidates, key=lambda item: item[0])
    ax.annotate(
        "", xy=a1, xytext=a0, xycoords=ax.transAxes, textcoords=ax.transAxes,
        arrowprops=dict(
            arrowstyle="<->", color=color, linewidth=2.8,
            shrinkA=0, shrinkB=0, mutation_scale=17,
            path_effects=[path_effects.Stroke(linewidth=4.5, foreground="#18323D"),
                          path_effects.Normal()],
        ), zorder=50,
    )
    midpoint = (a0 + a1) / 2
    ax.text2D(
        midpoint[0], midpoint[1] + 0.055, label, transform=ax.transAxes,
        color=color, fontsize=10, fontweight="bold", ha="center", va="bottom",
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.9, edgecolor=color),
        zorder=31,
    )


def draw_box(ax, lo, hi, color):
    x0, y0, z0 = lo
    x1, y1, z1 = hi
    corners = np.array([
        [x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],
        [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1],
    ])
    edges = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
             (0, 4), (1, 5), (2, 6), (3, 7)]
    for a, b in edges:
        ax.plot(*corners[[a, b]].T, color=color, linestyle="--", linewidth=1.5)


def main():
    with np.load(GEOMETRY_PATH) as data:
        scaffold = data["skeleton"].astype(bool)
    metadata = json.loads(METADATA_PATH.read_text(encoding="utf-8"))
    spacing_nm = np.asarray(metadata["voxel_size_A"], dtype=float) / 10.0
    lengths_nm = np.asarray(metadata["domain_size_A"], dtype=float) / 10.0

    coords, diameters = thickness_peaks(scaffold, spacing_nm)
    thin_center, thin_d, thin_target = select_region(coords, diameters, scaffold.shape, 10)
    thick_center, thick_d, thick_target = select_region(coords, diameters, scaffold.shape, 90)
    thin_crop = periodic_crop(scaffold, thin_center)
    thick_crop = periodic_crop(scaffold, thick_center)

    font_file = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
    if Path(font_file).exists():
        font_manager.fontManager.addfont(font_file)
        cjk_family = font_manager.FontProperties(fname=font_file).get_name()
    else:
        cjk_family = "DejaVu Sans"
    plt.rcParams.update({
        "font.family": [cjk_family, "DejaVu Sans"],
        "axes.unicode_minus": False,
        "figure.facecolor": "white",
        "font.size": 10,
    })

    fig = plt.figure(figsize=(13, 10), constrained_layout=False)
    gs = fig.add_gridspec(2, 2, left=0.07, right=0.97, top=0.88, bottom=0.09,
                          height_ratios=(1.05, 0.95), wspace=0.12, hspace=0.23)
    ax_main = fig.add_subplot(gs[0, :], projection="3d")
    points, faces = mesh_for_mask(scaffold, spacing_nm, reduction=0.85)
    points = smooth_surface(points, faces)
    add_mesh(ax_main, points, faces, SCAFFOLD_COLOR, alpha=1.0)
    ax_main.set_xlim(0, lengths_nm[0])
    ax_main.set_ylim(0, lengths_nm[1])
    ax_main.set_zlim(0, lengths_nm[2])
    ax_main.set_box_aspect(tuple(lengths_nm))
    ax_main.view_init(elev=25, azim=-55)
    ax_main.set_xlabel("x (nm)")
    ax_main.set_ylabel("y (nm)")
    ax_main.set_zlabel("z (nm)")
    ax_main.set_title("Ti₃SiC₂ 连续骨架：当前 Ag70 QSGS 计算域", loc="left",
                      color=INK, fontsize=14, fontweight="bold", pad=8)

    hx, hy, hz = HALF_WINDOW
    region_info = []
    for center, color in ((thin_center, BLUE), (thick_center, ORANGE)):
        x, y, z = center
        lo = ((x - hx) * spacing_nm[0], (y - hy) * spacing_nm[1], (z - hz) * spacing_nm[2])
        hi = ((x + hx + 1) * spacing_nm[0], (y + hy + 1) * spacing_nm[1],
              (min(scaffold.shape[2], z + hz + 1)) * spacing_nm[2])
        draw_box(ax_main, lo, hi, color)
        region_info.append((lo, hi))
    ax_main.text2D(
        0.015, 0.02,
        f"计算域 {lengths_nm[0]:.2f} × {lengths_nm[1]:.2f} × {lengths_nm[2]:.2f} nm；"
        f"体素 {spacing_nm[0]:.3f} × {spacing_nm[1]:.3f} × {spacing_nm[2]:.3f} nm",
        transform=ax_main.transAxes, fontsize=9, color=INK,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.9, edgecolor="#D8DEE5"),
    )

    ax_thin = fig.add_subplot(gs[1, 0], projection="3d")
    ax_thick = fig.add_subplot(gs[1, 1], projection="3d")
    draw_crop(ax_thin, thin_crop, spacing_nm, SCAFFOLD_COLOR, "(a) 较细分支局部放大（P10）", thin_d)
    draw_crop(ax_thick, thick_crop, spacing_nm, SCAFFOLD_COLOR, "(b) 较粗分支/交汇局部放大（P90）", thick_d)

    fig.suptitle("骨架局部放大图：两个窗口均来自同一纳米 QSGS 结构",
                 fontsize=17, fontweight="bold", color=INK, y=0.965)
    fig.text(0.07, 0.035,
             "说明：虚线框对应上方骨架中的实际取样窗口；小图双箭头标局部最大内切球直径代理（用于表示骨架局部厚度），"
             "P10/P90 为该代理的分位数，"
             "不是实验测得的分支直径。当前体素约 0.30 nm，最细分支只有约 2 个体素，尺寸受网格分辨率限制。",
             fontsize=9, color="#4B5563", ha="left", va="center", wrap=True)

    fig.canvas.draw()
    hx, hy, hz = HALF_WINDOW
    for center, diameter, ax_zoom, color in (
        (thin_center, thin_d, ax_thin, BLUE),
        (thick_center, thick_d, ax_thick, ORANGE),
    ):
        cx, cy, cz = center
        center_local_nm = np.array([
            (hx + 0.5) * spacing_nm[0],
            (hy + 0.5) * spacing_nm[1],
            (hz + 0.5) * spacing_nm[2],
        ])
        annotate_dimension(ax_zoom, center_local_nm, diameter, f"≈ {diameter:.2f} nm", color)

    # Dashed leaders associate the highlighted source windows with the two zoom panels.
    for center, ax_zoom, color in ((thin_center, ax_thin, BLUE), (thick_center, ax_thick, ORANGE)):
        xp, yp, _ = proj3d.proj_transform(*(center * spacing_nm), ax_main.get_proj())
        px, py = ax_main.transAxes.inverted().transform(ax_main.transData.transform((xp, yp)))
        leader = ConnectionPatch(
            xyA=(px, py), coordsA=ax_main.transAxes,
            xyB=(0.5, 1.0), coordsB=ax_zoom.transAxes,
            color=color, linestyle=(0, (4, 3)), linewidth=1.4,
            clip_on=False, zorder=10,
        )
        fig.add_artist(leader)

    fig.savefig(OUTPUT_PATH, dpi=220, bbox_inches="tight")
    plt.close(fig)

    notes = f"""# Ag70 scaffold zoom-in panels

Both inset windows are extracted from the same original `{GEOMETRY_PATH.name}` QSGS geometry.
They are not separate micro- and submicron structures.

- Geometry grid: `{scaffold.shape[0]} × {scaffold.shape[1]} × {scaffold.shape[2]}`.
- Domain: `{lengths_nm[0]:.4f} × {lengths_nm[1]:.4f} × {lengths_nm[2]:.4f} nm`.
- Voxel spacing: `{spacing_nm[0]:.5f} × {spacing_nm[1]:.5f} × {spacing_nm[2]:.5f} nm`.
- Thin-region center voxel (x,y,z): `{tuple(int(v) for v in thin_center)}`; selected local diameter proxy `{thin_d:.4f} nm` (P10 target `{thin_target:.4f} nm`).
- Thick/junction-region center voxel (x,y,z): `{tuple(int(v) for v in thick_center)}`; selected local diameter proxy `{thick_d:.4f} nm` (P90 target `{thick_target:.4f} nm`).
- Each crop is `{thin_crop.shape[0]} × {thin_crop.shape[1]} × {thin_crop.shape[2]}` voxels; x/y crops wrap periodically.

The inset double arrows show twice the Euclidean distance-transform value at a local maximum inside the scaffold. This is a voxel-based proxy for local scaffold thickness, not a direct measurement of branch diameter or pore spacing. At approximately 0.30 nm per voxel, the lower-tail feature is only about two voxels wide and is resolution-limited. Mesh decimation and smoothing are applied only to the rendering; all reported geometry statistics come from the original voxel mask.
"""
    NOTES_PATH.write_text(notes, encoding="utf-8")
    print(f"Wrote {OUTPUT_PATH}")
    print(f"P10 zoom center={tuple(thin_center)} diameter_proxy={thin_d:.4f} nm")
    print(f"P90 zoom center={tuple(thick_center)} diameter_proxy={thick_d:.4f} nm")


if __name__ == "__main__":
    main()
