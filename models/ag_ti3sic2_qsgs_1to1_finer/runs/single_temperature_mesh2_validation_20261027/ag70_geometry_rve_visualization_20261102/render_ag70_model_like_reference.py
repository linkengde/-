"""Render the actual Ag70 QSGS cell in an RVE-style explanatory figure.

This is a visualization of the existing voxel mask, not a regenerated geometry
and not a micron-scale RVE. The local thickness is a voxel-resolved inscribed-
sphere proxy, not a measured physical branch diameter.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/ag70_mplconfig")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import vtk
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from matplotlib.ticker import LogLocator, NullFormatter
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.colors import LightSource
from scipy.ndimage import distance_transform_edt, label, maximum_filter
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy


HERE = Path(__file__).resolve().parent
RUN_DIR = HERE.parent
GEOMETRY_PATH = RUN_DIR / "qsgs_morphology_ag70vol.npz"
METADATA_PATH = RUN_DIR / "composition_70vol_geometry.json"
OUTPUT_PATH = HERE / "ag70_current_qsgs_cell_rve_style.png"
NOTES_PATH = HERE / "GEOMETRY_VISUALIZATION_NOTES.md"

AG_COLOR = "#E9E4D9"
SCAFFOLD_COLOR = "#607891"
BLUE = "#1D4E89"
INK = "#1C2733"
LIGHTSOURCE = LightSource(azdeg=315, altdeg=40)


def mesh_for_mask(mask: np.ndarray, spacing_nm: np.ndarray, reduction: float = 0.95):
    """Extract an isosurface at half-voxel boundaries and decimate for plotting."""
    padded = np.pad(mask.astype(np.uint8), 1)
    image = vtk.vtkImageData()
    image.SetDimensions(*padded.shape)
    image.SetSpacing(*spacing_nm.tolist())
    image.SetOrigin(*(-0.5 * spacing_nm).tolist())
    scalars = numpy_to_vtk(padded.ravel(order="F"), deep=True, array_type=vtk.VTK_UNSIGNED_CHAR)
    scalars.SetName("phase_mask")
    image.GetPointData().SetScalars(scalars)

    surface = vtk.vtkFlyingEdges3D()
    surface.SetInputData(image)
    surface.SetValue(0, 0.5)
    surface.ComputeNormalsOff()
    surface.ComputeGradientsOff()
    surface.Update()

    # vtkDecimatePro barely reduced this highly branched voxel surface while
    # topology preservation was enabled. Quadric decimation is used strictly
    # for display geometry; phase fractions and thickness statistics retain
    # the original full-resolution mask.
    decimate = vtk.vtkQuadricDecimation()
    decimate.SetInputConnection(surface.GetOutputPort())
    decimate.SetTargetReduction(reduction)
    decimate.Update()
    poly = decimate.GetOutput()

    points = vtk_to_numpy(poly.GetPoints().GetData()).copy()
    packed = vtk_to_numpy(poly.GetPolys().GetData())
    faces = packed.reshape(-1, 4)[:, 1:].copy()
    return points, faces


def add_mesh(ax, points, faces, color, alpha=1.0):
    collection = Poly3DCollection(
        points[faces],
        facecolors=color,
        edgecolor="none",
        linewidth=0,
        alpha=alpha,
        shade=True,
        lightsource=LIGHTSOURCE,
        zsort="average",
    )
    ax.add_collection3d(collection)


def set_cube_axes(ax, lengths_nm, *, show_ticks=True):
    lx, ly, lz = lengths_nm
    ax.set_xlim(0, lx)
    ax.set_ylim(0, ly)
    ax.set_zlim(0, lz)
    ax.set_box_aspect((lx, ly, lz))
    ax.view_init(elev=24, azim=-56)
    ax.set_xlabel("x (nm)", labelpad=1)
    ax.set_ylabel("y (nm)", labelpad=1)
    ax.set_zlabel("z (nm)", labelpad=1)
    if show_ticks:
        ax.set_xticks([0, round(lx / 2, 1), round(lx, 1)])
        ax.set_yticks([0, round(ly / 2, 1), round(ly, 1)])
        ax.set_zticks([0, round(lz / 2, 1), round(lz, 1)])
    ax.tick_params(labelsize=8, pad=0)
    ax.xaxis.pane.set_facecolor((0.96, 0.97, 0.98, 0.5))
    ax.yaxis.pane.set_facecolor((0.96, 0.97, 0.98, 0.5))
    ax.zaxis.pane.set_facecolor((0.96, 0.97, 0.98, 0.5))
    ax.grid(True, linewidth=0.35, alpha=0.35)


def local_thickness_proxy(mask, spacing_nm):
    """Use local maxima of periodic-xy EDT as approximate local diameters."""
    tiled = np.tile(mask, (3, 3, 1))
    edt = distance_transform_edt(tiled, sampling=spacing_nm)
    nx, ny, nz = mask.shape
    center = edt[nx:2 * nx, ny:2 * ny, 2:-2]
    inner = mask[:, :, 2:-2]
    peaks = inner & (center >= maximum_filter(center, size=3, mode="nearest") - 1e-12)
    components, count = label(peaks, structure=np.ones((3, 3, 3), dtype=bool))
    radii = np.fromiter(
        (center[components == i].max() for i in range(1, count + 1)),
        dtype=float,
        count=count,
    )
    return 2 * radii


def main():
    with np.load(GEOMETRY_PATH) as data:
        scaffold = data["skeleton"].astype(bool)
    metadata = json.loads(METADATA_PATH.read_text(encoding="utf-8"))
    sizes_A = np.asarray(metadata["domain_size_A"], dtype=float)
    spacing_A = np.asarray(metadata["voxel_size_A"], dtype=float)
    lengths_nm = sizes_A / 10.0
    spacing_nm = spacing_A / 10.0
    ag = ~scaffold
    thickness_nm = local_thickness_proxy(scaffold, spacing_nm)
    p10, p50, p90, p95, tmax = np.percentile(thickness_nm, [10, 50, 90, 95, 100])

    scaffold_points, scaffold_faces = mesh_for_mask(scaffold, spacing_nm)
    ag_points, ag_faces = mesh_for_mask(ag, spacing_nm)

    cjk_font = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
    if Path(cjk_font).exists():
        font_manager.fontManager.addfont(cjk_font)
        cjk_family = font_manager.FontProperties(fname=cjk_font).get_name()
    else:
        cjk_family = "DejaVu Sans"
    plt.rcParams.update({
        "font.family": [cjk_family, "DejaVu Sans"],
        "axes.unicode_minus": False,
        "figure.facecolor": "white",
        "axes.titleweight": "bold",
        "font.size": 10,
    })
    fig = plt.figure(figsize=(16, 11), constrained_layout=False)
    grid = fig.add_gridspec(2, 2, left=0.055, right=0.97, top=0.90, bottom=0.12,
                            wspace=0.16, hspace=0.25)

    ax_a = fig.add_subplot(grid[0, 0], projection="3d")
    # Light Ag shell and blue scaffold surface; the mid-plane and isolated
    # scaffold panels provide unambiguous phase identification as well.
    add_mesh(ax_a, ag_points, ag_faces, AG_COLOR, alpha=0.28)
    add_mesh(ax_a, scaffold_points, scaffold_faces, SCAFFOLD_COLOR, alpha=0.96)
    set_cube_axes(ax_a, lengths_nm)
    ax_a.set_title("(a) 当前计算域：Ag / Ti₃SiC₂ 双连续结构", loc="left", color=BLUE, pad=9)
    ax_a.legend(handles=[
        Patch(facecolor=AG_COLOR, edgecolor="#8C887E", label=f"Ag  {ag.mean() * 100:.4f} vol.%"),
        Patch(facecolor=SCAFFOLD_COLOR, edgecolor="#34495E", label=f"Ti₃SiC₂  {scaffold.mean() * 100:.4f} vol.%"),
    ], loc="upper left", fontsize=8, frameon=True, framealpha=0.9)
    ax_a.text2D(0.03, 0.02, f"域尺寸 {lengths_nm[0]:.2f} × {lengths_nm[1]:.2f} × {lengths_nm[2]:.2f} nm",
                transform=ax_a.transAxes, fontsize=8, color=INK)

    ax_b = fig.add_subplot(grid[0, 1], projection="3d")
    add_mesh(ax_b, scaffold_points, scaffold_faces, SCAFFOLD_COLOR, alpha=1.0)
    set_cube_axes(ax_b, lengths_nm, show_ticks=False)
    ax_b.set_title("(b) Ti₃SiC₂ 连续骨架提取", loc="left", color=BLUE, pad=9)
    ax_b.text2D(0.02, 0.02,
                f"局部厚度代理：中位数 {p50:.2f} nm；P10–P90 {p10:.2f}–{p90:.2f} nm",
                transform=ax_b.transAxes, fontsize=8, color=INK)

    ax_c = fig.add_subplot(grid[1, 0])
    mid_z = scaffold.shape[2] // 2
    phase = scaffold[:, :, mid_z].T.astype(np.uint8)
    cmap = ListedColormap([AG_COLOR, SCAFFOLD_COLOR])
    ax_c.imshow(phase, origin="lower", interpolation="nearest", cmap=cmap,
                extent=(0, lengths_nm[0], 0, lengths_nm[1]), aspect="equal", vmin=0, vmax=1)
    ax_c.set_xlim(0, lengths_nm[0])
    ax_c.set_ylim(0, lengths_nm[1])
    ax_c.set_xlabel("x (nm)")
    ax_c.set_ylabel("y (nm)")
    ax_c.set_title(f"(c) 中央截面（z = {lengths_nm[2] / 2:.2f} nm）", loc="left", color=BLUE)
    ax_c.legend(handles=[
        Patch(facecolor=AG_COLOR, edgecolor="#8C887E", label="Ag"),
        Patch(facecolor=SCAFFOLD_COLOR, edgecolor="#34495E", label="Ti₃SiC₂"),
    ], loc="upper right", fontsize=8)
    bar_x, bar_y = 2.0, 2.0
    bar_len = 5.0
    ax_c.plot([bar_x, bar_x + bar_len], [bar_y, bar_y], color="black", lw=3)
    ax_c.text(bar_x + bar_len / 2, bar_y + 0.7, "5 nm", ha="center", fontsize=8)

    ax_d = fig.add_subplot(grid[1, 1])
    ax_d.set_xscale("log")
    ax_d.set_xlim(0.3, 5000)
    ax_d.set_ylim(-0.7, 2.7)
    ax_d.set_yticks([0, 1, 2])
    ax_d.set_yticklabels([
        "当前几何：P10–P90 厚度代理",
        "图示超细级：约 0.2 μm",
        "图示微米级：约 2 μm",
    ])
    ax_d.hlines(0, p10, p90, color=SCAFFOLD_COLOR, linewidth=8, alpha=0.8)
    ax_d.scatter([p50], [0], color="#D1495B", s=65, zorder=3, label=f"当前中位数 {p50:.2f} nm")
    ax_d.scatter([200], [1], color="#4C78A8", s=60, zorder=3)
    ax_d.scatter([2000], [2], color="#F28E2B", s=60, zorder=3)
    ax_d.axvline(0.6, color="#999999", linestyle="--", linewidth=0.8)
    ax_d.axvline(2.0, color="#999999", linestyle="--", linewidth=0.8)
    ax_d.xaxis.set_major_locator(LogLocator(base=10, numticks=7))
    ax_d.xaxis.set_minor_formatter(NullFormatter())
    ax_d.grid(True, which="both", axis="x", alpha=0.25)
    ax_d.set_xlabel("特征尺寸 / 局部厚度代理（nm，对数刻度）")
    ax_d.set_title("(d) 与示意图尺度对照", loc="left", color=BLUE)
    ax_d.text(0.02, 0.02,
              f"当前网格 {scaffold.shape[0]}×{scaffold.shape[1]}×{scaffold.shape[2]}；\n"
              f"体素 {spacing_nm[0]:.3f}×{spacing_nm[1]:.3f}×{spacing_nm[2]:.3f} nm；\n"
              f"局部最大值 {tmax:.2f} nm（可能是交汇处）",
              transform=ax_d.transAxes, fontsize=8, color=INK, va="bottom",
              bbox=dict(boxstyle="round,pad=0.35", facecolor="white", alpha=0.88, edgecolor="#D6DCE2"))

    fig.suptitle("Ag–Ti₃SiC₂ 70:30 现有 QSGS 几何：RVE 风格展示（实际为纳米计算域）",
                 fontsize=17, fontweight="bold", color=INK, y=0.965)
    fig.text(0.055, 0.055,
             "注：骨架厚度采用周期 x/y 距离变换的局部最大内切球直径代理；体素约 0.30 nm，最细特征仅数个体素，"
             "不能视为实验测得的分支尺寸。图示 0.2/2 μm 为参考图尺度，不是本模型尺寸。",
             fontsize=9, color="#4B5563", ha="left", va="center", wrap=True)
    fig.savefig(OUTPUT_PATH, dpi=220, bbox_inches="tight")
    plt.close(fig)

    notes = f"""# Ag70 current geometry visualization

This figure renders the existing `qsgs_morphology_ag70vol.npz`; it does not regenerate or edit the QSGS geometry.

- Grid: `{scaffold.shape[0]} × {scaffold.shape[1]} × {scaffold.shape[2]}`.
- Domain: `{lengths_nm[0]:.4f} × {lengths_nm[1]:.4f} × {lengths_nm[2]:.4f} nm`.
- Voxel spacing: `{spacing_nm[0]:.5f} × {spacing_nm[1]:.5f} × {spacing_nm[2]:.5f} nm`.
- Volume fractions: Ag `{ag.mean() * 100:.4f}%`; Ti₃SiC₂ `{scaffold.mean() * 100:.4f}%`.
- Local inscribed-sphere diameter proxy on periodic-x/y EDT maxima: P10 `{p10:.3f} nm`, median `{p50:.3f} nm`, P90 `{p90:.3f} nm`, P95 `{p95:.3f} nm`, max `{tmax:.3f} nm`.
- The EDT maxima are an approximate local-thickness statistic, not a uniquely segmented branch diameter or measured physical dimension. With approximately 0.30 nm voxels, a 0.60 nm estimate is only about two voxels and is resolution-limited.
- The 0.2 µm and 2 µm markers are reference scales from the user's supplied conceptual figure, not dimensions in this current nanoscale model.

The x/y periodicity is accounted for by tiling the mask 3×3 before the distance transform. The top and bottom two voxel planes are excluded from thickness statistics. The mesh in the 3-D panels is decimated only for display; all calculations use the original mask.
"""
    NOTES_PATH.write_text(notes, encoding="utf-8")
    print(f"Wrote {OUTPUT_PATH}")
    print(f"Thickness proxy nm: P10={p10:.3f}, median={p50:.3f}, P90={p90:.3f}, P95={p95:.3f}, max={tmax:.3f}")


if __name__ == "__main__":
    main()
