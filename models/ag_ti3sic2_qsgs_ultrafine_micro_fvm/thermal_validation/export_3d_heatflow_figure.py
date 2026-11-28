#!/usr/bin/env python3
"""用VTK提取相界面，再以Matplotlib绘制真实求解场的三维通路图。"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy


HERE = Path(__file__).resolve().parent
CASE_ID = "performance_n150_seed1105_ultrafine_500nm_R0e+00"
NPZ = HERE / "results" / "cases" / f"{CASE_ID}_axis0_fields.npz"
OUT = HERE / "figures" / "three_dimensional_heatflow_pathways.png"
N = 150
EDGE_UM = 10.0
SPACING_UM = EDGE_UM / N
COARSE_STEP = 3


def extract_surface(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """以1/3空间预览密度抽取骨架—Ag界面，只作三维显示。"""
    coarse = np.ascontiguousarray(mask[::COARSE_STEP, ::COARSE_STEP, ::COARSE_STEP], dtype=np.uint8)
    image = vtk.vtkImageData()
    image.SetDimensions(*coarse.shape)
    image.SetOrigin(*(0.5 * SPACING_UM,) * 3)
    image.SetSpacing(*(COARSE_STEP * SPACING_UM,) * 3)
    scalars = numpy_to_vtk(np.transpose(coarse, (2, 1, 0)).ravel(order="C"), deep=True,
                           array_type=vtk.VTK_UNSIGNED_CHAR)
    scalars.SetName("skeleton_id")
    image.GetPointData().SetScalars(scalars)

    surface = vtk.vtkFlyingEdges3D()
    surface.SetInputData(image)
    surface.SetValue(0, 0.5)
    surface.ComputeNormalsOff()
    surface.Update()
    simplify = vtk.vtkDecimatePro()
    simplify.SetInputConnection(surface.GetOutputPort())
    simplify.SetTargetReduction(0.72)
    simplify.PreserveTopologyOn()
    simplify.Update()
    triangles = vtk.vtkTriangleFilter()
    triangles.SetInputConnection(simplify.GetOutputPort())
    triangles.Update()
    poly = triangles.GetOutput()
    points = vtk_to_numpy(poly.GetPoints().GetData()).astype(np.float32, copy=True)
    packed = vtk_to_numpy(poly.GetPolys().GetData()).reshape(-1, 4)
    faces = packed[:, 1:].astype(np.int32, copy=True)
    return points, faces


def axes_setup(ax: plt.Axes, title: str) -> None:
    ax.set_title(title, fontsize=10, pad=2)
    ax.set_xlabel("x (μm)", labelpad=-2)
    ax.set_ylabel("y (μm)", labelpad=-2)
    ax.set_zlabel("z (μm)", labelpad=-1)
    ax.set_xlim(0, EDGE_UM)
    ax.set_ylim(0, EDGE_UM)
    ax.set_zlim(0, EDGE_UM)
    ax.set_xticks((0, 5, 10))
    ax.set_yticks((0, 5, 10))
    ax.set_zticks((0, 5, 10))
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=24, azim=-55)
    ax.grid(False)
    ax.xaxis.pane.set_facecolor((0.97, 0.97, 0.97, 1.0))
    ax.yaxis.pane.set_facecolor((0.97, 0.97, 0.97, 1.0))
    ax.zaxis.pane.set_facecolor((0.97, 0.97, 0.97, 1.0))


def main() -> None:
    with np.load(NPZ, allow_pickle=False) as archive:
        data = {key: np.asarray(archive[key]) for key in archive.files}
    skeleton = data["material_id"].astype(bool)
    points, faces = extract_surface(skeleton)
    q = np.sqrt(sum(data[f"heat_flux_{axis}_W_m2"].astype(np.float64) ** 2 for axis in range(3)))
    ag = ~skeleton
    cutoff = float(np.quantile(q[ag], 0.99))
    selected = np.flatnonzero((ag & (q >= cutoff)).ravel())
    # 限制点数以保持出版图可读；按通量排序保留最强的Ag单元。
    rank = np.argsort(q.ravel()[selected])[-9000:]
    flat = selected[rank]
    ijk = np.column_stack(np.unravel_index(flat, q.shape))
    xyz_um = (ijk + 0.5) * SPACING_UM
    values = q.ravel()[flat]

    fig = plt.figure(figsize=(10.2, 5.1), constrained_layout=True)
    for panel, overlay in enumerate((False, True), start=1):
        ax = fig.add_subplot(1, 2, panel, projection="3d")
        mesh = Poly3DCollection(points[faces], facecolor="#b9b0a7", edgecolor="none", alpha=0.20,
                                linewidth=0.0, zsort="average")
        ax.add_collection3d(mesh)
        if overlay:
            cloud = ax.scatter(xyz_um[:, 0], xyz_um[:, 1], xyz_um[:, 2], c=values,
                               cmap="inferno", s=2.8, alpha=0.9, linewidths=0,
                               depthshade=False, rasterized=True)
        axes_setup(ax, "Ti₃SiC₂ skeleton surface" if not overlay else "Skeleton + high-flux Ag cells")
    cb = fig.colorbar(cloud, ax=fig.axes[1], shrink=0.64, pad=0.06, aspect=24)
    cb.set_label("|q| (W/m²); Ag top 1% by cell flux")
    fig.suptitle("Computed 3D heat-flow pathways · 10 μm RVE · seed 1105 · 150³ · X-gradient", fontsize=12)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"saved={OUT} surface_triangles={len(faces)} ag_high_flux_cells={len(flat)} flux_cutoff_W_m2={cutoff:.6g}")


if __name__ == "__main__":
    main()
