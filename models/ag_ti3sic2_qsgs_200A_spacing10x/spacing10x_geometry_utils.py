"""固定尺寸枝状相场的体素连通、表面网格和测量工具。"""
from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np
from matplotlib.colors import to_rgb
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy import ndimage as ndi
from skimage.measure import marching_cubes
from skimage.morphology import skeletonize
import vtk
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray

VOXEL_A = 2.5
GRID = (80, 80, 80)
BOX_GRID = (80, 80, 120)
DISPLAY_ISOLEVEL = 0.55
DISPLAY_SIGMA_VOXELS = 0.7
SOURCE_GAP_A = 9.3719


def periodic_labels(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """对 x/y 周期、z 开放的六邻接体素区域编号。"""
    labels, count = ndi.label(mask, structure=ndi.generate_binary_structure(3, 1))
    parent = np.arange(count + 1, dtype=np.int64)

    def root(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = int(parent[value])
        return value

    for axis in (0, 1):
        low = np.take(labels, 0, axis=axis).ravel()
        high = np.take(labels, -1, axis=axis).ravel()
        for a, b in np.unique(np.column_stack((low, high)), axis=0):
            if a and b:
                ra, rb = root(int(a)), root(int(b))
                if ra != rb:
                    parent[rb] = ra

    roots = np.fromiter((root(i) for i in range(count + 1)), dtype=np.int64)
    labels = roots[labels]
    sizes = np.bincount(labels.ravel(), minlength=count + 1)
    sizes[0] = 0
    return labels, sizes


def sample_surface_field(field: np.ndarray, positions_a: np.ndarray) -> np.ndarray:
    """三线性采样周期 x/y、开放 z 的相场。"""
    indices = np.asarray(positions_a, dtype=float) / VOXEL_A - 0.5
    indices[0] %= GRID[0]
    indices[1] %= GRID[1]
    indices[2] = np.clip(indices[2], 0.0, GRID[2] - 1.0)
    return ndi.map_coordinates(field, indices, order=1, mode="nearest")


def find_display_gap(mask: np.ndarray) -> dict:
    """在平滑的骨架显示表面上测量接近原枝间距十倍的 Ag 通道。"""
    ag = ~mask
    distance, nearest = ndi.distance_transform_edt(
        ag, sampling=(VOXEL_A,) * 3, return_indices=True
    )
    ag_axis = skeletonize(ag, method="lee")
    z_index = np.indices(ag_axis.shape)[2]
    candidates = np.argwhere(ag_axis & (z_index > 15) & (z_index < 65))
    target_a = SOURCE_GAP_A * 10.0
    candidates = sorted(candidates, key=lambda q: abs(2.0 * distance[tuple(q)] - target_a))
    display_field = ndi.gaussian_filter(
        mask.astype(np.float32), sigma=DISPLAY_SIGMA_VOXELS,
        mode=("wrap", "wrap", "nearest"),
    ) - DISPLAY_ISOLEVEL
    steps = np.arange(0.0, 140.0, 0.5)

    best = None
    best_error = float("inf")
    for q_int in candidates:
        q = q_int.astype(float)
        closest = np.array([nearest[axis][tuple(q_int)] for axis in range(3)], dtype=float)
        direction = q - closest
        norm = float(np.linalg.norm(direction))
        if norm == 0:
            continue
        direction /= norm
        q_a = (q + 0.5) * VOXEL_A
        surface_points = []
        valid = True
        for sign in (-1.0, 1.0):
            ray = q_a[:, None] + sign * direction[:, None] * steps[None, :]
            values = sample_surface_field(display_field, ray)
            crossings = np.flatnonzero(values[1:] * values[:-1] <= 0.0)
            if not len(crossings):
                valid = False
                break
            i = int(crossings[0])
            denominator = abs(values[i]) + abs(values[i + 1])
            fraction = abs(values[i]) / denominator if denominator else 0.0
            travel = float(steps[i] + fraction * (steps[i + 1] - steps[i]))
            surface_points.append(q_a + sign * direction * travel)
        if not valid or any(p[2] < 10.0 or p[2] > 190.0 for p in surface_points):
            continue
        gap_a = float(np.linalg.norm(surface_points[0] - surface_points[1]))
        error = abs(gap_a - target_a)
        if error > 8.0 or error >= best_error:
            continue
        best_error = error
        best = {
            "center_index_xyz": q_int.astype(int).tolist(),
            "center_xyz_A": q_a.tolist(),
            "direction_xyz": direction.tolist(),
            "surface_points_xyz_A": [p.tolist() for p in surface_points],
            "gap_A": gap_a,
            "gap_nm": gap_a / 10.0,
            "target_from_original_gap_A": target_a,
            "ratio_to_original_gap": gap_a / SOURCE_GAP_A,
        }
    if best is None:
        raise RuntimeError("没有找到位于材料内部、可测量的枝间通道")
    return best


def write_vtk_field(path: Path, phase: np.ndarray, name: str) -> None:
    """写入 VTK STRUCTURED_POINTS 体素场，坐标单位为 Å。"""
    image = vtk.vtkImageData()
    image.SetDimensions(*(int(n) + 1 for n in phase.shape))
    image.SetSpacing(VOXEL_A, VOXEL_A, VOXEL_A)
    image.SetOrigin(0.0, 0.0, 0.0)
    data = np.asfortranarray(phase, dtype=np.uint8).ravel(order="F")
    scalar = numpy_to_vtk(data, deep=True, array_type=vtk.VTK_UNSIGNED_CHAR)
    scalar.SetName(name)
    image.GetCellData().SetScalars(scalar)
    writer = vtk.vtkStructuredPointsWriter()
    writer.SetFileName(str(path))
    writer.SetFileTypeToASCII()
    writer.SetInputData(image)
    if writer.Write() != 1:
        raise RuntimeError(f"VTK 写入失败：{path}")


def mesh_from_field(field: np.ndarray, origin_a: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """从平滑体素场提取三角面，返回 Å 坐标。"""
    padded = np.pad(field, 1, mode="constant", constant_values=0.0)
    vertices, faces, _, _ = marching_cubes(
        padded, level=DISPLAY_ISOLEVEL, spacing=(VOXEL_A,) * 3
    )
    vertices += np.asarray(origin_a, dtype=float) - 0.5 * VOXEL_A
    return vertices.astype(np.float32), faces.astype(np.int32)


def mesh_polydata(vertices: np.ndarray, faces: np.ndarray) -> vtk.vtkPolyData:
    """把 NumPy 三角面转为平滑的 VTK 表面网格。"""
    points = vtk.vtkPoints()
    points.SetData(numpy_to_vtk(np.ascontiguousarray(vertices), deep=True))
    cells = np.column_stack((np.full(len(faces), 3, dtype=np.int64), faces)).ravel()
    polygons = vtk.vtkCellArray()
    polygons.SetCells(len(faces), numpy_to_vtkIdTypeArray(cells, deep=True))
    poly = vtk.vtkPolyData()
    poly.SetPoints(points)
    poly.SetPolys(polygons)
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(poly)
    normals.ConsistencyOn()
    normals.AutoOrientNormalsOn()
    normals.SplittingOff()
    normals.Update()
    smooth = vtk.vtkWindowedSincPolyDataFilter()
    smooth.SetInputConnection(normals.GetOutputPort())
    smooth.SetNumberOfIterations(8)
    smooth.SetPassBand(0.12)
    smooth.BoundarySmoothingOn()
    smooth.FeatureEdgeSmoothingOff()
    smooth.NonManifoldSmoothingOn()
    smooth.NormalizeCoordinatesOn()
    smooth.Update()
    result = vtk.vtkPolyData()
    result.DeepCopy(smooth.GetOutput())
    return result


def write_stl(path: Path, poly: vtk.vtkPolyData) -> None:
    """写二进制 STL；单位为 Å。"""
    writer = vtk.vtkSTLWriter()
    writer.SetFileName(str(path))
    writer.SetFileTypeToBinary()
    writer.SetInputData(poly)
    if writer.Write() != 1:
        raise RuntimeError(f"STL 写入失败：{path}")


def shaded_collection(vertices: np.ndarray, faces: np.ndarray, color: str):
    """生成论文预览图用的带光照三角面。"""
    triangles = vertices[faces]
    normals = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    lengths = np.linalg.norm(normals, axis=1)
    lengths[lengths == 0] = 1.0
    normals /= lengths[:, None]
    light = np.asarray((0.34, -0.47, 0.81), dtype=float)
    light /= np.linalg.norm(light)
    intensity = 0.70 + 0.30 * np.maximum(normals @ light, 0.0)
    base = np.asarray(to_rgb(color))
    face_colors = np.clip(intensity[:, None] * base[None, :], 0.0, 1.0)
    return Poly3DCollection(triangles, facecolors=face_colors, edgecolors="none", linewidths=0.0)


def draw_box(ax, extent_a: tuple[float, float, float], color="#263746",
             linestyle="-", linewidth=1.0) -> None:
    """在 3D 坐标轴上绘制长方体边框。"""
    x, y, z = extent_a
    corners = np.array([
        [0,0,0],[x,0,0],[x,y,0],[0,y,0],
        [0,0,z],[x,0,z],[x,y,z],[0,y,z],
    ], dtype=float)
    edges = ((0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),
             (0,4),(1,5),(2,6),(3,7))
    for i, j in edges:
        ax.plot(*corners[[i,j]].T, color=color, linestyle=linestyle, linewidth=linewidth)


def set_3d_axes(ax, limits_a, camera=(23, -55)) -> None:
    """固定三维图范围和相机角度，隐藏坐标轴装饰。"""
    ax.set_xlim(limits_a[0])
    ax.set_ylim(limits_a[1])
    ax.set_zlim(limits_a[2])
    ax.set_box_aspect(tuple(max(1.0, hi - lo) for lo, hi in limits_a))
    ax.view_init(elev=camera[0], azim=camera[1])
    ax.set_axis_off()
