#!/usr/bin/env python3
"""Render the whole Ti3SiC2 skeleton and a linked local branch crop in OVITO."""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from ovito import scene
from ovito.data import DataCollection
from ovito.pipeline import Pipeline, StaticSource
from ovito.vis import StandardRenderer, Viewport

HERE = Path(__file__).resolve().parent
WHOLE_MESH = HERE / "Ti3SiC2_whole_surface_mesh.npz"
LOCAL_MESH = HERE / "Ti3SiC2_local_branch_surface_mesh.npz"

WHOLE_PNG = HERE / "Ti3SiC2_whole_skeleton_OVITO.png"
WHOLE_ROI_PNG = HERE / "Ti3SiC2_whole_skeleton_with_ROI.png"
LOCAL_PNG = HERE / "Ti3SiC2_local_branch_OVITO.png"
LOCAL_GAP_PNG = HERE / "Ti3SiC2_local_branch_spacing_9p37A_orange.png"
LOCAL_GAP_TIF = HERE / "Ti3SiC2_local_branch_spacing_9p37A_orange.tif"
LOCAL_GAP_PDF = HERE / "Ti3SiC2_local_branch_spacing_9p37A_orange.pdf"
COMPOSITE_PNG = HERE / "Ti3SiC2_whole_skeleton_with_local_branch.png"
COMPOSITE_GAP_PNG = HERE / "Ti3SiC2_whole_skeleton_local_spacing_9p37A_orange.png"
COMPOSITE_GAP_TIF = HERE / "Ti3SiC2_whole_skeleton_local_spacing_9p37A_orange.tif"
COMPOSITE_GAP_PDF = HERE / "Ti3SiC2_whole_skeleton_local_spacing_9p37A_orange.pdf"
COMPOSITE_TIF = HERE / "Ti3SiC2_whole_skeleton_with_local_branch.tif"
COMPOSITE_PDF = HERE / "Ti3SiC2_whole_skeleton_with_local_branch.pdf"
WHOLE_STATE = HERE / "Ti3SiC2_whole_skeleton.ovito"
LOCAL_STATE = HERE / "Ti3SiC2_local_branch.ovito"
REPORT = HERE / "render_report.json"

ROI = {"x": (67.5, 92.5), "y": (140.0, 165.0), "z": (102.5, 127.5)}
WHOLE_SIZE = (1600, 1600)
LOCAL_SIZE = (1100, 1100)
COLOR = (0.42, 0.53, 0.62)
BLUE = (30, 101, 192)
DIM_COLOR = (36, 48, 59)
ORANGE = (244, 112, 38)
# Explicit radian value equivalent to the existing rendered camera angle.
VIEW_FOV_RAD = 2.0 * math.atan(math.tan(38.0 / 2.0))


def load_mesh(path: Path, title: str) -> tuple[Pipeline, object]:
    mesh_data = np.load(path)
    data = DataCollection()
    mesh = data.surfaces.create(title)
    mesh.create_vertices(mesh_data["vertices"])
    mesh.create_faces(mesh_data["faces"])
    mesh.connect_opposite_halfedges()
    mesh.vis.enabled = True
    mesh.vis.surface_color = COLOR
    mesh.vis.smooth_shading = True
    mesh.vis.show_cap = True
    mesh.vis.highlight_edges = False
    pipeline = Pipeline(source=StaticSource(data=data))
    return pipeline, mesh


def make_renderer() -> StandardRenderer:
    renderer = StandardRenderer()
    renderer.antialiasing_level = 4
    renderer.ambient_occlusion = True
    renderer.ambient_occlusion_intensity = 0.62
    renderer.ambient_occlusion_radius = 2.5
    renderer.shadows = True
    renderer.shadow_intensity = 0.16
    renderer.order_independent_transparency = True
    return renderer


def make_viewport(center: tuple[float, float, float], direction: tuple[float, float, float], size: tuple[int, int]) -> Viewport:
    vp = Viewport(type=Viewport.Type.Perspective)
    center_v = np.asarray(center, dtype=float)
    direction_v = np.asarray(direction, dtype=float)
    vp.camera_pos = center_v - direction_v * 100.0
    vp.camera_dir = direction_v
    vp.camera_up = (0.0, 0.0, 1.0)
    vp.fov = VIEW_FOV_RAD
    vp.zoom_all(size=size)
    return vp


def projected_point(point: tuple[float, float, float], vp: Viewport, size: tuple[int, int]) -> tuple[float, float]:
    """Project a world-space point using the final OVITO perspective camera."""
    width, height = size
    camera = np.asarray(vp.camera_pos, dtype=float)
    forward = np.array(vp.camera_dir, dtype=float, copy=True)
    forward /= np.linalg.norm(forward)
    camera_up = np.array(vp.camera_up, dtype=float, copy=True)
    camera_up /= np.linalg.norm(camera_up)
    right = np.cross(forward, camera_up)
    right /= np.linalg.norm(right)
    up = np.cross(right, forward)
    relative = np.asarray(point, dtype=float) - camera
    depth = float(np.dot(relative, forward))
    # OVITO stores perspective fov in radians.
    tangent = math.tan(float(vp.fov) / 2.0)
    aspect = width / height
    px = width / 2.0 + (width / 2.0) * float(np.dot(relative, right)) / (depth * tangent * aspect)
    py = height / 2.0 - (height / 2.0) * float(np.dot(relative, up)) / (depth * tangent)
    return px, py


def surface_point_on_camera_ray(
    mesh_vertices: np.ndarray,
    mesh_faces: np.ndarray,
    pixel: tuple[float, float],
    vp: Viewport,
    size: tuple[int, int],
) -> np.ndarray:
    """Return the first displayed-mesh surface intersection under a pixel."""
    width, height = size
    camera = np.asarray(vp.camera_pos, dtype=float)
    forward = np.asarray(vp.camera_dir, dtype=float, copy=True)
    forward /= np.linalg.norm(forward)
    camera_up = np.asarray(vp.camera_up, dtype=float, copy=True)
    camera_up /= np.linalg.norm(camera_up)
    right = np.cross(forward, camera_up)
    right /= np.linalg.norm(right)
    up = np.cross(right, forward)
    tangent = math.tan(float(vp.fov) / 2.0)
    aspect = width / height
    px, py = pixel
    x = (px - width / 2.0) / (width / 2.0) * tangent * aspect
    y = (height / 2.0 - py) / (height / 2.0) * tangent
    ray_direction = forward + x * right + y * up
    ray_direction /= np.linalg.norm(ray_direction)

    triangles = np.asarray(mesh_vertices, dtype=float)[np.asarray(mesh_faces, dtype=int)]
    edge1 = triangles[:, 1] - triangles[:, 0]
    edge2 = triangles[:, 2] - triangles[:, 0]
    hvec = np.cross(np.broadcast_to(ray_direction, edge2.shape), edge2)
    determinant = np.einsum("ij,ij->i", edge1, hvec)
    valid = np.abs(determinant) > 1.0e-12
    inverse = np.zeros_like(determinant)
    inverse[valid] = 1.0 / determinant[valid]
    offset = camera - triangles[:, 0]
    u = inverse * np.einsum("ij,ij->i", offset, hvec)
    qvec = np.cross(offset, edge1)
    v = inverse * np.einsum("j,ij->i", ray_direction, qvec)
    distance = inverse * np.einsum("ij,ij->i", edge2, qvec)
    hit = valid & (u >= -1.0e-9) & (v >= -1.0e-9) & (u + v <= 1.0 + 1.0e-9) & (distance > 0.0)
    hit_ids = np.flatnonzero(hit)
    if len(hit_ids) == 0:
        raise RuntimeError(f"No displayed surface under dimension endpoint pixel {pixel}")
    first_id = hit_ids[np.argmin(distance[hit_ids])]
    return camera + ray_direction * distance[first_id]


def dashed_line(
    draw: ImageDraw.ImageDraw,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    fill: tuple[int, int, int],
    width: int,
    dash: float = 13.0,
    gap: float = 9.0,
    halo: bool = True,
) -> None:
    x0, y0 = start
    x1, y1 = end
    dx, dy = x1 - x0, y1 - y0
    length = math.hypot(dx, dy)
    if length <= 0:
        return
    ux, uy = dx / length, dy / length
    offset = 0.0
    while offset < length:
        end_offset = min(length, offset + dash)
        segment = (x0 + ux * offset, y0 + uy * offset, x0 + ux * end_offset, y0 + uy * end_offset)
        if halo:
            draw.line(segment, fill=(255, 255, 255), width=width + 5)
        draw.line(segment, fill=fill, width=width)
        offset += dash + gap


def draw_dimension(
    image: Image.Image,
    start: tuple[float, float],
    end: tuple[float, float],
    label: str,
    font: ImageFont.FreeTypeFont,
    *,
    offset: tuple[float, float] = (0.0, 0.0),
    label_offset: tuple[float, float] = (14.0, -16.0),
    extension_lines: bool = True,
    center_label: bool = False,
) -> None:
    """Draw a dimension arrow offset from a projected world-space edge."""
    draw = ImageDraw.Draw(image)
    dx, dy = offset
    a = (start[0] + dx, start[1] + dy)
    b = (end[0] + dx, end[1] + dy)
    vx, vy = b[0] - a[0], b[1] - a[1]
    norm = math.hypot(vx, vy)
    if norm == 0:
        return
    ux, uy = vx / norm, vy / norm
    px, py = -uy, ux
    # Use either extension lines (for a nearby measured edge) or short end
    # ticks (for a clean dimension marker in the panel's white margin).
    if extension_lines:
        for p0, p1 in ((start, a), (end, b)):
            draw.line((p0, p1), fill=(255, 255, 255), width=7)
            draw.line((p0, p1), fill=DIM_COLOR, width=2)
    else:
        for tip in (a, b):
            tick = (tip[0] + px * 13.0, tip[1] + py * 13.0)
            other = (tip[0] - px * 13.0, tip[1] - py * 13.0)
            draw.line((tick, other), fill=(255, 255, 255), width=12)
            draw.line((tick, other), fill=DIM_COLOR, width=4)
    draw.line((a, b), fill=(255, 255, 255), width=14)
    draw.line((a, b), fill=DIM_COLOR, width=8)

    def head(tip: tuple[float, float], other: tuple[float, float]) -> None:
        vx, vy = other[0] - tip[0], other[1] - tip[1]
        norm = math.hypot(vx, vy)
        if norm == 0:
            return
        ux, uy = vx / norm, vy / norm
        px, py = -uy, ux
        length, half_width = 24.0, 11.0
        base = (tip[0] + ux * length, tip[1] + uy * length)
        white = [
            tip,
            (base[0] + px * (half_width + 3), base[1] + py * (half_width + 3)),
            (base[0] - px * (half_width + 3), base[1] - py * (half_width + 3)),
        ]
        dark = [
            tip,
            (base[0] + px * half_width, base[1] + py * half_width),
            (base[0] - px * half_width, base[1] - py * half_width),
        ]
        draw.polygon(white, fill=(255, 255, 255))
        draw.polygon(dark, fill=DIM_COLOR)

    head(a, b)
    head(b, a)

    box = draw.textbbox((0, 0), label, font=font)
    text_w, text_h = box[2] - box[0], box[3] - box[1]
    mid_x, mid_y = (a[0] + b[0]) / 2, (a[1] + b[1]) / 2
    if center_label:
        mid = (mid_x - text_w / 2 + label_offset[0], mid_y + label_offset[1])
    else:
        mid = (mid_x + label_offset[0], mid_y + label_offset[1])
    rect = (mid[0] - 6, mid[1] - 4, mid[0] + text_w + 6, mid[1] + text_h + 5)
    draw.rounded_rectangle(rect, radius=5, fill=(255, 255, 255), outline=(188, 198, 207), width=1)
    draw.text(mid, label, font=font, fill=DIM_COLOR)


def draw_branch_gap_annotation(
    image: Image.Image,
    start: tuple[float, float],
    end: tuple[float, float],
    label: str,
    font: ImageFont.FreeTypeFont,
) -> None:
    """Draw a prominent orange, double-headed marker across a branch gap."""
    draw = ImageDraw.Draw(image)
    dx, dy = end[0] - start[0], end[1] - start[1]
    length = math.hypot(dx, dy)
    if length == 0:
        return
    ux, uy = dx / length, dy / length
    px, py = -uy, ux
    # White keyline keeps the orange annotation visible over both the white
    # pore and the dark surface where arrowheads meet each branch.
    draw.line((start, end), fill=(255, 255, 255), width=16)
    draw.line((start, end), fill=ORANGE, width=9)

    def head(tip: tuple[float, float], direction: tuple[float, float]) -> None:
        vx, vy = direction
        perp = (-vy, vx)
        base = (tip[0] + vx * 27.0, tip[1] + vy * 27.0)
        outer = [
            tip,
            (base[0] + perp[0] * 15.0, base[1] + perp[1] * 15.0),
            (base[0] - perp[0] * 15.0, base[1] - perp[1] * 15.0),
        ]
        inner = [
            tip,
            (base[0] + perp[0] * 11.0, base[1] + perp[1] * 11.0),
            (base[0] - perp[0] * 11.0, base[1] - perp[1] * 11.0),
        ]
        draw.polygon(outer, fill=(255, 255, 255))
        draw.polygon(inner, fill=ORANGE)

    head(start, (ux, uy))
    head(end, (-ux, -uy))

    box = draw.textbbox((0, 0), label, font=font, stroke_width=2)
    text_w, text_h = box[2] - box[0], box[3] - box[1]
    text_x = (start[0] + end[0] - text_w) / 2.0
    text_y = (start[1] + end[1]) / 2.0 - text_h - 29.0
    draw.text(
        (text_x, text_y),
        label,
        font=font,
        fill=ORANGE,
        stroke_width=3,
        stroke_fill=(255, 255, 255),
    )


def draw_roi_box(image: Image.Image, vp: Viewport, size: tuple[int, int]) -> list[tuple[float, float]]:
    draw = ImageDraw.Draw(image)
    xs, ys, zs = ROI["x"], ROI["y"], ROI["z"]
    corners = []
    for ix in (0, 1):
        for iy in (0, 1):
            for iz in (0, 1):
                world = (xs[ix], ys[iy], zs[iz])
                corners.append(projected_point(world, vp, size))
    # Use a clear projected selection frame, as in the reference figure.
    left = min(p[0] for p in corners)
    right = max(p[0] for p in corners)
    top = min(p[1] for p in corners)
    bottom = max(p[1] for p in corners)
    frame = ((left, top), (right, top), (right, bottom), (left, bottom))
    for i, start in enumerate(frame):
        dashed_line(draw, start, frame[(i + 1) % len(frame)], fill=BLUE, width=5, dash=11, gap=7)
    return corners


def render_view(vp: Viewport, path: Path, size: tuple[int, int], renderer: StandardRenderer) -> Image.Image:
    rendered = vp.render_image(
        size=size,
        filename=str(path),
        background=(1.0, 1.0, 1.0),
        alpha=False,
        renderer=renderer,
        stop_on_error=True,
    )
    rendered.save(str(path), "PNG")
    with Image.open(path) as source:
        rgb = source.convert("RGB")
        rgb.save(path, "PNG", dpi=(300, 300))
    return rgb


def save_print_formats(image: Image.Image) -> None:
    image.save(COMPOSITE_PNG, format="PNG", dpi=(300, 300), optimize=True)
    image.save(COMPOSITE_TIF, format="TIFF", dpi=(300, 300), compression="tiff_lzw")
    image.save(COMPOSITE_PDF, format="PDF", resolution=300.0)


def main() -> None:
    renderer = make_renderer()
    font_path = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
    title_font = ImageFont.truetype(font_path, 58, index=0)
    panel_font = ImageFont.truetype(font_path, 42, index=0)
    subtitle_font = ImageFont.truetype(font_path, 29, index=0)
    small_font = ImageFont.truetype(font_path, 23, index=0)
    gap_font = ImageFont.truetype(font_path, 36, index=0)

    whole_pipeline, whole_mesh = load_mesh(WHOLE_MESH, "Whole Ti3SiC2 skeleton from QSGS mask")
    whole_pipeline.add_to_scene()
    whole_vp = make_viewport((100.0, 100.0, 100.0), (-1.0, -1.08, -0.74), WHOLE_SIZE)
    scene.viewports.active_vp = whole_vp
    whole_img = render_view(whole_vp, WHOLE_PNG, WHOLE_SIZE, renderer)
    whole_with_roi = whole_img.copy()
    roi_pixels = draw_roi_box(whole_with_roi, whole_vp, WHOLE_SIZE)
    whole_with_roi.save(WHOLE_ROI_PNG, "PNG", dpi=(300, 300))

    whole_pipeline.remove_from_scene()
    local_pipeline, local_mesh = load_mesh(LOCAL_MESH, "Connected local Ti3SiC2 branch crop")
    local_pipeline.add_to_scene()
    local_center = (80.0, 152.5, 115.0)
    local_vp = make_viewport(local_center, (-1.0, -1.0, -0.57), LOCAL_SIZE)
    scene.viewports.active_vp = local_vp
    local_img = render_view(local_vp, LOCAL_PNG, LOCAL_SIZE, renderer)
    mesh_vertices = np.load(LOCAL_MESH)["vertices"].astype(float)
    mesh_faces = np.load(LOCAL_MESH)["faces"].astype(int)
    # The first endpoint is the user's indicated point on the left branch.
    # The right endpoint is placed on the facing surface of the right branch;
    # the indicated screen location itself was still in the pore.
    gap_start_px = (410.0, 360.0)
    gap_end_px = (652.0, 340.0)
    gap_world = [
        surface_point_on_camera_ray(mesh_vertices, mesh_faces, pixel, local_vp, LOCAL_SIZE)
        for pixel in (gap_start_px, gap_end_px)
    ]
    gap_distance_A = float(np.linalg.norm(gap_world[0] - gap_world[1]))
    gap_distance_nm = gap_distance_A / 10.0
    draw_branch_gap_annotation(
        local_img,
        gap_start_px,
        gap_end_px,
        f"{gap_distance_A:.2f} Å",
        gap_font,
    )
    local_img.save(LOCAL_PNG, "PNG", dpi=(300, 300))
    local_img.save(LOCAL_GAP_PNG, "PNG", dpi=(300, 300), optimize=True)
    local_img.save(LOCAL_GAP_TIF, "TIFF", dpi=(300, 300), compression="tiff_lzw")
    local_img.save(LOCAL_GAP_PDF, "PDF", resolution=300.0)

    # Paper-style plate: the ROI box and two dashed leaders link the inset to
    # the same coordinates in the full 200 A skeleton.
    canvas = Image.new("RGB", (2920, 1900), "white")
    main_origin = (20, 175)
    detail_origin = (1780, 455)
    detail_size = (1000, 1000)
    canvas.paste(whole_with_roi, main_origin)
    detail = local_img.resize(detail_size, Image.Resampling.LANCZOS)
    canvas.paste(detail, detail_origin)
    draw = ImageDraw.Draw(canvas)
    draw.rectangle(
        (detail_origin[0] - 4, detail_origin[1] - 4,
         detail_origin[0] + detail_size[0] + 3, detail_origin[1] + detail_size[1] + 3),
        outline=BLUE,
        width=5,
    )

    # Connect the two corners facing the inset to matching inset edges.
    # The selected ROI is centered around the closest (+x,+y,+z) model corner.
    point_a = projected_point((ROI["x"][1], ROI["y"][1], ROI["z"][1]), whole_vp, WHOLE_SIZE)
    point_b = projected_point((ROI["x"][1], ROI["y"][0], ROI["z"][0]), whole_vp, WHOLE_SIZE)
    point_a = (point_a[0] + main_origin[0], point_a[1] + main_origin[1])
    point_b = (point_b[0] + main_origin[0], point_b[1] + main_origin[1])
    dashed_line(draw, point_a, (detail_origin[0] - 6, detail_origin[1] + 120), fill=BLUE, width=2, dash=14, gap=10, halo=False)
    dashed_line(draw, point_b, (detail_origin[0] - 6, detail_origin[1] + 880), fill=BLUE, width=2, dash=14, gap=10, halo=False)

    draw.text((70, 32), "Ti3SiC2 连通骨架与局部枝状结构", font=title_font, fill=(28, 42, 55), stroke_width=1)
    draw.text((72, 108), "蓝色虚线框标示从整体模型中截取的局部范围", font=subtitle_font, fill=(75, 91, 106))
    draw.text((95, 185), "整体骨架", font=panel_font, fill=(28, 42, 55), stroke_width=1)
    draw.text((95, 238), "200 × 200 × 200 Å  ·  20 nm", font=subtitle_font, fill=(75, 91, 106))
    draw.text((1785, 325), "局部主连通分叉", font=panel_font, fill=(28, 42, 55), stroke_width=1)
    draw.text((1788, 380), "25 × 25 × 25 Å  ·  2.5 nm", font=subtitle_font, fill=(75, 91, 106))
    save_print_formats(canvas)
    canvas.save(COMPOSITE_GAP_PNG, format="PNG", dpi=(300, 300), optimize=True)
    canvas.save(COMPOSITE_GAP_TIF, format="TIFF", dpi=(300, 300), compression="tiff_lzw")
    canvas.save(COMPOSITE_GAP_PDF, format="PDF", resolution=300.0)

    report = {
        "status": "OVITO_WHOLE_SKELETON_WITH_LOCAL_BRANCH_FIGURE_EXPORTED",
        "renderer": "OVITO 3.17.2 StandardRenderer with Mesa lavapipe",
        "model": "200 x 200 x 200 A material region from the existing Ag-Ti3SiC2 QSGS model",
        "local_roi_A": ROI,
        "whole_render_size_px": list(WHOLE_SIZE),
        "local_render_size_px": list(LOCAL_SIZE),
        "whole_camera": {
            "position": np.asarray(whole_vp.camera_pos).tolist(),
            "direction": np.asarray(whole_vp.camera_dir).tolist(),
            "up": np.asarray(whole_vp.camera_up).tolist(),
            "fov_deg": math.degrees(float(whole_vp.fov)),
        },
        "local_camera": {
            "position": np.asarray(local_vp.camera_pos).tolist(),
            "direction": np.asarray(local_vp.camera_dir).tolist(),
            "up": np.asarray(local_vp.camera_up).tolist(),
            "fov_deg": math.degrees(float(local_vp.fov)),
        },
        "roi_projected_pixels": [list(map(float, p)) for p in roi_pixels],
        "dimension_annotations": {
            "quantity": "3D surface-to-surface separation between the two user-indicated facing branches",
            "gap_A": round(gap_distance_A, 4),
            "gap_nm": round(gap_distance_nm, 5),
            "label": f"{gap_distance_A:.2f} Å",
            "endpoint_coordinates_A": [
                np.asarray(gap_world[0]).tolist(),
                np.asarray(gap_world[1]).tolist(),
            ],
            "endpoint_image_pixels_local": [list(gap_start_px), list(gap_end_px)],
            "method": "Euclidean distance between the first displayed-mesh surface intersections of camera rays through the two arrow tips; coordinates are in Å",
            "endpoint_method": "Ray intersections on the facing surfaces in the rendered view",
            "endpoint_pixels_local": [list(gap_start_px), list(gap_end_px)],
        },
        "geometry_note": "The surface meshes visualize the existing QSGS skeleton mask and are smoothed for display only; the source atomic model is unchanged. The labelled length is measured in the model's Å coordinate system between the two displayed branch surfaces.",
        "outputs": [
            WHOLE_PNG.name,
            WHOLE_ROI_PNG.name,
            LOCAL_PNG.name,
            LOCAL_GAP_PNG.name,
            LOCAL_GAP_TIF.name,
            LOCAL_GAP_PDF.name,
            COMPOSITE_PNG.name,
            COMPOSITE_GAP_PNG.name,
            COMPOSITE_GAP_TIF.name,
            COMPOSITE_GAP_PDF.name,
            COMPOSITE_TIF.name,
            COMPOSITE_PDF.name,
            WHOLE_STATE.name,
            LOCAL_STATE.name,
        ],
    }
    REPORT.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
