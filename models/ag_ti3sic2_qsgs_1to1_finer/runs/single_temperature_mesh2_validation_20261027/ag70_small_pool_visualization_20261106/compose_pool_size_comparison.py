"""Compose the existing early-small-pool and late-large-pool 3D snapshots."""
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


HERE = Path(__file__).resolve().parent
RUNS = HERE.parent
EARLY = HERE / "ag70_caseA_small_pool_4p52ps_3d.png"
LATE = RUNS / "ag70vol_1x_continue_to_2300K_20261028/ag70_phase_map_3d_18p636ps.png"
OUT = HERE / "ag70_caseA_small_and_large_pool_comparison.png"
NOTES = HERE / "POOL_SIZE_COMPARISON.md"

FONT_PATH = Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
FONT_BOLD_PATH = Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc")


def load_font(size: int, *, bold: bool = False):
    path = FONT_BOLD_PATH if bold and FONT_BOLD_PATH.exists() else FONT_PATH
    if path.exists():
        return ImageFont.truetype(str(path), size=size)
    return ImageFont.load_default()


def fit_image(path: Path, max_size: tuple[int, int]) -> Image.Image:
    image = Image.open(path).convert("RGB")
    image.thumbnail(max_size, Image.Resampling.LANCZOS)
    return image


def main() -> None:
    panel_w = 1900
    outer_pad = 70
    gutter = 50
    header_h = 190
    panel_title_h = 100
    image_h = 1560
    metrics_h = 210
    footer_h = 140
    canvas_w = outer_pad * 2 + panel_w * 2 + gutter
    canvas_h = header_h + panel_title_h + image_h + metrics_h + footer_h
    canvas = Image.new("RGB", (canvas_w, canvas_h), "white")
    draw = ImageDraw.Draw(canvas)

    title_font = load_font(57, bold=True)
    panel_font = load_font(38, bold=True)
    metric_font = load_font(27)
    footer_font = load_font(24)

    title = "Ag70 Case A：早期小熔池与后期大主熔池"
    draw.text((outer_pad, 35), title, font=title_font, fill="#1c2733")
    draw.text((outer_pad, 111), "同一 1× 模型的两个实际三维检查点（并非同一时刻）",
              font=metric_font, fill="#4b5563")

    panels = [
        {
            "x": outer_pad,
            "title": "早期小熔池｜t = 4.52 ps",
            "path": EARLY,
            "metrics": [
                "65 个 f_l >= 0.5 连通体素；液态体积积分 1.839 nm³",
                "包围尺寸 2.40 × 1.81 × 0.95 nm；全域峰温 1264.4 K",
            ],
        },
        {
            "x": outer_pad + panel_w + gutter,
            "title": "后期大主熔池｜t = 18.636 ps",
            "path": LATE,
            "metrics": [
                "76,392 个 f_l >= 0.5 连通体素；液态体积积分 2183.181 nm³",
                "包围尺寸 24.57 × 24.15 × 12.37 nm；全域峰温 2296.5 K",
            ],
        },
    ]

    top = header_h
    for panel in panels:
        x = panel["x"]
        draw.text((x, top + 16), panel["title"], font=panel_font, fill="#1c2733")
        image = fit_image(panel["path"], (panel_w - 20, image_h - 10))
        ix = x + (panel_w - image.width) // 2
        iy = top + panel_title_h + (image_h - image.height) // 2
        canvas.paste(image, (ix, iy))
        metric_y = top + panel_title_h + image_h + 16
        for line in panel["metrics"]:
            draw.text((x + 8, metric_y), line, font=metric_font, fill="#374151")
            metric_y += 42

    footer_y = header_h + panel_title_h + image_h + metrics_h + 13
    draw.text((outer_pad, footer_y),
              "液态团簇按 f_l >= 0.5 的六邻面连通体识别；体积积分为团簇内 sum(f_l) × 单体素体积。",
              font=footer_font, fill="#4b5563")
    draw.text((outer_pad, footer_y + 40),
              "两个面板对应不同时刻；模型不包含电弧、银蒸发、液体流动或飞溅。",
              font=footer_font, fill="#4b5563")

    canvas.save(OUT, format="PNG", optimize=True)
    NOTES.write_text(
        "# Ag70 Case A pool-size comparison\n\n"
        "This composite places two real snapshots from the same native 1× Case A model side by side. "
        "They are at different simulation times and must not be interpreted as simultaneous pools.\n\n"
        "- Early small pool: 4.52 ps; 65 connected native voxels at `f_l >= 0.5`; fraction-integrated liquid volume 1.838613 nm³; axis-aligned extent 2.397080 × 1.811145 × 0.951669 nm; peak temperature 1264.383 K.\n"
        "- Late large main pool: 18.636 ps; 76,392 connected native voxels at `f_l >= 0.5`; fraction-integrated liquid volume 2183.181305 nm³; axis-aligned extent 24.570070 × 24.148599 × 12.371697 nm; peak temperature 2296.532 K.\n\n"
        "The liquid-volume integral is `sum(f_l) × voxel_volume` over the thresholded connected component. "
        "The thermal model does not resolve an arc, evaporation, liquid flow, or splashing.\n",
        encoding="utf-8",
    )
    print(f"wrote {OUT} ({canvas.width}x{canvas.height})")


if __name__ == "__main__":
    main()
