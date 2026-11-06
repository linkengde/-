"""Render disconnected main and satellite liquid clusters from one checkpoint."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import ConnectionPatch, Patch
import numpy as np
from mpl_toolkits.mplot3d import proj3d

HERE = Path(__file__).resolve().parent
RUN_DIR = HERE.parent
MODEL_DIR = RUN_DIR.parents[1]
PHASE_DIR = MODEL_DIR / "phase_field_qsgs"
sys.path.insert(0, str(PHASE_DIR))
from qsgs_liquid_connectivity import liquid_cluster_labels

GEOMETRY = RUN_DIR / "qsgs_morphology_ag70vol.npz"
CHECKPOINT = RUN_DIR / "ag70vol_1x_continue_to_2300K_20261028/checkpoints/caseA_ag70_1x_heat_15p29ps.npz"
OUTPUT = HERE / "ag70_caseA_main_and_satellite_pools_15p29ps_3d.png"
REPORT = HERE / "SIMULTANEOUS_POOL_COMPONENTS_15P29PS.md"

COLORS = {
    "scaffold": "#263f51",
    "solid": "#83b7d5",
    "partial": "#ffd166",
    "main": "#f28e2b",
    "satellite": "#cc3d7a",
}


def draw_box(ax, lo, hi, color, lw=1.2):
    x0, y0, z0 = lo
    x1, y1, z1 = hi
    c = np.array([[x0,y0,z0],[x1,y0,z0],[x1,y1,z0],[x0,y1,z0],
                  [x0,y0,z1],[x1,y0,z1],[x1,y1,z1],[x0,y1,z1]])
    edges = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),
             (0,4),(1,5),(2,6),(3,7)]
    for i, j in edges:
        ax.plot(*c[[i,j]].T, color=color, linewidth=lw, linestyle=(0,(4,3)))


def phase_masks(scaffold, lf):
    ag = ~scaffold
    liquid = ag & (lf >= 0.5)
    partial = ag & (lf > 1e-12) & (lf < 0.5)
    solid = ag & ~(liquid | partial)
    return liquid, partial, solid


def main():
    with np.load(GEOMETRY, allow_pickle=False) as d:
        scaffold = np.asarray(d["skeleton"], dtype=bool)
        size_nm = np.asarray(d["size_A"], dtype=float) / 10.0
    with np.load(CHECKPOINT, allow_pickle=False) as d:
        temperature = np.asarray(d["temperature_K"], dtype=float)
        lf = np.asarray(d["ag_liquid_fraction"], dtype=float)
        metadata = json.loads(str(d["metadata_json"].item()))

    spacing = size_nm / np.asarray(scaffold.shape)
    voxel_volume = float(np.prod(spacing))
    liquid, partial, solid = phase_masks(scaffold, lf)
    labels, ncomp = liquid_cluster_labels(liquid)
    counts = np.bincount(labels.ravel())
    counts[0] = 0
    ids = np.argsort(counts)[::-1]
    ids = [int(i) for i in ids if counts[i] > 0]
    if len(ids) < 2:
        raise ValueError(f"Expected a main pool and separate satellite(s), found {len(ids)} component(s).")
    main_id = ids[0]
    satellite_ids = ids[1:]
    main_mask = labels == main_id
    satellite_mask = np.isin(labels, satellite_ids)
    satellite_largest_id = satellite_ids[0]
    satellite_largest = labels == satellite_largest_id

    def metrics(mask):
        idx = np.where(mask)
        lo = np.array([v.min() for v in idx])
        hi = np.array([v.max() for v in idx])
        dims = (hi - lo + 1) * spacing
        centroid = np.array([v.mean() for v in idx]) * spacing
        cells = int(mask.sum())
        liquid_volume = float(np.sum(lf[mask]) * voxel_volume)
        return {"lo": lo, "hi": hi, "dims": dims, "centroid": centroid,
                "cells": cells, "mask_volume": cells * voxel_volume,
                "liquid_volume": liquid_volume}

    main_info = metrics(main_mask)
    sat_info = metrics(satellite_largest)
    t_ps = float(metadata["time_ps"])

    # A native-resolution window around the largest satellite pool. The full
    # volume remains in the left panel, where the satellite stays at true scale.
    pad_xy = 8
    x0 = max(0, int(sat_info["lo"][0]) - pad_xy)
    x1 = min(scaffold.shape[0], int(sat_info["hi"][0]) + 1 + pad_xy)
    y0 = max(0, int(sat_info["lo"][1]) - pad_xy)
    y1 = min(scaffold.shape[1], int(sat_info["hi"][1]) + 1 + pad_xy)
    z0 = max(0, int(sat_info["lo"][2]) - 8)
    z1 = scaffold.shape[2]
    crop = (x0,x1,y0,y1,z0,z1)
    crop_sc = scaffold[x0:x1,y0:y1,z0:z1]
    crop_lf = lf[x0:x1,y0:y1,z0:z1]
    crop_liq, crop_partial, crop_solid = phase_masks(crop_sc, crop_lf)
    crop_main = main_mask[x0:x1,y0:y1,z0:z1]
    crop_sat = satellite_mask[x0:x1,y0:y1,z0:z1]
    crop_size = np.array([x1-x0,y1-y0,z1-z0]) * spacing

    font_path = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
    if Path(font_path).exists():
        font_manager.fontManager.addfont(font_path)
        font_name = font_manager.FontProperties(fname=font_path).get_name()
    else:
        font_name = "DejaVu Sans"
    plt.rcParams.update({"font.family": [font_name, "DejaVu Sans"],
                         "axes.unicode_minus": False, "figure.facecolor": "white"})

    fig = plt.figure(figsize=(16, 9.5), constrained_layout=False)
    gs = fig.add_gridspec(1, 2, left=0.025, right=0.985, top=0.80, bottom=0.16,
                          width_ratios=(1.08, 0.92), wspace=0.025)
    ax = fig.add_subplot(gs[0,0], projection="3d")

    # Coarse full-domain material context: one representative cell per 4^3 block.
    stride = 4
    edges = [np.r_[np.arange(0,n,stride), n] for n in scaffold.shape]
    centers = [(e[:-1] + e[1:] - 1)//2 for e in edges]
    phase = np.zeros(scaffold.shape, dtype=np.uint8)
    phase[solid] = 1
    phase[partial] = 2
    phase[liquid] = 3
    sampled = phase[np.ix_(*centers)]
    display = np.ones(sampled.shape, dtype=bool)
    for i in range(len(centers[0])):
        for j in range(len(centers[1])):
            for k in range(len(centers[2])):
                block_lo = (edges[0][i], edges[1][j], edges[2][k])
                block_hi = (edges[0][i+1], edges[1][j+1], edges[2][k+1])
                if all(block_lo[a] < crop[2*a+1] and block_hi[a] > crop[2*a] for a in range(3)):
                    display[i,j,k] = False
    xe, ye, ze = np.meshgrid(edges[0]*spacing[0], edges[1]*spacing[1],
                             edges[2]*spacing[2], indexing="ij")
    for val, color in ((0,COLORS["scaffold"]),(1,COLORS["solid"]),(2,COLORS["partial"])):
        ax.voxels(xe,ye,ze,(sampled==val)&display,facecolors=color,
                  edgecolors=(0.12,0.19,0.23,0.14),linewidth=0.035,shade=True)

    # Exact liquid components are overlaid at native resolution.
    nx, ny, nz = scaffold.shape
    gx, gy, gz = np.meshgrid(np.arange(nx+1)*spacing[0],
                             np.arange(ny+1)*spacing[1],
                             np.arange(nz+1)*spacing[2], indexing="ij")
    ax.voxels(gx,gy,gz,main_mask,
              facecolors=matplotlib.colors.to_rgba(COLORS["main"],0.94),
              edgecolors=(0.34,0.15,0.04,0.28),linewidth=0.045,shade=True)
    ax.voxels(gx,gy,gz,satellite_mask,
              facecolors=matplotlib.colors.to_rgba(COLORS["satellite"],1.0),
              edgecolors=(0.24,0.04,0.12,0.72),linewidth=0.09,shade=True)
    draw_box(ax, np.array([x0,y0,z0])*spacing,
             np.array([x1,y1,z1])*spacing, COLORS["satellite"], lw=1.3)
    ax.set_xlim(0,size_nm[0]); ax.set_ylim(0,size_nm[1]); ax.set_zlim(0,size_nm[2])
    ax.set_box_aspect(size_nm); ax.set_proj_type("ortho"); ax.view_init(elev=28,azim=-56)
    ax.set_xlabel("x (nm)"); ax.set_ylabel("y (nm)"); ax.set_zlabel("z (nm)")
    ax.grid(True,alpha=0.20)
    ax.set_title("同一时刻的全域位置",loc="left",fontweight="bold",pad=6)

    # Native-resolution detail around the satellite; main-pool cells within
    # the crop remain orange to show that the pink cluster is disconnected.
    axz = fig.add_subplot(gs[0,1], projection="3d")
    nxc, nyc, nzc = crop_sc.shape
    xx, yy, zz = np.meshgrid(np.arange(nxc+1)*spacing[0],
                             np.arange(nyc+1)*spacing[1],
                             np.arange(nzc+1)*spacing[2], indexing="ij")
    visible = np.ones(crop_sc.shape, dtype=bool)
    visible[nxc//2:, :nyc//2, :max(1,nzc//2)] = False
    for mask, color in ((crop_sc,COLORS["scaffold"]),
                        (crop_solid,COLORS["solid"]),
                        (crop_partial,COLORS["partial"])):
        axz.voxels(xx,yy,zz,mask&visible,facecolors=color,
                   edgecolors=(0.12,0.19,0.23,0.13),linewidth=0.04,shade=True)
    axz.voxels(xx,yy,zz,crop_main,
               facecolors=matplotlib.colors.to_rgba(COLORS["main"],0.94),
               edgecolors=(0.34,0.15,0.04,0.34),linewidth=0.06,shade=True)
    axz.voxels(xx,yy,zz,crop_sat,
               facecolors=COLORS["satellite"],
               edgecolors=(0.24,0.04,0.12,0.85),linewidth=0.13,shade=True)
    axz.set_xlim(0,crop_size[0]); axz.set_ylim(0,crop_size[1]); axz.set_zlim(0,crop_size[2])
    axz.set_box_aspect(crop_size); axz.set_proj_type("ortho"); axz.view_init(elev=30,azim=-56)
    axz.set_xlabel("局部 x (nm)"); axz.set_ylabel("局部 y (nm)"); axz.set_zlabel("局部 z (nm)")
    axz.grid(True,alpha=0.20)
    axz.set_title("主熔池边缘与独立小熔池｜原生体素放大",loc="left",fontweight="bold",pad=6)

    fig.canvas.draw()
    center = sat_info["centroid"]
    xp,yp,_ = proj3d.proj_transform(*center,ax.get_proj())
    point = ax.transAxes.inverted().transform(ax.transData.transform((xp,yp)))
    fig.add_artist(ConnectionPatch(xyA=point,coordsA=ax.transAxes,xyB=(0.42,1.0),
                                   coordsB=axz.transAxes,color=COLORS["satellite"],
                                   linestyle=(0,(4,3)),linewidth=1.35,clip_on=False))

    fig.suptitle(f"Ag70 Case A：同一时刻的大主熔池与独立小熔池｜t = {t_ps:.2f} ps",
                 fontsize=17,fontweight="bold",y=0.965)
    legend = [Patch(facecolor=COLORS["scaffold"],label="Ti₃SiC₂ 骨架"),
              Patch(facecolor=COLORS["solid"],label="固态 Ag"),
              Patch(facecolor=COLORS["partial"],label="部分液化 Ag"),
              Patch(facecolor=COLORS["main"],label=f"主熔池 ({main_info['cells']:,} 体素)"),
              Patch(facecolor=COLORS["satellite"],label=f"独立小熔池 ({sat_info['cells']:,} 体素)")]
    fig.legend(handles=legend,loc="upper center",bbox_to_anchor=(0.5,0.91),
               ncol=5,frameon=False,fontsize=9)

    metrics_text = (
        f"峰温 {temperature.max():.1f} K；液态连通分量 {ncomp} 个。\n"
        f"主熔池：{main_info['cells']:,} 体素，液态体积积分 {main_info['liquid_volume']:.3f} nm³，"
        f"包围尺寸 {main_info['dims'][0]:.2f}×{main_info['dims'][1]:.2f}×{main_info['dims'][2]:.2f} nm。\n"
        f"最大独立小熔池：{sat_info['cells']} 体素，液态体积积分 {sat_info['liquid_volume']:.3f} nm³，"
        f"质心 ({sat_info['centroid'][0]:.2f}, {sat_info['centroid'][1]:.2f}, {sat_info['centroid'][2]:.2f}) nm。"
    )
    fig.text(0.035,0.075,metrics_text,fontsize=9.5,ha="left",va="center",
             bbox=dict(boxstyle="round,pad=0.45",facecolor="white",alpha=0.94,edgecolor="#d8dee5"))
    fig.text(0.035,0.018,
             "左图背景按 4×4×4 体素代表性取样，液相团簇叠加原生网格；右图为小熔池邻域原生体素放大。"
             "模型不含电弧、蒸发、液体流动或飞溅。",
             fontsize=8.8,color="#4b5563",ha="left",va="center")
    fig.savefig(OUTPUT,dpi=210,bbox_inches="tight")
    plt.close(fig)

    report = f"""# Same-time main and satellite Ag pools at {t_ps:.6f} ps

- Model: Ag70 vol%, native 1× QSGS geometry, single-temperature enthalpy method, Case A.
- Checkpoint: `{CHECKPOINT.name}`; this image uses one checkpoint only.
- Peak temperature: {temperature.max():.6f} K.
- Connected liquid components at `f_l >= 0.5` with face-neighbor connectivity: {ncomp}.
- Main pool: {main_info['cells']} voxels; threshold-mask volume {main_info['mask_volume']:.6f} nm³; fraction-integrated liquid volume {main_info['liquid_volume']:.6f} nm³; dimensions {main_info['dims'][0]:.6f} × {main_info['dims'][1]:.6f} × {main_info['dims'][2]:.6f} nm; centroid {main_info['centroid'].tolist()} nm.
- Largest separate satellite: {sat_info['cells']} voxels; threshold-mask volume {sat_info['mask_volume']:.6f} nm³; fraction-integrated liquid volume {sat_info['liquid_volume']:.6f} nm³; dimensions {sat_info['dims'][0]:.6f} × {sat_info['dims'][1]:.6f} × {sat_info['dims'][2]:.6f} nm; centroid {sat_info['centroid'].tolist()} nm.
- Other separate components: {len(satellite_ids)-1}.

The left panel locates both clusters in the same full-domain state. Its matrix background is sampled once per 4×4×4 native-voxel block, while the liquid clusters are shown at native resolution. The right panel shows the satellite neighborhood at native resolution. These rendering choices do not alter the saved thermal fields.

This thermal model does not include arc-root motion, evaporation, liquid flow, or splashing. The spatially fixed Gaussian source is not a moving arc model.
"""
    REPORT.write_text(report,encoding="utf-8")
    print(f"saved={OUTPUT}")
    print(f"t={t_ps:.6f} ps components={ncomp} main_voxels={main_info['cells']} satellite_voxels={sat_info['cells']}")


if __name__ == "__main__":
    main()
