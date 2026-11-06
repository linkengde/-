"""Render the first clearly resolved small Ag pool from the existing Case A run."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import ConnectionPatch, Patch
import numpy as np
from mpl_toolkits.mplot3d import proj3d

import sys
PHASE_DIR = Path(__file__).resolve().parents[3] / "phase_field_qsgs"
sys.path.insert(0, str(PHASE_DIR))
from qsgs_liquid_connectivity import liquid_cluster_labels


OUT_DIR = Path(__file__).resolve().parent
RUN_DIR = OUT_DIR.parent
GEOMETRY = RUN_DIR / "qsgs_morphology_ag70vol.npz"
CHECKPOINT = RUN_DIR / "ag70vol_resume_1x_from1ps_20261028/checkpoints/caseA_ag70_1x_final_4p52ps.npz"
PNG = OUT_DIR / "ag70_caseA_small_pool_4p52ps_3d.png"
NOTES = OUT_DIR / "SMALL_POOL_3D_NOTES.md"

COLORS = {
    "scaffold": "#263f51",
    "solid_ag": "#83b7d5",
    "partial_ag": "#ffd166",
    "liquid_ag": "#f28e2b",
}


def draw_box(ax, lo, hi, color="#f28e2b", lw=1.2, ls="--"):
    x0, y0, z0 = lo
    x1, y1, z1 = hi
    c = np.array([[x0,y0,z0],[x1,y0,z0],[x1,y1,z0],[x0,y1,z0],
                  [x0,y0,z1],[x1,y0,z1],[x1,y1,z1],[x0,y1,z1]])
    edges = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]
    for i,j in edges:
        ax.plot(*c[[i,j]].T, color=color, linewidth=lw, linestyle=ls)


def make_edges(n, spacing, stride=1):
    e = np.r_[np.arange(0, n, stride), n]
    return e * spacing


def set_view(ax, size_nm, *, labels=True):
    ax.set_xlim(0, size_nm[0]); ax.set_ylim(0, size_nm[1]); ax.set_zlim(0, size_nm[2])
    ax.set_box_aspect(size_nm)
    ax.set_proj_type("ortho")
    ax.view_init(elev=28, azim=-56)
    if labels:
        ax.set_xlabel("x (nm)", labelpad=3)
        ax.set_ylabel("y (nm)", labelpad=3)
        ax.set_zlabel("z (nm)", labelpad=3)
    ax.grid(True, alpha=0.22)


def phase_masks(scaffold, lf):
    ag = ~scaffold
    liquid = ag & (lf >= 0.5)
    partial = ag & (lf > 1e-12) & (lf < 0.5)
    solid = ag & ~(liquid | partial)
    return ag, liquid, partial, solid


def add_coarse_context(ax, scaffold, lf, spacing_nm, crop_bounds, stride=4):
    nx, ny, nz = scaffold.shape
    edges = [np.r_[np.arange(0, n, stride), n] for n in scaffold.shape]
    centers = [(e[:-1] + e[1:] - 1) // 2 for e in edges]
    ag, liquid, partial, solid = phase_masks(scaffold, lf)
    phase = np.zeros(scaffold.shape, dtype=np.uint8)
    phase[solid] = 1
    phase[partial] = 2
    phase[liquid] = 3
    sampled = phase[np.ix_(*centers)]
    masks = (sampled == 0, sampled == 1, sampled == 2)

    # A lower-front corner is removed to expose interior phases. The small
    # pool at the upper surface remains intact. Its larger ROI is opened in
    # the coarse context and shown at native resolution in the right inset.
    display = np.ones(sampled.shape, dtype=bool)
    display[len(centers[0]) // 2:, :len(centers[1]) // 2, :len(centers[2]) // 2] = False
    x0, x1, y0, y1, z0, z1 = crop_bounds
    for i, x in enumerate(centers[0]):
        for j, y in enumerate(centers[1]):
            for k, z in enumerate(centers[2]):
                block_lo = (edges[0][i], edges[1][j], edges[2][k])
                block_hi = (edges[0][i+1], edges[1][j+1], edges[2][k+1])
                intersects_roi = all(block_lo[a] < crop_bounds[2*a+1] and block_hi[a] > crop_bounds[2*a] for a in range(3))
                if intersects_roi:
                    display[i,j,k] = False

    x, y, z = np.meshgrid(edges[0]*spacing_nm[0], edges[1]*spacing_nm[1],
                          edges[2]*spacing_nm[2], indexing="ij")
    specs = ((masks[0], COLORS["scaffold"]), (masks[1], COLORS["solid_ag"]),
             (masks[2], COLORS["partial_ag"]))
    for mask, color in specs:
        ax.voxels(x, y, z, mask & display,
                  facecolors=color, edgecolors=(0.12,0.19,0.23,0.13),
                  linewidth=0.035, shade=True)

    # Overlay the exact 65-cell liquid cluster at native resolution.
    native_pool = liquid[x0:x1, y0:y1, z0:z1]
    xe = np.arange(x0, x1 + 1) * spacing_nm[0]
    ye = np.arange(y0, y1 + 1) * spacing_nm[1]
    ze = np.arange(z0, z1 + 1) * spacing_nm[2]
    xe3, ye3, ze3 = np.meshgrid(xe, ye, ze, indexing="ij")
    ax.voxels(xe3, ye3, ze3, native_pool,
              facecolors=matplotlib.colors.to_rgba(COLORS["liquid_ag"], 0.96),
              edgecolors=(0.35,0.16,0.06,0.5), linewidth=0.07, shade=True)

    box_lo = np.array([x0,y0,z0]) * spacing_nm
    box_hi = np.array([x1,y1,z1]) * spacing_nm
    draw_box(ax, box_lo, box_hi, color=COLORS["liquid_ag"], lw=1.15)
    return display.shape


def main():
    with np.load(GEOMETRY, allow_pickle=False) as data:
        scaffold = np.asarray(data["skeleton"], dtype=bool)
        size_nm = np.asarray(data["size_A"], dtype=float) / 10.0
    with np.load(CHECKPOINT, allow_pickle=False) as data:
        temperature = np.asarray(data["temperature_K"], dtype=float)
        lf = np.asarray(data["ag_liquid_fraction"], dtype=float)
        metadata = json.loads(str(data["metadata_json"].item()))

    spacing_nm = size_nm / np.asarray(scaffold.shape, dtype=float)
    ag, liquid, partial, solid = phase_masks(scaffold, lf)
    labels, ncomp = liquid_cluster_labels(liquid)
    counts = np.bincount(labels.ravel())
    counts[0] = 0
    pool_id = int(np.argmax(counts))
    pool = labels == pool_id
    pool_voxels = int(counts[pool_id])
    voxel_volume = float(np.prod(spacing_nm))
    pool_volume = pool_voxels * voxel_volume
    pool_liquid_volume = float(np.sum(lf[pool]) * voxel_volume)
    xyz = np.where(pool)
    lo = np.array([int(c.min()) for c in xyz])
    hi = np.array([int(c.max()) for c in xyz])
    dims_nm = (hi - lo + 1) * spacing_nm

    # Native-resolution local window around the top-surface pool.
    pad_xy, pad_z = 8, 8
    crop = (max(0,lo[0]-pad_xy), min(scaffold.shape[0],hi[0]+1+pad_xy),
            max(0,lo[1]-pad_xy), min(scaffold.shape[1],hi[1]+1+pad_xy),
            max(0,lo[2]-pad_z), scaffold.shape[2])
    x0,x1,y0,y1,z0,z1 = crop
    crop_sc = scaffold[x0:x1,y0:y1,z0:z1]
    crop_lf = lf[x0:x1,y0:y1,z0:z1]
    crop_ag, crop_liq, crop_part, crop_solid = phase_masks(crop_sc,crop_lf)
    crop_size = np.array([x1-x0,y1-y0,z1-z0]) * spacing_nm
    ex=np.arange(x1-x0+1)*spacing_nm[0]
    ey=np.arange(y1-y0+1)*spacing_nm[1]
    ez=np.arange(z1-z0+1)*spacing_nm[2]
    X,Y,Z=np.meshgrid(ex,ey,ez,indexing="ij")

    fontfile="/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
    if Path(fontfile).exists():
        font_manager.fontManager.addfont(fontfile)
        cjk=font_manager.FontProperties(fname=fontfile).get_name()
    else:
        cjk="DejaVu Sans"
    plt.rcParams.update({"font.family":[cjk,"DejaVu Sans"],"axes.unicode_minus":False,
                         "figure.facecolor":"white","font.size":10})

    fig=plt.figure(figsize=(15,9.2), constrained_layout=False)
    gs=fig.add_gridspec(1,2,left=0.035,right=0.98,top=0.82,bottom=0.15,
                        width_ratios=(1.05,0.95),wspace=0.06)
    ax=fig.add_subplot(gs[0,0],projection="3d")
    add_coarse_context(ax, scaffold, lf, spacing_nm, crop, stride=4)
    set_view(ax,size_nm)
    ax.set_title("Ag70 Case A 全域位置（背景每 4×4×4 体素取样）",loc="left",
                 color="#1c2733",fontweight="bold",pad=7)
    ax.text2D(0.02,0.02,
              f"t = {metadata['time_ps']:.2f} ps；全域峰温 {temperature.max():.1f} K\n"
              f"橙色液态团簇在全域中很小，右侧为原生体素放大",
              transform=ax.transAxes,fontsize=9,
              bbox=dict(boxstyle="round,pad=0.35",facecolor="white",alpha=0.92,edgecolor="#d8dee5"))

    axz=fig.add_subplot(gs[0,1],projection="3d")
    crop_phase=np.zeros(crop_sc.shape,dtype=np.uint8)
    crop_phase[~crop_sc]=1
    crop_phase[crop_solid]=1; crop_phase[crop_part]=2; crop_phase[crop_liq]=3
    visible=np.ones(crop_phase.shape,dtype=bool)
    # Clip away the near/front lower corner to expose the small pool and scaffold.
    visible[crop_phase.shape[0]//2:, :crop_phase.shape[1]//2,
            :max(1, crop_phase.shape[2]//2)] = False
    specs=((crop_phase==0,COLORS["scaffold"]),(crop_phase==1,COLORS["solid_ag"]),
           (crop_phase==2,COLORS["partial_ag"]),(crop_phase==3,COLORS["liquid_ag"]))
    for idx,(mask,color) in enumerate(specs):
        rgba=matplotlib.colors.to_rgba(color,0.58 if idx==3 else 1.0)
        axz.voxels(X,Y,Z,mask & visible,facecolors=rgba,
                   edgecolors=(0.12,0.19,0.23,0.15),linewidth=0.04,shade=True)
    # Put exact pool outline on top so the small cluster remains easy to find.
    pool_crop=pool[x0:x1,y0:y1,z0:z1]
    axz.voxels(X, Y, Z, pool_crop,
               facecolors=matplotlib.colors.to_rgba(COLORS["liquid_ag"],0.94),
               edgecolors=(0.35,0.16,0.06,0.65),linewidth=0.11,shade=True)
    # Local axes use the same physical lengths as the exact native grid.
    axz.set_xlim(0,crop_size[0]); axz.set_ylim(0,crop_size[1]); axz.set_zlim(0,crop_size[2])
    axz.set_box_aspect(crop_size); axz.set_proj_type("ortho"); axz.view_init(elev=30,azim=-56)
    axz.set_xlabel("局部 x (nm)",labelpad=2); axz.set_ylabel("局部 y (nm)",labelpad=2); axz.set_zlabel("局部 z (nm)",labelpad=2)
    axz.grid(True,alpha=0.2)
    axz.set_title("(b) 小熔池局部放大：原生体素",loc="left",color="#1c2733",fontweight="bold",pad=7)
    axz.text2D(0.02,0.02,
               f"{pool_voxels} 个 fₗ ≥ 0.5 体素；液态体积积分 {pool_liquid_volume:.3f} nm³\n"
               f"体素掩膜体积 {pool_volume:.3f} nm³；"
               f"包围尺寸约 {dims_nm[0]:.2f} × {dims_nm[1]:.2f} × {dims_nm[2]:.2f} nm",
               transform=axz.transAxes,fontsize=9,
               bbox=dict(boxstyle="round,pad=0.35",facecolor="white",alpha=0.92,edgecolor="#d8dee5"))

    # Connect the highlighted ROI in the full cell to the local crop.
    fig.canvas.draw()
    pool_center_nm=((lo+hi+1)/2)*spacing_nm
    xp,yp,_=proj3d.proj_transform(*pool_center_nm,ax.get_proj())
    p=ax.transAxes.inverted().transform(ax.transData.transform((xp,yp)))
    fig.add_artist(ConnectionPatch(xyA=p,coordsA=ax.transAxes,xyB=(0.5,1.0),coordsB=axz.transAxes,
                                   color=COLORS["liquid_ag"],linestyle=(0,(4,3)),linewidth=1.5,clip_on=False))

    fig.suptitle("Ag70 局部小熔池形成初期｜单温焓法 Case A",fontsize=17,fontweight="bold",color="#1c2733",y=0.96)
    fig.legend(handles=[Patch(facecolor=COLORS["scaffold"],label="Ti₃SiC₂ 骨架"),
                        Patch(facecolor=COLORS["solid_ag"],label="固态 Ag"),
                        Patch(facecolor=COLORS["partial_ag"],label="部分液化 Ag (0 < fₗ < 0.5)"),
                        Patch(facecolor=COLORS["liquid_ag"],label="液态 Ag (fₗ ≥ 0.5)")],
               loc="upper center",bbox_to_anchor=(0.5,0.885),ncol=4,frameon=False,fontsize=9)
    fig.text(0.035,0.065,
             "这是 4.52 ps 实际检查点的热相变快照；不是电弧、液体流动、银蒸发或飞溅模拟。全域背景仅用于定位，"
             "池体数量、体积和尺寸均由未降采样的原生网格计算。",
             fontsize=9,color="#4b5563",ha="left",va="center")
    fig.savefig(PNG,dpi=220,bbox_inches="tight")
    plt.close(fig)

    note=f"""# Ag70 small-pool visualization at 4.52 ps

- Source: existing Ag70, 1×, Case A checkpoint `{CHECKPOINT.name}`; no new heat simulation was run.
- Time: `{metadata['time_ps']:.6f} ps`; peak temperature `{temperature.max():.3f} K`.
- Main connected cluster uses the solver's `f_liquid >= 0.5` threshold: `{pool_voxels}` native voxels.
- Threshold-mask geometric volume: `{pool_volume:.6f} nm³`; liquid-volume integral within that cluster (`sum(f_liquid) × voxel volume`): `{pool_liquid_volume:.6f} nm³`.
- Axis-aligned voxel extent: `{dims_nm[0]:.6f} × {dims_nm[1]:.6f} × {dims_nm[2]:.6f} nm`.
- There is one `f_l >= 0.5` connected cluster at this snapshot. The 4×4×4 global context is only for location; the right inset and reported metrics use native voxels.
- This enthalpy model contains no arc-root plasma, material evaporation, liquid flow, or splashing.
"""
    NOTES.write_text(note,encoding="utf-8")
    print(f"wrote {PNG}")
    print(f"t={metadata['time_ps']:.3f} ps, components={ncomp}, voxels={pool_voxels}, mask_volume={pool_volume:.6f} nm^3, liquid_volume={pool_liquid_volume:.6f} nm^3, dims_nm={dims_nm.tolist()}")


if __name__=="__main__":
    main()
