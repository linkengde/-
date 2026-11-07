#!/usr/bin/env python3
"""Render A/B/C near-peak pool fields from each completed hold-start snapshot."""
import csv
import importlib.util
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy import ndimage

HERE = Path(__file__).resolve().parent
PARSER = HERE / "analyze_abc_full.py"
spec = importlib.util.spec_from_file_location("full_analyzer", PARSER)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

CASES = (("A", "A · pressure off", 0),
         ("B", "B · 501 MPa", 1),
         ("C", "C · 501 MPa + evaporation", 1))
COLORS = {"skeleton": "#183645", "solid Ag": "#83b8d2",
          "partial Ag": "#ffe000", "liquid Ag": "#f18424"}


def last_row(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))[-1]


def find_row_at_time(path, t):
    with path.open(newline="") as f:
        data = list(csv.DictReader(f))
    return min(data, key=lambda r: abs(float(r["time_ps"]) - t))


def draw_case(ax, case, title, crop=False):
    folder = HERE / "abc_full" / case
    paths = sorted(folder.glob("p1h_marangoni_on_hold_start_*_t*.bin"))
    if len(paths) != 1:
        raise RuntimeError(f"{case}: expected one hold-start snapshot, found {len(paths)}")
    snap = module.read_snapshot(paths[0])
    a = snap["fields"]
    x = (snap["x0"] + (np.arange(snap["nx"]) + .5)*snap["dx_m"])*1e9
    y = (snap["y0"] + (np.arange(snap["ny"]) + .5)*snap["dx_m"])*1e9
    z = (snap["z0"] + (np.arange(snap["nz"]) + .5)*snap["dx_m"])*1e9
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    code, ag, sc, lf = a[19], a[2], a[3], a[4]
    skeleton = (code > 1.5) | (sc > .5)
    agdom = (~skeleton) & ((code > .5) | (ag > .2)) & (ag > .2)
    solid = agdom & (lf <= 1.e-3)
    partial = agdom & (lf > 1.e-3) & (lf < .5)
    liquid = agdom & (lf >= .5)
    # For the full-model view, draw the external material skin and the
    # skeleton/solid-Ag immediately bordering the melt. Scattering every
    # interior voxel overfills the projection and hides the phase colors.
    if not crop:
        material = code > .5
        inner = ndimage.binary_erosion(material,
                    structure=ndimage.generate_binary_structure(3, 1),
                    border_value=0)
        shell = material & ~inner
        near_melt = ndimage.binary_dilation(liquid | partial,
                    structure=ndimage.generate_binary_structure(3, 1),
                    iterations=1)
        skeleton = skeleton & (shell | near_melt)
        solid = solid & (shell | near_melt)
    roi = np.ones_like(skeleton, dtype=bool)
    top = max(z)
    if crop:
        i0 = int(np.argmin(np.abs(x - 9.58832)))
        j0 = int(np.argmin(np.abs(y - 9.65943950771472)))
        iz = np.where(a[0, i0, j0, :] > .5)[0]
        top = z[int(iz.max())] + snap["dx_m"]*1e9/2 if iz.size else max(z)
        roi = ((np.abs(X-9.58832) <= 6.) & (np.abs(Y-9.65943950771472) <= 6.) &
               (Z >= top-6.) & (Z <= top+snap["dx_m"]*1e9))
    # Keep the cutaway readable: render only the near-top solid-Ag skin in
    # the close view, then show the full liquid pool and dark skeleton through it.
    solid_display = solid & (Z >= top - .75) if crop else solid
    for mask, key, size, alpha, stride_xy, stride_z in (
        (skeleton, "skeleton", 5.2, .95, 1, 1),
        (solid_display, "solid Ag", 4.4, .48, 2, 2),
        (partial, "partial Ag", 6.3, .98, 1, 1),
        (liquid, "liquid Ag", 6.5, .78, 1, 1),
    ):
        mask = mask & roi
        if crop:
            size = {"skeleton": 8., "solid Ag": 8., "partial Ag": 17., "liquid Ag": 18.}[key]
            alpha = {"skeleton": .62, "solid Ag": .16, "partial Ag": .98, "liquid Ag": .90}[key]
            stride_xy, stride_z = (2, 1) if key == "skeleton" else (3, 2) if key == "solid Ag" else (1, 1)
        else:
            size = {"skeleton": 4.0, "solid Ag": 3.8, "partial Ag": 6.0, "liquid Ag": 6.3}[key]
            alpha = {"skeleton": .84, "solid Ag": .38, "partial Ag": .98, "liquid Ag": .86}[key]
        if stride_xy > 1 or stride_z > 1:
            m = mask[::stride_xy, ::stride_xy, ::stride_z]
            xp, yp, zp = X[::stride_xy, ::stride_xy, ::stride_z], Y[::stride_xy, ::stride_xy, ::stride_z], Z[::stride_xy, ::stride_xy, ::stride_z]
        else:
            m, xp, yp, zp = mask, X, Y, Z
        i, j, k = np.where(m)
        if i.size:
            ax.scatter(xp[i,j,k], yp[i,j,k], zp[i,j,k], s=size, c=COLORS[key],
                       alpha=alpha, marker="o", linewidths=0, depthshade=False)
    row = find_row_at_time(folder / "p1h_history.csv", snap["time_ps"])
    ax.set_title(f"{title}\nt={snap['time_ps']:.3f} ps · skeleton Tmax={float(row['Tmax_skeleton_K']):.0f} K · liquid Ag={float(row['liquid_Ag_volume_nm3']):.1f} nm³", fontsize=10, pad=8)
    ax.set_xlabel("x (nm)", labelpad=-2); ax.set_ylabel("y (nm)", labelpad=-2); ax.set_zlabel("z (nm)", labelpad=-3)
    if crop:
        i0 = int(np.argmin(np.abs(x-9.58832))); j0 = int(np.argmin(np.abs(y-9.65943950771472)))
        iz = np.where(a[0,i0,j0,:] > .5)[0]
        top = z[int(iz.max())] + snap["dx_m"]*1e9/2 if iz.size else max(z)
        ax.set_xlim(3.58832,15.58832); ax.set_ylim(3.65944,15.65944); ax.set_zlim(top-6.,top+snap["dx_m"]*1e9)
        ax.set_box_aspect((12,12,6.3)); ax.view_init(elev=32,azim=-56)
    else:
        ax.set_xlim(min(x),max(x)); ax.set_ylim(min(y),max(y)); ax.set_zlim(min(z),max(z))
        ax.set_box_aspect((len(x),len(y),len(z))); ax.view_init(elev=24,azim=-56)
    ax.tick_params(labelsize=7,pad=0); ax.grid(False)
    ax.xaxis.pane.set_alpha(.035); ax.yaxis.pane.set_alpha(.035); ax.zaxis.pane.set_alpha(.035)


def draw_cluster_case(ax, case, title, secondary_only=False):
    folder = HERE / "abc_full" / case
    paths = sorted(folder.glob("p1h_marangoni_on_hold_start_*_t*.bin"))
    snap = module.read_snapshot(paths[0])
    a = snap["fields"]
    x = (snap["x0"] + (np.arange(snap["nx"]) + .5)*snap["dx_m"])*1e9
    y = (snap["y0"] + (np.arange(snap["ny"]) + .5)*snap["dx_m"])*1e9
    z = (snap["z0"] + (np.arange(snap["nz"]) + .5)*snap["dx_m"])*1e9
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    sample, ag, sc, lf, code = a[0], a[2], a[3], a[4], a[19]
    theta = sample.astype(float)*ag.astype(float)*lf.astype(float)
    liquid_mask = theta >= .5
    labels, count = module.label_periodic_xy(liquid_mask)
    dx_nm = snap["dx_m"]*1e9
    cell_nm3 = dx_nm**3
    volumes = np.bincount(labels.ravel(), weights=theta.ravel(), minlength=count+1)*cell_nm3
    ranked = sorted(range(1, count+1), key=lambda lab: volumes[lab], reverse=True)
    # Limit context to the heated top region, but retain every thresholded
    # liquid component inside it. Components are colored independently.
    i0, j0 = int(np.argmin(np.abs(x-9.58832))), int(np.argmin(np.abs(y-9.65943950771472)))
    iz = np.where(sample[i0,j0,:] > .5)[0]
    top = z[int(iz.max())] + dx_nm/2 if iz.size else max(z)
    roi = ((np.abs(X-9.58832) <= 6.) & (np.abs(Y-9.65943950771472) <= 6.) &
           (Z >= min(z)) & (Z <= top+dx_nm))
    shown_labels = ranked[1:] if secondary_only else ranked
    shown_mask = np.isin(labels, shown_labels)
    near_pool = ndimage.binary_dilation(shown_mask,
                 structure=ndimage.generate_binary_structure(3,1),iterations=1)
    skeleton = ((code > 1.5) | (sc > .5)) & near_pool & roi
    solid = (code == 1) & (ag > .2) & (lf <= 1.e-3) & (Z >= top-.45) & roi
    partial = (code == 1) & (ag > .2) & (lf > 1.e-3) & (lf < .5) & roi
    context_layers = [(skeleton, COLORS["skeleton"], 5., .28, 3, 2)]
    if not secondary_only:
        context_layers += [
            (solid, COLORS["solid Ag"], 5., .08, 4, 2),
            (partial, COLORS["partial Ag"], 12., .62, 1, 1),
        ]
    for mask, color, size, alpha, sx, sz in context_layers:
        m = mask[::sx,::sx,::sz]
        xp,yp,zp = X[::sx,::sx,::sz],Y[::sx,::sx,::sz],Z[::sx,::sx,::sz]
        ii,jj,kk = np.where(m)
        if ii.size:
            ax.scatter(xp[ii,jj,kk],yp[ii,jj,kk],zp[ii,jj,kk],s=size,c=color,
                       alpha=alpha,linewidths=0,depthshade=False)
    palette = plt.get_cmap("tab20")
    for rank, lab in enumerate(ranked):
        if secondary_only and rank == 0:
            continue
        m = (labels == lab) & roi
        ii,jj,kk = np.where(m)
        if not ii.size:
            continue
        color = COLORS["liquid Ag"] if rank == 0 else palette((rank-1)%20)
        ax.scatter(X[ii,jj,kk],Y[ii,jj,kk],Z[ii,jj,kk],s=22,
                   c=[color],alpha=.22 if rank == 0 else .98,linewidths=0,depthshade=False)
    if secondary_only:
        secondary_volume = sum(volumes[lab] for lab in ranked[1:])
        label = f"{count-1} secondary pools · combined={secondary_volume:.1f} nm³ · main omitted"
    else:
        label = f"{count} six-face periodic pools · main={volumes[ranked[0]]:.1f} nm³"
    ax.set_title(f"{title}\n{label}",fontsize=10,pad=8)
    ax.set_xlabel("x (nm)",labelpad=-2); ax.set_ylabel("y (nm)",labelpad=-2); ax.set_zlabel("z (nm)",labelpad=-3)
    ax.set_xlim(3.58832,15.58832); ax.set_ylim(3.65944,15.65944); ax.set_zlim(min(z),top+dx_nm)
    ax.set_box_aspect((12,12,16)); ax.view_init(elev=24,azim=-56)
    ax.tick_params(labelsize=7,pad=0); ax.grid(False)
    ax.xaxis.pane.set_alpha(.035); ax.yaxis.pane.set_alpha(.035); ax.zaxis.pane.set_alpha(.035)


def main():
    legend = [Line2D([0],[0],marker="o",color="w",label=k,markerfacecolor=v,markersize=9) for k,v in COLORS.items()]
    for crop, name, suptitle in ((False,"P1P_ABC_FULL_HOLD_START_3D.png","A/B/C · heat-source region and connected melt pools"),
                                 (True,"P1P_ABC_FULL_HOLD_START_ZOOM.png","A/B/C · top-surface melt-pool detail")):
        fig = plt.figure(figsize=(16,7.4),facecolor="white")
        for i,(case,title,_) in enumerate(CASES,1):
            draw_case(fig.add_subplot(1,3,i,projection="3d"),case,title,crop)
        fig.legend(handles=legend,loc="upper center",ncol=4,frameon=False,bbox_to_anchor=(.5,1.01),fontsize=11)
        fig.suptitle("Ag–Ti₃SiC₂ · full-run hold-start fields (near peak pool)",y=.96,fontsize=15)
        fig.text(.5,.015,"Exploratory continuum VOF · liquid Ag orange · partial Ag yellow · solid Ag blue · Ti₃SiC₂ dark",ha="center",fontsize=9,color="#555")
        fig.tight_layout(rect=[0,.035,1,.91],w_pad=.5)
        out=HERE/name; fig.savefig(out,dpi=180,bbox_inches="tight"); plt.close(fig); print(out)
    for secondary_only, name, heading, note in (
        (False,"P1P_ABC_LIQUID_COMPONENTS.png","individual liquid-pool components","Main liquid component is translucent orange to reveal smaller pools, each shown in a distinct color."),
        (True,"P1P_ABC_SECONDARY_POOLS.png","secondary liquid pools only","The dominant component is omitted so the smaller disconnected pools are visible; color identifies each pool."),
    ):
        fig = plt.figure(figsize=(16,7.2),facecolor="white")
        for i,(case,title,_) in enumerate(CASES,1):
            draw_cluster_case(fig.add_subplot(1,3,i,projection="3d"),case,title,secondary_only)
        fig.suptitle(f"Ag–Ti₃SiC₂ · {heading} at the peak/hold-start snapshot",y=.97,fontsize=15)
        fig.text(.5,.02,note+" Dark = Ti₃SiC₂; pale blue = solid Ag; yellow = partial Ag.",ha="center",fontsize=9,color="#555")
        fig.tight_layout(rect=[0,.04,1,.92],w_pad=.5)
        out=HERE/name; fig.savefig(out,dpi=180,bbox_inches="tight"); plt.close(fig); print(out)


if __name__ == "__main__": main()
