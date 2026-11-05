#!/usr/bin/env python3
"""Plot saved P1H free-surface height changes at four existing event snapshots."""
from __future__ import annotations

import argparse
import csv
import gzip
import struct
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager

HEADER = struct.Struct("<8s6I5d")
SOURCE_X_NM = 9.58832
SOURCE_Y_NM = 9.65943950771472
RUN_REL = "pressure_sensitivity_pi1/marangoni_on_gaussian_167MPa"
EVENTS = [
    ("压力首次实际作用后（场输出晚 0.116 ps）", "regular_t000.263509"),
    ("最大熔池 / 关热时", "heat_off_t000.797446"),
    ("压力已关闭后的首个保存场（晚 0.197 ps）", "regular_t001.014300"),
    ("Ag 完全凝固", "solidified_t001.112870"),
]


def read_frame(path: Path, magic_expected: bytes, fields_expected: int):
    with gzip.open(path, "rb") as stream:
        raw = stream.read()
    magic, version, nx, ny, nz, nf, flags, time_ps, dx_m, x0, y0, z0 = HEADER.unpack_from(raw)
    if magic != magic_expected or version != 1 or nf != fields_expected:
        raise ValueError(f"invalid P1H frame: {path}")
    values = np.frombuffer(raw, dtype="<f4", offset=HEADER.size)
    expected = nx*ny*nz*nf
    if values.size != expected:
        raise ValueError(f"truncated P1H frame: {path}")
    return values.reshape(nf,nx,ny,nz), dict(time_ps=time_ps, dx_m=dx_m,
                                             origin_m=(x0,y0,z0),shape=(nx,ny,nz))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-root", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    for candidate in ("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
                      "/usr/share/fonts/truetype/noto/NotoSansCJK-Regular.ttc"):
        if Path(candidate).exists():
            font_manager.fontManager.addfont(candidate)
            plt.rcParams["font.family"] = font_manager.FontProperties(fname=candidate).get_name()
            break
    plt.rcParams["axes.unicode_minus"] = False

    run_dir = args.run_root/RUN_REL
    initial, imeta = read_frame(run_dir/"p1h_marangoni_on_initial_t000.000000_ps.bin.gz",
                               b"P1HSNAP1",20)
    initial_vof = initial[0]
    fig, axes = plt.subplots(2,2,figsize=(14,11),constrained_layout=True)
    cmap = plt.get_cmap("RdBu_r")
    last_im = None
    report_rows = {r["time_ps"]:r for r in csv.DictReader((run_dir/"p1h_history.csv").open())}
    for ax,(title,suffix) in zip(axes.flat,EVENTS):
        main_path=run_dir/f"p1h_marangoni_on_{suffix}_ps.bin.gz"
        pressure_path=run_dir/f"p1h_pressure_marangoni_on_{suffix}_ps.bin.gz"
        fields,meta=read_frame(main_path,b"P1HSNAP1",20)
        pfields,pmeta=read_frame(pressure_path,b"P1HPRS01",4)
        if meta["time_ps"]!=pmeta["time_ps"] or meta["shape"]!=imeta["shape"]:
            raise ValueError(f"main/pressure/initial snapshot mismatch: {suffix}")
        nx,ny,nz=meta["shape"]; dx_nm=meta["dx_m"]*1e9
        x0,y0,z0=(v*1e9 for v in meta["origin_m"])
        zface=z0+(np.arange(nz)+1.)*dx_nm
        h0=np.max(np.where(initial_vof>.5,zface[None,None,:],-np.inf),axis=2)
        h=np.max(np.where(fields[0]>.5,zface[None,None,:],-np.inf),axis=2)
        delta_h=h0-h
        xc=x0+(np.arange(nx)+.5)*dx_nm; yc=y0+(np.arange(ny)+.5)*dx_nm
        extent=(xc[0]-.5*dx_nm,xc[-1]+.5*dx_nm,yc[0]-.5*dx_nm,yc[-1]+.5*dx_nm)
        last_im=ax.imshow(delta_h.T,origin="lower",extent=extent,cmap=cmap,
                          vmin=-.65,vmax=.65,interpolation="nearest",aspect="equal")
        pressure_xy=np.max(pfields[2],axis=2)/1.e6
        contour_levels=[10,25,50,100,150]
        levels=[v for v in contour_levels if pressure_xy.max()>v]
        if levels:
            ax.contour(xc,yc,pressure_xy.T,levels=levels,colors="#222222",
                       linewidths=.75,alpha=.75)
        # The center circle indicates the analysis region used in the CSV.
        circle=plt.Circle((SOURCE_X_NM,SOURCE_Y_NM),5.,fill=False,
                          color="#111111",linewidth=.7,linestyle="--")
        ax.add_patch(circle)
        pmax=float(pfields[2].max())/1.e6
        trow=min(report_rows.values(),key=lambda r:abs(float(r["time_ps"])-meta["time_ps"]))
        ax.set_title(f"{title}\n保存场 t={meta['time_ps']:.6f} ps | p_eff,max={pmax:.3g} MPa | "
                     f"液态 Ag={float(trow['liquid_Ag_volume_nm3']):.2f} nm³",fontsize=10)
        ax.set_xlabel("x (nm)"); ax.set_ylabel("y (nm)")
        ax.set_xlim(extent[0],extent[1]); ax.set_ylim(extent[2],extent[3])
    cb=fig.colorbar(last_im,ax=list(axes.flat),shrink=.83,pad=.025)
    cb.set_label("初始表面高度 − 当前表面高度 (nm)；正值为下凹，负值为上移")
    fig.suptitle("Ag70 P1H：现有 167 MPa 算例的自由表面响应\n"
                 "俯视高度差图；黑色虚线为热源中心半径 5 nm；压力等值线单位 MPa\n"
                 "固定网格投影（Δx≈0.302 nm），用于查看表面响应，不代表原子尺度形貌",
                 fontsize=14)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.output,dpi=180,facecolor="white")
    plt.close(fig)
    print(f"wrote {args.output}")


if __name__=="__main__":
    main()
