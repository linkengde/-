#!/usr/bin/env python3
"""统计初态和 5.5 ps 末帧中的 Ag-only 最近邻距离。"""

import gzip
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


HERE = Path(__file__).resolve().parent
EARLY = HERE.parent / "round7_347k_3ps"
LX = 200.0
LY = 200.0
# z 是开放边界；给 KDTree 一个远大于盒高的周期长度，避免 z 向镜像参与最近邻。
OPEN_Z_SENTINEL = 1.0e6


def load_ag_positions(path):
    """读取压缩 LAMMPS dump 的 Ag 坐标，并将周期 xy 坐标折回盒内。"""
    with gzip.open(path, "rt") as stream:
        # 轨迹头固定为 9 行，列顺序为 id type x y z vx vy vz fx fy fz。
        values = np.loadtxt(stream, skiprows=9, usecols=(0, 1, 2, 3, 4))
    ag = values[values[:, 1] == 4]
    positions = ag[:, 2:5].copy()
    positions[:, 0] %= LX
    positions[:, 1] %= LY
    return positions


def summarize(path):
    """用 xy 周期、z 开放的近似盒子计算每个 Ag 最近的 Ag 邻居。"""
    positions = load_ag_positions(path)
    tree = cKDTree(
        positions,
        boxsize=np.array([LX, LY, OPEN_Z_SENTINEL]),
    )
    distances, _ = tree.query(positions, k=2, workers=4)
    nearest = distances[:, 1]
    print(
        f"{path.name}: nAg={len(positions)} "
        f"min={nearest.min():.6f} Å "
        f"p01={np.percentile(nearest, 1):.6f} Å "
        f"p05={np.percentile(nearest, 5):.6f} Å "
        f"median={np.median(nearest):.6f} Å "
        f"count(nn<2.4 Å)={(nearest < 2.4).sum()}"
    )


if __name__ == "__main__":
    summarize(EARLY / "frame_00000_step0.dump.gz")
    summarize(HERE / "frame_05.500ps_step11000.dump.gz")
