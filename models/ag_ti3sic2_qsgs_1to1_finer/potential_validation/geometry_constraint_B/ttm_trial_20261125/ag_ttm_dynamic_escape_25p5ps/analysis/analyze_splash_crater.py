#!/usr/bin/env python3
"""对动态脱离分支的 Ag 喷溅候选团簇和表面凹陷做同口径筛查。"""
import gzip
import glob
import json
import os
import numpy as np
from scipy.spatial import cKDTree

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.dirname(__file__)


def analyze(path):
    """读取一帧 LAMMPS dump，统计高于初始表面的脱离 Ag 团簇。"""
    with gzip.open(path, "rt") as stream:
        while not stream.readline().startswith("ITEM: ATOMS"):
            pass
        atoms = np.loadtxt(stream)

    # dump 列为 id,type,x,y,z,vx,vy,vz,fx,fy,fz；type 4 是 Ag。
    is_ag = atoms[:, 1].astype(int) == 4
    xyz = atoms[is_ag, 2:5].copy()
    ag_vz = atoms[is_ag, 7].copy()
    xyz[:, 0] %= 200.0
    xyz[:, 1] %= 200.0

    # 对全体 Ag 以 4 Å 邻接构造团簇。x/y 周期；z 非周期，模拟中 z 最大值远小于 300 Å。
    tree = cKDTree(xyz, boxsize=[200.0, 200.0, np.inf])
    pairs = tree.query_pairs(4.0, output_type="ndarray")
    parent = np.arange(len(xyz))

    def root(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    for first, second in pairs:
        first_root = root(int(first))
        second_root = root(int(second))
        if first_root != second_root:
            parent[second_root] = first_root

    labels = np.array([root(i) for i in range(len(xyz))])
    _, inverse, sizes = np.unique(labels, return_inverse=True, return_counts=True)
    minimum_z = np.full(len(sizes), np.inf)
    np.minimum.at(minimum_z, inverse, xyz[:, 2])
    # 仅将整个团簇都高于 z=150 Å 的连通体记作“脱离表面团簇候选”。
    detached = minimum_z > 150.0
    detached_sizes = sorted((int(n) for n in sizes[detached]), reverse=True)
    isolated = detached[inverse] & (sizes[inverse] == 1)

    # 取 100<z<160 Å 的 Ag，避免空中原子抬高表面；按 10 Å XY 格算 z90。
    surface = xyz[(xyz[:, 2] > 100.0) & (xyz[:, 2] < 160.0)]
    bins = {}
    for x, y, z in surface:
        bins.setdefault((int(x // 10), int(y // 10)), []).append(float(z))
    radial = []
    for (ix, iy), values in bins.items():
        radius = np.hypot((ix + 0.5) * 10.0 - 100.0,
                          (iy + 0.5) * 10.0 - 100.0)
        radial.append((radius, float(np.percentile(values, 90))))
    radial = np.asarray(radial)
    heights = []
    for low, high in ((0, 20), (20, 40), (40, 60), (60, 80)):
        region = radial[(radial[:, 0] >= low) & (radial[:, 0] < high), 1]
        heights.append(float(np.median(region)) if len(region) else None)

    return {
        "frame": os.path.basename(path),
        "Ag_atom_count": int(is_ag.sum()),
        "Ag_z_max_A": round(float(xyz[:, 2].max()), 3),
        "Ag_count_z_gt_150_A": int((xyz[:, 2] > 150.0).sum()),
        "Ag_count_z_gt_160_A": int((xyz[:, 2] > 160.0).sum()),
        "fully_detached_cluster_count_min_z_gt_150_A": int(detached.sum()),
        "fully_detached_cluster_atom_count": int(sizes[detached].sum()),
        "largest_detached_clusters_atoms": detached_sizes[:12],
        "detached_single_atom_count": int(np.count_nonzero(isolated)),
        "detached_single_atoms_with_positive_vz": int(np.count_nonzero(ag_vz[isolated] > 0.0)),
        "detached_clusters_size_le_5_count": int(np.count_nonzero((sizes[detached] <= 5))),
        "detached_single_atoms": [{"id": int(atom_id), "z_A": round(float(z), 3), "vz_A_per_ps": round(float(vz), 4)} for atom_id, z, vz in zip(atoms[is_ag, 0].astype(int)[isolated], xyz[isolated, 2], ag_vz[isolated])],
        "surface_z90_median_by_radius_A": {
            "r0_20": round(heights[0], 3),
            "r20_40": round(heights[1], 3),
            "r40_60": round(heights[2], 3),
            "r60_80": round(heights[3], 3),
        },
        "center_minus_r20_40_A": round(heights[0] - heights[1], 3),
    }


patterns = (
    "continued_9p5_to_10p0ps/frame*",
    "ttm_dynamic_escape_11p5_to_12p0ps/frame*",
    "ttm_dynamic_escape_12p5_to_13p0ps/frame*",
    "ttm_dynamic_escape_13p0_to_13p5ps/frame*",
    "ttm_dynamic_escape_13p5_to_14p0ps/frame*",
    "ttm_dynamic_escape_14p0_to_14p5ps/frame*",
    "ttm_dynamic_escape_14p5_to_15p0ps/frame*",
    "ttm_dynamic_escape_15p0_to_15p5ps/frame*",
    "ttm_dynamic_escape_15p5_to_16p0ps/frame*",
    "ttm_dynamic_escape_16p0_to_16p5ps/frame*",
    "ttm_dynamic_escape_16p5_to_17p0ps/frame*",
    "ttm_dynamic_escape_17p0_to_17p5ps/frame*",
    "ttm_dynamic_escape_17p5_to_18p0ps/frame*",
    "ttm_dynamic_escape_18p0_to_18p5ps/frame*",
    "ttm_dynamic_escape_18p5_to_19p0ps/frame*",
    "ttm_dynamic_escape_19p0_to_19p5ps/frame*",
    "ttm_dynamic_escape_19p5_to_20p0ps/frame*",
    "ttm_dynamic_escape_20p0_to_20p5ps/frame*",
    "ttm_dynamic_escape_20p5_to_21p0ps/frame*",
    "ttm_dynamic_escape_21p0_to_21p5ps/frame*",
    "ttm_dynamic_escape_extend_z400_21p5_to_22p0ps/frame*",
    "ttm_dynamic_escape_extend_z400_22p0_to_22p5ps/frame*",
    "ttm_dynamic_escape_extend_z400_22p5_to_23p0ps/frame*",
    "ttm_dynamic_escape_extend_z400_23p0_to_23p5ps/frame*",
    "ttm_dynamic_escape_extend_z400_23p5_to_24p0ps/frame*",
    "ttm_dynamic_escape_extend_z400_24p0_to_24p5ps/frame*",
    "ttm_dynamic_escape_extend_z400_24p5_to_25p0ps/frame*",
    "ttm_dynamic_escape_extend_z400_25p0_to_25p5ps/frame*",
)
paths = []
for pattern in patterns:
    paths.extend(glob.glob(os.path.join(BASE, pattern)))
results = [analyze(path) for path in paths]
with open(os.path.join(OUT, "splash_crater_checkpoint_analysis.json"), "w") as stream:
    json.dump(results, stream, indent=2)
for row in results:
    print(json.dumps(row, ensure_ascii=False))
