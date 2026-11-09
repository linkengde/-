#!/usr/bin/env python3
"""Corrected radial-core reanalysis of the Round7 1.6 ps TTMD trajectory."""
from __future__ import annotations

import csv
import gzip
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Patch
import numpy as np
from ovito.io import import_file
from ovito.modifiers import PolyhedralTemplateMatchingModifier, SelectTypeModifier, SteinhardtParametersModifier
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import ConvexHull, QhullError, cKDTree

OUT = Path(os.environ.get("ROUND7_EXTRACT_OUT", str(Path(__file__).resolve().parent))).resolve()
MODEL = OUT.parents[1]
RUN_ORIGINAL = MODEL / "runs" / "ttm_mtp_round7_continue_to_1p5ps_20261109"
RUN = Path(os.environ.get("ROUND7_RUN_DIR", str(RUN_ORIGINAL))).resolve()
RUN_05 = MODEL / "runs" / "ttm_mtp_round7_continue_to_0p5ps_20261109"
END_PS = float(os.environ.get("ROUND7_END_PS", "0.90"))
WINDOW_PS = float(os.environ.get("ROUND7_WINDOW_PS", "0.4"))
SAMPLE_INTERVAL_PS = float(os.environ.get("ROUND7_SAMPLE_INTERVAL_PS", "0.05"))
TIMES = np.round(np.arange(END_PS - WINDOW_PS, END_PS + SAMPLE_INTERVAL_PS * 0.25, SAMPLE_INTERVAL_PS), 5)
END_TAG = f"{END_PS:.2f}".rstrip("0").rstrip(".").replace(".", "p")
BOX_XY = 200.0
Q6_CUTOFF_A = 4.0
Q6_SCREEN = 0.4081069282833649
PTM_FCC = 1
PTM_NONFCC_FRACTION_SCREEN = 0.8
MSD_1PS_REFERENCE_A2 = 0.5171679303519052
MSD_SCREEN_A2 = MSD_1PS_REFERENCE_A2 * WINDOW_PS
CONNECTIVITY_A = 3.5
font_manager.fontManager.addfont("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
plt.rcParams["font.family"] = "Noto Sans CJK SC"


def load_initial():
    path = MODEL / "ag_ti3sic2_qsgs.data.gz"
    with gzip.open(path, "rt") as stream:
        lines = stream.readlines()
    natoms = int(next(line.split()[0] for line in lines if "atoms" in line))
    start = next(i for i, line in enumerate(lines) if line.startswith("Atoms #")) + 2
    atoms = np.loadtxt(lines[start:start + natoms])
    atoms = atoms[np.argsort(atoms[:, 0])]
    ids = atoms[:, 0].astype(np.int64)
    types = atoms[:, 1].astype(np.int32)
    xyz = atoms[:, 2:5]
    newton = (xyz[:, 2] >= 40.0) & (xyz[:, 2] < 150.0)
    dx = (xyz[:, 0] - 100.0 + BOX_XY / 2) % BOX_XY - BOX_XY / 2
    dy = (xyz[:, 1] - 100.0 + BOX_XY / 2) % BOX_XY - BOX_XY / 2
    depth = np.maximum(0.0, 150.0 - xyz[:, 2])
    # Keep the historical exponent cutoff (<= 2) but use the full radial x/y term.
    core = newton & (types == 4) & ((dx / 30.0) ** 2 + (dy / 30.0) ** 2 + depth / 20.0 <= 2.0)
    return ids, types, newton, core


def read_frame(time_ps: float):
    tag = f"{time_ps:.2f}".rstrip("0").rstrip(".").replace(".", "p")
    candidates = [RUN_05] if time_ps <= 0.500001 else [RUN, RUN_ORIGINAL]
    path = next((run / "checkpoints" / f"atoms_{tag}ps.dump" for run in candidates
                 if (run / "checkpoints" / f"atoms_{tag}ps.dump").is_file()), None)
    if path is None:
        raise FileNotFoundError(f"No saved snapshot for {time_ps:.5f} ps in {[str(x) for x in candidates]}")
    data = np.loadtxt(path, skiprows=9)
    data = data[np.argsort(data[:, 0])]
    return path, data[:, 0].astype(np.int64), data[:, 1].astype(np.int32), data[:, 2:5]


def q6_ptm(path: Path):
    pipeline = import_file(str(path))
    pipeline.modifiers.append(SelectTypeModifier(types={4}))
    q6 = SteinhardtParametersModifier()
    q6.max_l = 6
    q6.neighbor_mode = SteinhardtParametersModifier.NeighborMode.Cutoff
    q6.cutoff = Q6_CUTOFF_A
    q6.only_selected = True
    pipeline.modifiers.append(q6)
    ptm = PolyhedralTemplateMatchingModifier()
    ptm.only_selected = True
    pipeline.modifiers.append(ptm)
    particles = pipeline.compute().particles
    order = np.argsort(np.asarray(particles["Particle Identifier"], dtype=np.int64))
    ids = np.asarray(particles["Particle Identifier"], dtype=np.int64)[order]
    q6_values = np.asarray(particles["Steinhardt"], dtype=float)[order, 4]
    ptm_values = np.asarray(particles["Structure Type"], dtype=np.int32)[order]
    return ids, q6_values, ptm_values


def drift_corrected_displacement(frames, newton):
    net = np.zeros_like(frames[0])
    for before, after in zip(frames[:-1], frames[1:]):
        delta = after - before
        delta[:, :2] -= np.rint(delta[:, :2] / BOX_XY) * BOX_XY
        delta -= np.mean(delta[newton], axis=0)
        net += delta
    return np.einsum("ij,ij->i", net, net)


def clusters(mask, xyz, ids):
    selected = np.flatnonzero(mask)
    if not len(selected):
        return [], {"candidate_atoms": 0, "component_count": 0, "largest_component_atoms": 0,
                    "largest_component_dimensions_A": None, "largest_component_convex_hull_volume_A3": None}
    points = xyz[selected].copy()
    points[:, 0] %= BOX_XY
    points[:, 1] %= BOX_XY
    zlo = float(points[:, 2].min())
    zspan = float(np.ptp(points[:, 2]))
    zperiod = max(zspan + 2 * CONNECTIVITY_A + 1.0e-6, 4 * CONNECTIVITY_A)
    points[:, 2] = points[:, 2] - zlo + CONNECTIVITY_A
    pairs = cKDTree(points, boxsize=[BOX_XY, BOX_XY, zperiod]).query_pairs(CONNECTIVITY_A, output_type="ndarray")
    if len(pairs):
        graph = coo_matrix((np.ones(2 * len(pairs)),
                            (np.r_[pairs[:, 0], pairs[:, 1]], np.r_[pairs[:, 1], pairs[:, 0]])),
                           shape=(len(selected), len(selected))).tocsr()
        _, labels = connected_components(graph, directed=False, return_labels=True)
    else:
        labels = np.arange(len(selected))
    counts = np.bincount(labels)
    order = np.argsort(counts)[::-1]
    clusters_out = []
    for label in order:
        local = np.flatnonzero(labels == label)
        original = selected[local]
        component = points[local].copy()
        for axis in (0, 1):
            wrapped = component[:, axis] % BOX_XY
            ordered = np.sort(wrapped)
            gaps = np.diff(np.r_[ordered, ordered[0] + BOX_XY])
            gap_index = int(np.argmax(gaps))
            arc_start = ordered[(gap_index + 1) % len(ordered)]
            component[:, axis] = arc_start + (wrapped - arc_start) % BOX_XY
        hull_volume = None
        try:
            if len(component) >= 4:
                hull_volume = float(ConvexHull(component).volume)
        except QhullError:
            pass
        clusters_out.append({
            "rank": len(clusters_out) + 1,
            "atoms": int(len(original)),
            "atom_ids": np.asarray(ids)[original].astype(int).tolist(),
            "dimensions_A": np.ptp(component, axis=0).tolist(),
            "convex_hull_volume_A3": hull_volume,
        })
    largest = clusters_out[0]
    summary = {"candidate_atoms": int(len(selected)), "component_count": int(len(counts)),
               "largest_component_atoms": largest["atoms"],
               "largest_component_dimensions_A": largest["dimensions_A"],
               "largest_component_convex_hull_volume_A3": largest["convex_hull_volume_A3"],
               "connectivity_cutoff_A": CONNECTIVITY_A, "periodic_xy": True}
    return clusters_out, summary


def write_dump(path: Path, timestep: int, ids, xyz, q6, ptm):
    with path.open("w", encoding="utf-8") as stream:
        stream.write(f"ITEM: TIMESTEP\n{timestep}\nITEM: NUMBER OF ATOMS\n{len(ids)}\n")
        stream.write("ITEM: BOX BOUNDS pp pp ff\n0 200\n0 200\n0 250\n")
        stream.write("ITEM: ATOMS id type x y z q6 ptm\n")
        for atom_id, point, q, p in zip(ids, xyz, q6, ptm):
            stream.write(f"{int(atom_id)} 4 {point[0]:.8f} {point[1]:.8f} {point[2]:.8f} {q:.8f} {int(p)}\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ids0, types0, newton, core = load_initial()
    frames, q6_frames, ptm_frames = [], [], []
    for n, t in enumerate(TIMES, start=1):
        path, ids, types, xyz = read_frame(float(t))
        if not np.array_equal(ids, ids0) or not np.array_equal(types, types0):
            raise RuntimeError(f"Atom identities differ in {path}")
        qids, q6, ptm = q6_ptm(path)
        if not np.array_equal(qids, ids0):
            raise RuntimeError(f"OVITO ID order differs in {path}")
        frames.append(xyz)
        q6_frames.append(q6)
        ptm_frames.append(ptm)
        print(f"processed Q6/PTM snapshot {n}/{len(TIMES)} at {t:.2f} ps", flush=True)
    q6_mean = np.mean(q6_frames, axis=0)
    nonfcc = np.mean(np.asarray(ptm_frames) != PTM_FCC, axis=0)
    msd_i = drift_corrected_displacement(frames, newton)
    joint = core & (q6_mean < Q6_SCREEN) & (nonfcc >= PTM_NONFCC_FRACTION_SCREEN) & (msd_i >= MSD_SCREEN_A2)
    endpoint = frames[-1]
    endpoint_q6 = q6_frames[-1]
    endpoint_ptm = ptm_frames[-1]
    clusters_json, cluster_summary = clusters(joint, endpoint, ids0)
    candidate_idx = np.flatnonzero(joint)
    largest_ids = set(clusters_json[0]["atom_ids"]) if clusters_json else set()
    largest_idx = np.array([i for i in candidate_idx if int(ids0[i]) in largest_ids], dtype=int)

    report = {
        "status": "PROVISIONAL_ROUND7_Q6_PTM_MSD_SCREEN",
        "source_run": str(RUN), "snapshot_endpoint_ps": float(END_PS), "window_ps": [float(TIMES[0]), float(TIMES[-1])],
        "window_duration_ps": WINDOW_PS, "sample_interval_ps": SAMPLE_INTERVAL_PS,
        "potential": "Round7 energy-force pilot MTP; not independently validated for liquid Ag/interface behavior",
        "same_atoms": "initial Newton-layer Ag atoms selected by the legacy x-z spatial prefilter (dy omitted; not a radial source core)",
        "spatial_prefilter_note": "The legacy predicate computes dy but does not use it; selected atoms span y wherever the x-z condition passes.",
        "source_core_ag_atoms": int(core.sum()), "joint_candidate_atoms": int(joint.sum()),
        "rule": {
            "q6": f"mean per-atom Ag q6 over this window < {Q6_SCREEN:.9f}; cutoff {Q6_CUTOFF_A} A",
            "ptm": f"non-FCC fraction over snapshots >= {PTM_NONFCC_FRACTION_SCREEN:.2f}; FCC code {PTM_FCC}",
            "msd_i": f"drift-corrected per-atom net displacement squared over {WINDOW_PS:.2f} ps >= {MSD_SCREEN_A2:.9f} A^2",
            "msd_threshold_note": f"linearly scaled from the Round5 1 ps exploratory threshold {MSD_1PS_REFERENCE_A2:.9f} A^2; not Round7 calibrated",
        },
        "cluster_summary": cluster_summary, "clusters_ranked": clusters_json,
        "threshold_provenance": "Q6/PTM/MSD numerical cutoffs borrowed from Round5 pure-Ag development-MTP exploratory references; scaling and transfer to Round7 are provisional.",
        "limitations": [
            (f"This {WINDOW_PS:.2f} ps screen is shorter than the reference 1 ps window; it is not a final melt-pool determination."
             if WINDOW_PS < 0.999 else
             "The window spans 1 ps, but the thresholds still come from Round5 and are not Round7-calibrated."),
            "Round7 has a different energy-force MTP and its own pure-Ag thresholds have not been calibrated.",
            "Q6/PTM/MSD jointly identify mobile non-FCC Ag candidates, but do not validate the underlying interfacial potential.",
            "Clusters use 3.5 A Ag-Ag connectivity, periodic x/y and open z; convex-hull volume is only a geometric envelope.",
        ],
    }
    (OUT / f"round7_{END_TAG}ps_joint_screen.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    with (OUT / f"round7_{END_TAG}ps_source_core_ag.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        window_tag = f"{WINDOW_PS:g}ps"
        writer.writerow(["atom_id", f"mean_q6_{window_tag}", f"ptm_nonfcc_fraction_{window_tag}", f"msd_i_{window_tag}_A2", "joint_candidate"])
        for i in np.flatnonzero(core):
            writer.writerow([int(ids0[i]), float(q6_mean[i]), float(nonfcc[i]), float(msd_i[i]), int(joint[i])])
    timestep = int(round(END_PS * 10000))
    write_dump(OUT / f"round7_{END_TAG}ps_joint_candidates.dump", timestep, ids0[candidate_idx], endpoint[candidate_idx], endpoint_q6[candidate_idx], endpoint_ptm[candidate_idx])
    write_dump(OUT / f"round7_{END_TAG}ps_largest_candidate_cluster.dump", timestep, ids0[largest_idx], endpoint[largest_idx], endpoint_q6[largest_idx], endpoint_ptm[largest_idx])

    colors = np.full(len(ids0), "#8ABBD2", dtype=object)
    skeleton = np.isin(types0, (1, 2, 3))
    colors[skeleton] = "#263F50"
    colors[core & ~joint] = "#D5DDE0"
    colors[joint] = "#F28C28"
    fig = plt.figure(figsize=(14, 10), dpi=160, facecolor="white")
    ax = fig.add_axes([0.07, 0.04, 0.87, 0.80], projection="3d")
    ax.scatter(endpoint[:, 0], endpoint[:, 1], endpoint[:, 2], c=colors, s=2.2, alpha=0.9,
               linewidths=0, depthshade=False, rasterized=True)
    corners = np.array([[0,0,0],[200,0,0],[200,200,0],[0,200,0],[0,0,150],[200,0,150],[200,200,150],[0,200,150]], dtype=float)
    for a, b in [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]:
        ax.plot(*corners[[a,b]].T, color="#AFB9BF", linewidth=0.65, alpha=0.65)
    ax.set_xlim(0,200); ax.set_ylim(0,200); ax.set_zlim(0,200)
    ax.set_box_aspect((200,200,200)); ax.view_init(elev=27, azim=-58); ax.set_axis_off()
    handles = [Patch(facecolor="#263F50", label=f"Ti3SiC2 骨架 ({int(skeleton.sum()):,})"),
               Patch(facecolor="#8ABBD2", label=f"Ag ({int((types0 == 4).sum()):,})"),
               Patch(facecolor="#D5DDE0", label=f"热源核心区 Ag，未通过筛选 ({int((core & ~joint).sum()):,})"),
               Patch(facecolor="#F28C28", label=f"Q6+PTM+MSD 临时候选 ({int(joint.sum()):,})")]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.025,0.89), frameon=False, fontsize=10)
    fig.text(0.5, 0.965, f"Round7 Ag–Ti3SiC2 三维快照 | t = {END_PS:.3f} ps", ha="center", va="top", fontsize=19)
    fig.text(0.5, 0.925, f"{TIMES[0]:.2f}–{END_PS:.2f} ps 临时联合筛查 | 候选 {int(joint.sum()):,} Ag | 最大团簇 {cluster_summary['largest_component_atoms']:,}",
             ha="center", va="top", fontsize=13)
    fig.text(0.5, 0.012, f"仅为{WINDOW_PS:.1f} ps窗口的探索性筛查；阈值借用Round5，Round7尚未独立标定，不能视作最终熔化确认。",
             ha="center", va="bottom", fontsize=9, color="#555555")
    fig.savefig(OUT / f"round7_{END_TAG}ps_joint_candidate_3d.png", dpi=160, facecolor="white")
    plt.close(fig)
    print(json.dumps({"joint_candidates": int(joint.sum()), "clusters": cluster_summary["component_count"],
                      "largest_cluster": cluster_summary["largest_component_atoms"],
                      "dimensions_A": cluster_summary["largest_component_dimensions_A"],
                      "image": str(OUT / f"round7_{END_TAG}ps_joint_candidate_3d.png")}, indent=2), flush=True)


if __name__ == "__main__":
    main()
