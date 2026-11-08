#!/usr/bin/env python3
"""Apply the calibrated same-Ag Q6/PTM/MSD rule over rolling 1 ps windows."""
from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ovito.io import import_file
from ovito.modifiers import PolyhedralTemplateMatchingModifier, SelectTypeModifier, SteinhardtParametersModifier
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import ConvexHull, QhullError, cKDTree

HERE = Path(__file__).resolve().parent
MODEL = HERE.parents[1]
OLD_RUN = MODEL / "runs" / "ttm_mtp_round5_continue_to_0p5ps_20261107"
MID_RUN = MODEL / "runs" / "ttm_mtp_round5_joint_classifier_to_1p15ps_20261107"
THRESHOLD_FILE = MODEL.parent / "ag_ti3sic2_qsgs_1to1_finer" / "potential_validation" / "mtp_round5_joint_melt_classifier_20261107" / "calibrated_thresholds.json"
TIMES_PS = np.round(np.arange(0.15, 1.5501, 0.05), 5)
WINDOW_INTERVALS = 20
BOX_XY = 200.0
FCC_CODE = 1
Q6_CUTOFF_A = 4.0
CONNECTIVITY_A = 3.5


def initial_atoms():
    with gzip.open(MODEL / "ag_ti3sic2_qsgs.data.gz", "rt") as stream:
        lines = stream.readlines()
    natoms = int(next(line.split()[0] for line in lines if "atoms" in line))
    start = next(i for i, line in enumerate(lines) if line.startswith("Atoms #")) + 2
    rows = np.loadtxt(lines[start:start + natoms])
    rows = rows[np.argsort(rows[:, 0])]
    return rows[:, 0].astype(np.int64), rows[:, 1].astype(np.int32), rows[:, 2:5]


def dump_path(time_ps: float) -> Path:
    tag = f"{time_ps:.4f}".rstrip("0").rstrip(".").replace(".", "p")
    if time_ps <= 0.500001:
        folder = OLD_RUN
    elif time_ps <= 1.150001:
        folder = MID_RUN
    else:
        folder = HERE
    return folder / "checkpoints" / f"atoms_{tag}ps.dump"


def parse_one_frame(path: Path):
    with path.open() as stream:
        lines = stream.readlines()
    natoms = int(lines[3])
    bounds = np.array([[float(v) for v in line.split()[:2]] for line in lines[5:8]])
    headers = lines[8].split()[2:]
    col = {name: i for i, name in enumerate(headers)}
    data = np.loadtxt(lines[9:9 + natoms])
    data = data[np.argsort(data[:, col["id"]])]
    ids = data[:, col["id"]].astype(np.int64)
    types = data[:, col["type"]].astype(np.int32)
    xyz = data[:, [col[a] for a in ("x", "y", "z")]]
    return {"path": path, "ids": ids, "types": types, "xyz": xyz, "time_ps": int(lines[1]) * 0.0001}


def q6_ptm_at(path: Path):
    pipeline = import_file(str(path))
    pipeline.modifiers.append(SelectTypeModifier(types={4}))
    q6_mod = SteinhardtParametersModifier()
    q6_mod.max_l = 6
    q6_mod.neighbor_mode = SteinhardtParametersModifier.NeighborMode.Cutoff
    q6_mod.cutoff = Q6_CUTOFF_A
    q6_mod.only_selected = True
    pipeline.modifiers.append(q6_mod)
    ptm_mod = PolyhedralTemplateMatchingModifier()
    ptm_mod.only_selected = True
    pipeline.modifiers.append(ptm_mod)
    data = pipeline.compute()
    ids = np.asarray(data.particles["Particle Identifier"], dtype=np.int64)
    order = np.argsort(ids)
    types = np.asarray(data.particles["Particle Type"], dtype=np.int32)[order]
    q6 = np.asarray(data.particles["Steinhardt"], dtype=float)[order, 4]
    ptm = np.asarray(data.particles["Structure Type"], dtype=np.int32)[order]
    return ids[order], types, q6, ptm


def displacement_window(frames, newton_mask):
    net = np.zeros_like(frames[0]["xyz"])
    for previous, current in zip(frames[:-1], frames[1:]):
        delta = current["xyz"] - previous["xyz"]
        delta[:, :2] -= np.rint(delta[:, :2] / BOX_XY) * BOX_XY
        delta -= np.mean(delta[newton_mask], axis=0)
        net += delta
    return np.linalg.norm(net, axis=1)


def cluster_stats(mask, xyz, atom_ids):
    selected = np.flatnonzero(mask)
    if not len(selected):
        return {"candidate_atoms": 0, "component_count": 0, "largest_component_atoms": 0,
                "largest_component_dimensions_A": None, "largest_component_convex_hull_volume_A3": None,
                "largest_component_ids": []}
    points = xyz[selected].copy()
    points[:, 0] %= BOX_XY
    points[:, 1] %= BOX_XY
    zlo, zhi = float(points[:, 2].min()), float(points[:, 2].max())
    z_period = max(zhi - zlo + 2 * CONNECTIVITY_A + 1.0e-6, 4 * CONNECTIVITY_A)
    points[:, 2] = points[:, 2] - zlo + CONNECTIVITY_A
    pairs = cKDTree(points, boxsize=[BOX_XY, BOX_XY, z_period]).query_pairs(CONNECTIVITY_A, output_type="ndarray")
    if len(pairs):
        graph = coo_matrix((np.ones(2 * len(pairs)),
                            (np.r_[pairs[:, 0], pairs[:, 1]], np.r_[pairs[:, 1], pairs[:, 0]])),
                           shape=(len(selected), len(selected))).tocsr()
        _, labels = connected_components(graph, directed=False, return_labels=True)
    else:
        labels = np.arange(len(selected))
    counts = np.bincount(labels)
    largest_label = int(counts.argmax())
    ids = selected[labels == largest_label]
    component = points[labels == largest_label].copy()
    anchor = component[0].copy()
    component[:, 0] = anchor[0] + (component[:, 0] - anchor[0] + BOX_XY / 2) % BOX_XY - BOX_XY / 2
    component[:, 1] = anchor[1] + (component[:, 1] - anchor[1] + BOX_XY / 2) % BOX_XY - BOX_XY / 2
    hull_volume = None
    try:
        if len(component) >= 4:
            hull_volume = float(ConvexHull(component).volume)
    except QhullError:
        pass
    return {"candidate_atoms": int(len(selected)), "component_count": int(len(counts)),
            "largest_component_atoms": int(counts[largest_label]),
            "largest_component_dimensions_A": np.ptp(component, axis=0).tolist(),
            "largest_component_convex_hull_volume_A3": hull_volume,
            "largest_component_ids": np.asarray(atom_ids)[ids].astype(int).tolist(),
            "connectivity_cutoff_A": CONNECTIVITY_A, "periodic_xy": True}


def main():
    thresholds = json.loads(THRESHOLD_FILE.read_text())
    if thresholds["status"] != "SEPARATED_FOR_EXPLORATORY_STRUCTURAL_CLASSIFIER":
        raise RuntimeError("Reference thresholds did not separate; refusing composite classification")
    ids0, types0, xyz0 = initial_atoms()
    newton = (xyz0[:, 2] >= 40.0) & (xyz0[:, 2] < 150.0)
    rxy = np.linalg.norm((xyz0[:, :2] - np.array([100.0, 100.0]) + BOX_XY / 2) % BOX_XY - BOX_XY / 2, axis=1)
    depth = np.maximum(0.0, 150.0 - xyz0[:, 2])
    core = ((rxy / 30.0) ** 2 + depth / 20.0 <= 2.0) & newton & (types0 == 4)
    idx_core = np.flatnonzero(core)

    frames, q6_all, ptm_all = [], [], []
    for t in TIMES_PS:
        path = dump_path(float(t))
        if not path.is_file():
            raise FileNotFoundError(path)
        frame = parse_one_frame(path)
        if not np.array_equal(frame["ids"], ids0) or not np.array_equal(frame["types"], types0):
            raise RuntimeError(f"Atom identity mismatch in {path}")
        qids, qtypes, q6, ptm = q6_ptm_at(path)
        if not np.array_equal(qids, ids0) or not np.array_equal(qtypes, types0):
            raise RuntimeError(f"OVITO atom order mismatch in {path}")
        frames.append(frame)
        q6_all.append(q6)
        ptm_all.append(ptm)
        if len(frames) % 5 == 0:
            print(f"structural snapshots: {len(frames)}/{len(TIMES_PS)} at {t:.2f} ps", flush=True)

    q6_all = np.asarray(q6_all)
    ptm_all = np.asarray(ptm_all)
    q6_threshold = float(thresholds["q6_liquid_if_mean_over_1ps_below"])
    ptm_threshold = float(thresholds["ptm_nonfcc_fraction_over_1ps_min"])
    msd_threshold = float(thresholds["per_atom_msd_i1ps_min_A2"])
    windows = []
    latest_core_metrics = None
    for end_index in range(WINDOW_INTERVALS, len(frames)):
        start_index = end_index - WINDOW_INTERVALS
        end_t = frames[end_index]["time_ps"]
        start_t = frames[start_index]["time_ps"]
        if not (abs((end_t - start_t) - 1.0) < 1e-6):
            raise RuntimeError(f"Window length is {end_t - start_t}, expected 1 ps")
        q6_mean = np.mean(q6_all[start_index:end_index + 1, idx_core], axis=0)
        ptm_nonfcc = np.mean(ptm_all[start_index:end_index + 1, idx_core] != FCC_CODE, axis=0)
        displacement = displacement_window(frames[start_index:end_index + 1], newton)
        msd_i = displacement[idx_core] ** 2
        joint = (q6_mean < q6_threshold) & (ptm_nonfcc >= ptm_threshold) & (msd_i >= msd_threshold)
        endpoint_xyz = frames[end_index]["xyz"][idx_core]
        clusters = cluster_stats(joint, endpoint_xyz, ids0[idx_core])
        windows.append({
            "window_ps": [start_t, end_t],
            "joint_liquid_candidate_atoms": int(joint.sum()),
            "q6_joint_candidates_quantiles": np.percentile(q6_mean[joint], [5, 50, 95]).tolist() if np.any(joint) else None,
            "ptm_nonfcc_fraction_joint_candidates_quantiles": np.percentile(ptm_nonfcc[joint], [5, 50, 95]).tolist() if np.any(joint) else None,
            "per_atom_msd_i_A2_joint_candidates_quantiles": np.percentile(msd_i[joint], [5, 50, 95]).tolist() if np.any(joint) else None,
            "source_core_ag_group_msd_i_A2": float(np.mean(msd_i)),
            "joint_candidate_group_msd_i_A2": float(np.mean(msd_i[joint])) if np.any(joint) else None,
            "cluster": clusters,
        })
        if end_index == len(frames) - 1:
            latest_core_metrics = (ids0[idx_core], q6_mean, ptm_nonfcc, msd_i, joint)

    report = {
        "status": "EXPLORATORY_ROLLING_Q6_PTM_MSD_CLASSIFICATION",
        "potential": "Round-5 force-only development MTP; energies/virials uncalibrated",
        "same_atom_set": "initial Newton-layer Ag atoms in the central TTM source core",
        "source_core_ag_atoms": int(core.sum()),
        "structural_sample_interval_ps": 0.05,
        "msd_window_ps": 1.0,
        "q6_definition": "Ag-only Steinhardt q6, cutoff 4.0 A, per-atom time average within each rolling window",
        "ptm_definition": "FCC code 1; non-FCC fraction per Ag over the same 21 snapshots",
        "msd_definition": "unwrapped XY-periodic per-atom squared net displacement; Newton-layer drift removed at each 0.05 ps increment",
        "thresholds": thresholds,
        "windows": windows,
        "caveats": [
            "Classifier is calibrated only against the same development MTP pure-Ag references and remains exploratory.",
            "Large forces and uncalibrated MTP energy/virial prevent quantitative physical melt-pool claims.",
            "Clusters use 3.5 A Ag-Ag connectivity, periodic x/y, open z; convex hull is a geometric envelope.",
        ],
    }
    (HERE / "joint_melt_rolling_report_1p55ps.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    fields = ["window_start_ps", "window_end_ps", "joint_candidate_atoms", "component_count", "largest_component_atoms",
              "largest_component_x_A", "largest_component_y_A", "largest_component_z_A", "source_core_group_msd_i_A2"]
    with (HERE / "joint_melt_rolling_timeseries_to_1p55ps.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for item in windows:
            dims = item["cluster"].get("largest_component_dimensions_A") or [None, None, None]
            writer.writerow({"window_start_ps": item["window_ps"][0], "window_end_ps": item["window_ps"][1],
                             "joint_candidate_atoms": item["joint_liquid_candidate_atoms"],
                             "component_count": item["cluster"]["component_count"],
                             "largest_component_atoms": item["cluster"]["largest_component_atoms"],
                             "largest_component_x_A": dims[0], "largest_component_y_A": dims[1],
                             "largest_component_z_A": dims[2], "source_core_group_msd_i_A2": item["source_core_ag_group_msd_i_A2"]})
    if latest_core_metrics is not None:
        ids, q6_mean, ptm_nonfcc, msd_i, joint = latest_core_metrics
        with (HERE / "joint_classifier_source_core_ag_1p55ps.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(["atom_id", "mean_q6_1ps", "ptm_nonfcc_fraction_1ps", "msd_i_1ps_A2", "joint_candidate"])
            writer.writerows(zip(ids.astype(int), q6_mean, ptm_nonfcc, msd_i, joint.astype(int)))

    times = [x["window_ps"][1] for x in windows]
    liquid = [x["joint_liquid_candidate_atoms"] for x in windows]
    largest = [x["cluster"]["largest_component_atoms"] for x in windows]
    components = [x["cluster"]["component_count"] for x in windows]
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True, facecolor="white")
    axes[0].plot(times, liquid, "o-", label="joint liquid-like Ag candidates", color="#E67E22")
    axes[0].plot(times, largest, "s-", label="largest connected cluster", color="#9C2D20")
    axes[0].set_ylabel("Ag atoms")
    axes[0].legend(frameon=False)
    axes[1].plot(times, components, "o-", color="#305D82")
    axes[1].set_ylabel("connected components")
    axes[1].set_xlabel("1 ps window endpoint (ps)")
    axes[1].set_xticks(times)
    axes[1].grid(alpha=0.25)
    fig.suptitle("Rolling 1 ps Q6/PTM/MSD Ag cluster analysis through 1.55 ps", fontsize=14)
    fig.tight_layout()
    fig.savefig(HERE / "joint_melt_rolling_timeseries_to_1p55ps.png", dpi=180)
    plt.close(fig)
    print(json.dumps({"status": report["status"], "windows": len(windows),
                      "latest": windows[-1], "timeseries_csv": "joint_melt_rolling_timeseries_to_1p55ps.csv"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
