#!/usr/bin/env python3
"""Round7 interfacial contact, partial RDF, and MAX-skeleton analysis."""
from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/round7-interface-mpl")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter1d
from ovito.io import import_file


HERE = Path(__file__).resolve().parent
RUNS = Path("/workspace/round7-source/models/ag_ti3sic2_qsgs_200x200x150_md_finer_skeleton/runs")
TRAJ_DIR = RUNS / "round7_fiveframe_joint_analysis_to_1p6ps_20261109"
CORRECTED_DIR = RUNS / "round7_corrected_radial_core_analysis_1p6ps_20261109"
TRAJ = TRAJ_DIR / "round7_5frames_full_system_0p6_to_1p6ps.dump.gz"
CORE_CSV = CORRECTED_DIR / "round7_1p6ps_source_core_ag.csv"
TIMES = np.array([0.60, 0.85, 1.10, 1.35, 1.60])
N_TYPES = {1: "Ti", 2: "Si", 3: "C", 4: "Ag"}
XY_BOX_A = 200.0
MATERIAL_Z_A = 150.0
CORE_CONTACT_A = 3.5
TI_C_CUTOFF_A = 3.0
TI_SI_CUTOFF_A = 3.85
SCAFFOLD_CONNECTIVITY_A = 3.25
RDF_MAX_A = 6.0
RDF_DR_A = 0.05

INK = "#1f2d38"
MUTED = "#61717b"
PALETTE = {"Ti": "#386cb0", "Si": "#f28e2b", "C": "#4e9f71", "Ag": "#9aa6b2"}
FIG_DIRS = (1.1, 0.3, 0.14)
FONT = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"


def read_core_ids(path: Path) -> np.ndarray:
    with path.open(newline="", encoding="utf-8") as stream:
        ids = np.array([int(row["atom_id"]) for row in csv.DictReader(stream)], dtype=np.int64)
    if len(ids) != 1787 or len(np.unique(ids)) != len(ids):
        raise ValueError("Expected the 1,787 unique corrected Round7 source-core Ag IDs")
    return np.sort(ids) - 1


def periodic_tree(positions: np.ndarray) -> cKDTree:
    """Minimum-image x/y only; z is open."""
    points = np.asarray(positions, dtype=float).copy()
    points[:, :2] %= XY_BOX_A
    return cKDTree(points, boxsize=(XY_BOX_A, XY_BOX_A, 0.0))


def frame_data() -> tuple[np.ndarray, np.ndarray, list[int]]:
    pipeline = import_file(str(TRAJ), sort_particles=True)
    if pipeline.source.num_frames != len(TIMES):
        raise ValueError(f"Expected five saved trajectory frames, got {pipeline.source.num_frames}")
    positions, element_types = [], None
    for frame in range(len(TIMES)):
        data = pipeline.compute(frame)
        ids = np.asarray(data.particles["Particle Identifier"], dtype=np.int64)
        current_types = np.asarray(data.particles["Particle Type"], dtype=np.int32)
        if not np.array_equal(ids, np.arange(1, len(ids) + 1)):
            raise ValueError(f"Particle IDs are not sorted/contiguous at frame {frame}")
        if element_types is None:
            element_types = current_types.copy()
        elif not np.array_equal(element_types, current_types):
            raise ValueError(f"Element/type map changed at frame {frame}")
        positions.append(np.asarray(data.particles["Position"], dtype=float).copy())
    return np.stack(positions), element_types, [352637] * len(TIMES)


def pair_matrix(centers: np.ndarray, targets: np.ndarray, cutoff: float) -> np.ndarray:
    return periodic_tree(centers).sparse_distance_matrix(
        periodic_tree(targets), cutoff, output_type="ndarray")


def collect_interface_metrics(frames: np.ndarray, types: np.ndarray,
                             core_indices: np.ndarray) -> tuple[list[dict], np.ndarray, np.ndarray]:
    scaffold_indices = np.flatnonzero(types <= 3)
    phase_indices = {t: np.flatnonzero(types == t) for t in (1, 2, 3, 4)}
    records: list[dict] = []
    initial_core_interface = None
    skeleton_interface = None
    skeleton_displacements: dict[tuple[float, int, str], float] = {}

    increments = np.diff(frames, axis=0)
    increments[:, :, :2] -= XY_BOX_A * np.rint(increments[:, :, :2] / XY_BOX_A)
    mobile_scaffold = ((types <= 3) & (frames[0, :, 2] >= 40.0)
                       & (frames[0, :, 2] < MATERIAL_Z_A))
    drift = np.mean(increments[:, mobile_scaffold, :], axis=1)
    cumulative = np.zeros((len(types), 3), dtype=float)

    for fi, time in enumerate(TIMES):
        pos = frames[fi]
        core_pos = pos[core_indices]
        scaffold_pos = pos[scaffold_indices]
        d_nearest, _ = periodic_tree(scaffold_pos).query(core_pos, workers=2)
        core_interface = d_nearest <= CORE_CONTACT_A
        if fi == 0:
            initial_core_interface = core_interface.copy()

        by_species: dict[str, np.ndarray] = {}
        for t in (1, 2, 3):
            target_pos = pos[phase_indices[t]]
            by_species[N_TYPES[t]] = periodic_tree(target_pos).query_ball_point(
                core_pos, CORE_CONTACT_A, return_length=True, workers=2)

        row = {
            "time_ps": float(time),
            "source_core_ag_atoms": int(len(core_indices)),
            "source_core_ag_nearest_scaffold_median_A": float(np.median(d_nearest)),
            "source_core_ag_nearest_scaffold_p10_A": float(np.quantile(d_nearest, 0.10)),
            "source_core_ag_nearest_scaffold_p90_A": float(np.quantile(d_nearest, 0.90)),
            "source_core_ag_within_3p5A_of_scaffold": int(core_interface.sum()),
            "source_core_ag_within_3p5A_fraction": float(core_interface.mean()),
            "baseline_interface_ag_contact_fraction": float(core_interface[initial_core_interface].mean()),
            "baseline_interface_ag_atoms": int(initial_core_interface.sum()),
        }
        for symbol, values in by_species.items():
            row[f"mean_{symbol}_neighbors_within_3p5A_per_core_Ag"] = float(values.mean())
            row[f"total_{symbol}_contacts_within_3p5A"] = int(values.sum())
        records.append(row)

        if fi > 0:
            cumulative += increments[fi - 1] - drift[fi - 1]
        if fi == 0:
            skeleton_interface = np.zeros(len(types), dtype=bool)
            for t in (1, 2, 3):
                idx = phase_indices[t]
                nearest_ag = periodic_tree(pos[phase_indices[4]]).query(pos[idx], workers=2)[0]
                skeleton_interface[idx] = nearest_ag <= CORE_CONTACT_A
        for t in (1, 2, 3):
            idx = phase_indices[t]
            for group, selected in (("interface", skeleton_interface[idx]),
                                    ("interior", ~skeleton_interface[idx])):
                value = float(np.sqrt(np.mean(np.sum(cumulative[idx[selected]] ** 2, axis=1))))
                skeleton_displacements[(float(time), t, group)] = value

    records_by_time = {row["time_ps"]: row for row in records}
    for time in TIMES:
        row = records_by_time[float(time)]
        for t in (1, 2, 3):
            for group in ("interface", "interior"):
                row[f"{N_TYPES[t]}_{group}_skeleton_rms_displacement_A"] = skeleton_displacements[(float(time), t, group)]
    return records, initial_core_interface, skeleton_interface


def accessible_ball_volume(radius: float, center_z: np.ndarray,
                           zmin: float = 0.0, zmax: float = MATERIAL_Z_A) -> np.ndarray:
    """Volume of a radius-r sphere clipped by the open-z material slab."""
    lo = np.maximum(-radius, zmin - center_z)
    hi = np.minimum(radius, zmax - center_z)
    valid = hi > lo
    out = np.zeros_like(center_z, dtype=float)
    primitive = lambda z: np.pi * (radius * radius * z - z ** 3 / 3.0)
    out[valid] = primitive(hi[valid]) - primitive(lo[valid])
    return out


def collect_partial_rdfs(frames: np.ndarray, types: np.ndarray,
                         core_indices: np.ndarray) -> list[dict]:
    edges = np.linspace(0.0, RDF_MAX_A, int(round(RDF_MAX_A / RDF_DR_A)) + 1)
    centers_r = 0.5 * (edges[:-1] + edges[1:])
    material_volume = XY_BOX_A * XY_BOX_A * MATERIAL_Z_A
    records: list[dict] = []
    pair_targets = [(4, "Ag-Ag"), (1, "Ag-Ti"), (2, "Ag-Si"), (3, "Ag-C")]
    target_indices = {t: np.flatnonzero((types == t)) for t, _ in pair_targets}
    for fi, time in enumerate(TIMES):
        pos = frames[fi]
        core_pos = pos[core_indices]
        z_core = core_pos[:, 2]
        if np.any(z_core < 0) or np.any(z_core > MATERIAL_Z_A + 2.0):
            raise ValueError("A tracked source-core Ag moved far outside the material slab")
        shell_volumes = np.stack([
            accessible_ball_volume(float(hi), z_core) - accessible_ball_volume(float(lo), z_core)
            for lo, hi in zip(edges[:-1], edges[1:])
        ], axis=1).sum(axis=0)

        for target_type, pair_name in pair_targets:
            target_indices_in_material = target_indices[target_type][
                (pos[target_indices[target_type], 2] >= 0.0)
                & (pos[target_indices[target_type], 2] < MATERIAL_Z_A)]
            target_pos = pos[target_indices_in_material]
            matrix = pair_matrix(core_pos, target_pos, RDF_MAX_A)
            distances = matrix["v"] if len(matrix) else np.empty(0)
            if target_type == 4 and len(matrix):
                same_atom = core_indices[matrix["i"]] == target_indices_in_material[matrix["j"]]
                distances = distances[~same_atom]
            counts, _ = np.histogram(distances, bins=edges)
            rho = len(target_indices_in_material) / material_volume
            expected = rho * shell_volumes
            g = np.divide(counts, expected, out=np.zeros_like(expected), where=expected > 0)
            for r, value, n in zip(centers_r, g, counts):
                records.append({"time_ps": float(time), "pair": pair_name,
                                "r_A": float(r), "conditional_partial_g": float(value),
                                "neighbor_pairs_in_bin": int(n),
                                "core_centers": int(len(core_indices)),
                                "target_atoms_in_material_slab": int(len(target_indices_in_material))})
    return records


def pair_shell_data(frames: np.ndarray, types: np.ndarray,
                    skeleton_interface: np.ndarray) -> tuple[list[dict], dict]:
    idx_type = {t: np.flatnonzero(types == t) for t in (1, 2, 3)}
    definitions = [(1, 3, "Ti-C", TI_C_CUTOFF_A), (1, 2, "Ti-Si", TI_SI_CUTOFF_A)]
    baseline_pairs: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for ta, tb, name, cutoff in definitions:
        mat = pair_matrix(frames[0, idx_type[ta]], frames[0, idx_type[tb]], cutoff)
        first = idx_type[ta][mat["i"]]
        second = idx_type[tb][mat["j"]]
        is_interface = skeleton_interface[first] | skeleton_interface[second]
        baseline_pairs[name] = (first, second, is_interface)

    records: list[dict] = []
    connectivity = []
    mobile_ag = np.flatnonzero(types == 4)
    for fi, time in enumerate(TIMES):
        pos = frames[fi]
        for ta, tb, name, cutoff in definitions:
            ia, ib = idx_type[ta], idx_type[tb]
            matrix = pair_matrix(pos[ia], pos[ib], cutoff)
            counts_ti = np.bincount(matrix["i"], minlength=len(ia))
            interface_mask = skeleton_interface[ia]
            pair_first, pair_second, pair_interface = baseline_pairs[name]
            delta = pos[pair_first] - pos[pair_second]
            delta[:, :2] -= XY_BOX_A * np.rint(delta[:, :2] / XY_BOX_A)
            lengths = np.linalg.norm(delta, axis=1)
            retained = lengths <= cutoff
            for group, atom_mask, pair_mask in (
                ("interface", interface_mask, pair_interface),
                ("interior", ~interface_mask, ~pair_interface),
            ):
                group_counts = counts_ti[atom_mask]
                group_lengths = lengths[pair_mask]
                records.append({
                    "time_ps": float(time), "pair": name, "group": group,
                    "first_shell_cutoff_A": float(cutoff),
                    "mean_central_Ti_coordination": float(group_counts.mean()),
                    "median_baseline_pair_length_A": float(np.median(group_lengths)),
                    "q25_baseline_pair_length_A": float(np.quantile(group_lengths, 0.25)),
                    "q75_baseline_pair_length_A": float(np.quantile(group_lengths, 0.75)),
                    "baseline_first_shell_pairs": int(pair_mask.sum()),
                    "baseline_pairs_retained_within_cutoff": int(np.sum(retained[pair_mask])),
                    "baseline_pair_retention_fraction": float(retained[pair_mask].mean()),
                    "central_Ti_atoms": int(atom_mask.sum()),
                })

        scaffold_idx = np.flatnonzero(types <= 3)
        scaffold_tree = periodic_tree(pos[scaffold_idx])
        edge_pairs = scaffold_tree.query_pairs(SCAFFOLD_CONNECTIVITY_A, output_type="ndarray")
        graph = coo_matrix((np.ones(len(edge_pairs), dtype=np.uint8),
                            (edge_pairs[:, 0], edge_pairs[:, 1])),
                           shape=(len(scaffold_idx), len(scaffold_idx))).tocsr()
        ncomp, labels = connected_components(graph, directed=False)
        sizes = np.bincount(labels)
        connectivity.append({"time_ps": float(time), "scaffold_atoms": int(len(scaffold_idx)),
                             "same_phase_distance_cutoff_A": SCAFFOLD_CONNECTIVITY_A,
                             "connected_components": int(ncomp),
                             "largest_component_atoms": int(sizes.max()),
                             "largest_component_fraction": float(sizes.max() / len(scaffold_idx)),
                             "geometric_edges": int(len(edge_pairs))})
        del graph, edge_pairs, labels, sizes, scaffold_tree
    return records, {"frames": connectivity, "baseline_pairs_by_type": {
        name: int(len(values[0])) for name, values in baseline_pairs.items()}}


def save_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path}")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def figure_style() -> None:
    if Path(FONT).exists():
        from matplotlib import font_manager
        font_manager.fontManager.addfont(FONT)
        plt.rcParams["font.family"] = font_manager.FontProperties(fname=FONT).get_name()
    plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False,
                         "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": MUTED, "text.color": INK,
                         "axes.titleweight": "bold", "figure.facecolor": "white",
                         "axes.facecolor": "white", "savefig.facecolor": "white"})


def save_figure(fig: plt.Figure, stem: str) -> None:
    fig.savefig(HERE / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(HERE / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(HERE / f"{stem}.tif", dpi=300, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def plot_interface(metrics: list[dict]) -> None:
    figure_style()
    fig, axes = plt.subplots(1, 3, figsize=(15.2, 4.9), constrained_layout=True)
    t = np.array([row["time_ps"] for row in metrics])
    median = np.array([row["source_core_ag_nearest_scaffold_median_A"] for row in metrics])
    p10 = np.array([row["source_core_ag_nearest_scaffold_p10_A"] for row in metrics])
    p90 = np.array([row["source_core_ag_nearest_scaffold_p90_A"] for row in metrics])
    ax = axes[0]
    ax.fill_between(t, p10, p90, color="#86aebf", alpha=0.25, label="第10–90百分位")
    ax.plot(t, median, "o-", color="#356f86", lw=2.2, label="中位最近距离")
    ax.axhline(CORE_CONTACT_A, color="#d17c45", ls="--", lw=1.3, label="接触统计阈值 3.5 Å")
    ax.set(title="(a) 热源 Ag 到 MAX 骨架的距离", xlabel="时间 / ps", ylabel="最近 Ti/Si/C 距离 / Å")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    for symbol in ("Ti", "Si", "C"):
        key = f"mean_{symbol}_neighbors_within_3p5A_per_core_Ag"
        ax.plot(t, [row[key] for row in metrics], "o-", lw=2, color=PALETTE[symbol], label=f"Ag–{symbol}")
    ax.set(title="(b) Ag–骨架近邻接触", xlabel="时间 / ps", ylabel="每个核心 Ag 的平均近邻数（≤3.5 Å）")
    ax.legend(frameon=False, fontsize=9)

    ax = axes[2]
    for symbol in ("Ti", "Si", "C"):
        key_i = f"{symbol}_interface_skeleton_rms_displacement_A"
        key_b = f"{symbol}_interior_skeleton_rms_displacement_A"
        ax.plot(t, [row[key_i] for row in metrics], "o-", lw=2, color=PALETTE[symbol], label=f"{symbol} 界面邻域")
        ax.plot(t, [row[key_b] for row in metrics], "o--", lw=1.7, color=PALETTE[symbol], alpha=0.65, label=f"{symbol} 骨架内部")
    ax.set(title="(c) 骨架原子的漂移修正迁移", xlabel="时间 / ps", ylabel="RMS 位移 / Å")
    ax.legend(frameon=False, fontsize=7, ncol=2)
    fig.suptitle("Ag/Ti3SiC2 界面接触与骨架迁移", x=0.04, y=1.06, ha="left", fontsize=17, fontweight="bold")
    fig.text(0.04, -0.02, "固定跟踪 1,787 个热源核心 Ag；骨架界面邻域按 0.60 ps 时到最近 Ag ≤3.5 Å 定义。", fontsize=9, color=MUTED)
    save_figure(fig, "Figure11_Round7_interface_contact_migration")


def plot_rdf(records: list[dict], rdf_center_count: int) -> None:
    figure_style()
    fig, axes = plt.subplots(2, 2, figsize=(11.3, 8.1), constrained_layout=True, sharex=True)
    pairs = ("Ag-Ag", "Ag-Ti", "Ag-Si", "Ag-C")
    cmap = plt.get_cmap("viridis")
    time_colors = [cmap(x) for x in np.linspace(0.08, 0.9, len(TIMES))]
    for ax, pair in zip(axes.flat, pairs):
        for time, color in zip(TIMES, time_colors):
            subset = [row for row in records if row["pair"] == pair and row["time_ps"] == float(time)]
            raw = np.array([row["conditional_partial_g"] for row in subset])
            ax.plot([row["r_A"] for row in subset], gaussian_filter1d(raw, sigma=1.2),
                    color=color, lw=1.7, label=f"{time:.2f} ps")
        ax.set_title(pair, fontsize=12, fontweight="bold")
        ax.grid(axis="y", color="#e7ecef", lw=0.7)
        ax.set_xlim(0, RDF_MAX_A)
        ax.set_ylabel("条件部分 g(r)")
    for ax in axes[-1, :]:
        ax.set_xlabel("原子间距 r / Å")
    axes[0, 0].legend(frameon=False, fontsize=8, ncol=2)
    fig.suptitle("热源核心 Ag 周围的局部部分径向分布", x=0.07, y=1.04, ha="left", fontsize=17, fontweight="bold")
    fig.text(0.07, -0.015,
             f"Ag 中心为同一组留在 z<150 Å 的 {rdf_center_count:,} 个核心 Ag；目标原子限制在材料层。按材料层密度及开放 z 壳层体积归一化。",
             fontsize=8.5, color=MUTED)
    save_figure(fig, "Figure12_Round7_interfacial_partial_RDF")


def plot_max_order(order: list[dict], continuity: dict) -> None:
    figure_style()
    fig, axes = plt.subplots(2, 2, figsize=(11.8, 8.4), constrained_layout=True)
    t = np.array([row["time_ps"] for row in continuity["frames"]])
    pair_colors = {"Ti-C": "#277da1", "Ti-Si": "#f8961e"}
    styles = {"interface": "-", "interior": "--"}

    ax = axes[0, 0]
    for pair in ("Ti-C", "Ti-Si"):
        for group in ("interface", "interior"):
            rows = [r for r in order if r["pair"] == pair and r["group"] == group]
            ax.plot([r["time_ps"] for r in rows], [r["mean_central_Ti_coordination"] for r in rows],
                    marker="o", lw=1.8, color=pair_colors[pair], ls=styles[group], label=f"{pair}｜{group}")
    ax.set(title="(a) Ti 的第一近邻配位", xlabel="时间 / ps", ylabel="平均近邻数 / Ti")
    ax.legend(frameon=False, fontsize=8, ncol=2)

    ax = axes[0, 1]
    for pair in ("Ti-C", "Ti-Si"):
        for group in ("interface", "interior"):
            rows = [r for r in order if r["pair"] == pair and r["group"] == group]
            med = np.array([r["median_baseline_pair_length_A"] for r in rows])
            q25 = np.array([r["q25_baseline_pair_length_A"] for r in rows])
            q75 = np.array([r["q75_baseline_pair_length_A"] for r in rows])
            x = np.array([r["time_ps"] for r in rows])
            ax.plot(x, med, marker="o", lw=1.8, color=pair_colors[pair], ls=styles[group], label=f"{pair}｜{group}")
            ax.fill_between(x, q25, q75, color=pair_colors[pair], alpha=0.08)
    ax.set(title="(b) 0.60 ps 第一近邻对的键长代理", xlabel="时间 / ps", ylabel="中位近邻距离 / Å")
    ax.legend(frameon=False, fontsize=8, ncol=2)

    ax = axes[1, 0]
    for pair in ("Ti-C", "Ti-Si"):
        for group in ("interface", "interior"):
            rows = [r for r in order if r["pair"] == pair and r["group"] == group]
            ax.plot([r["time_ps"] for r in rows],
                    [100 * r["baseline_pair_retention_fraction"] for r in rows],
                    marker="o", lw=1.8, color=pair_colors[pair], ls=styles[group], label=f"{pair}｜{group}")
    ax.set(title="(c) 初始近邻对保留率", xlabel="时间 / ps", ylabel="仍在固定截断内 / %", ylim=(0, 102))
    ax.legend(frameon=False, fontsize=8, ncol=2)

    ax = axes[1, 1]
    largest = np.array([100 * row["largest_component_fraction"] for row in continuity["frames"]])
    components = np.array([row["connected_components"] for row in continuity["frames"]])
    ax.plot(t, largest, "o-", color="#2a7f62", lw=2.2, label="最大 Ti/Si/C 连通分量")
    ax.set(title="(d) MAX 骨架几何连通性", xlabel="时间 / ps", ylabel="最大连通分量 / 全骨架 %", ylim=(0, 102))
    ax2 = ax.twinx()
    ax2.spines["top"].set_visible(False)
    ax2.plot(t, components, "s--", color="#c2544a", lw=1.6, label="连通分量数")
    ax2.set_ylabel("连通分量数", color="#c2544a")
    handles = ax.get_legend_handles_labels()[0] + ax2.get_legend_handles_labels()[0]
    labels = ax.get_legend_handles_labels()[1] + ax2.get_legend_handles_labels()[1]
    ax.legend(handles, labels, frameon=False, fontsize=8, loc="lower left")

    fig.suptitle("Ti3SiC2 MAX 骨架的近邻结构与连续性", x=0.06, y=1.04, ha="left", fontsize=17, fontweight="bold")
    fig.text(0.06, -0.015,
             f"固定首帧阈值：Ti–C ≤{TI_C_CUTOFF_A:.2f} Å、Ti–Si ≤{TI_SI_CUTOFF_A:.2f} Å（首帧径向分布首谷）；连通性为任意 Ti/Si/C 原子间距 ≤{SCAFFOLD_CONNECTIVITY_A:.2f} Å。",
             fontsize=8.3, color=MUTED)
    save_figure(fig, "Figure13_Round7_MAX_skeleton_local_order")


def main() -> None:
    for path in (TRAJ, CORE_CSV):
        if not path.exists():
            raise FileNotFoundError(path)
    frames, types, counts = frame_data()
    if len(types) != 352637 or not np.array_equal(np.unique(types), np.array([1, 2, 3, 4])):
        raise ValueError("Unexpected Round7 atom count or element type labels")
    core_indices = read_core_ids(CORE_CSV)
    if not np.all(types[core_indices] == 4):
        raise ValueError("Corrected source-core IDs include non-Ag atoms")

    interface_metrics, initial_core_interface, skeleton_interface = collect_interface_metrics(frames, types, core_indices)
    # Keep a fixed center population that remains inside the material slab at every saved frame.
    # Core Ag IDs that cross z=150 Å are excluded from this fixed material-slab RDF only.
    rdf_core_indices = core_indices[np.all(frames[:, core_indices, 2] < MATERIAL_Z_A, axis=0)]
    if len(rdf_core_indices) < len(core_indices) - 10:
        raise ValueError("Too many source-core Ag atoms left the material slab for a fixed-center RDF")
    rdf_records = collect_partial_rdfs(frames, types, rdf_core_indices)
    max_order, continuity = pair_shell_data(frames, types, skeleton_interface)

    save_csv(HERE / "Figure11_round7_interface_metrics.csv", interface_metrics)
    save_csv(HERE / "Figure12_round7_interfacial_partial_RDF.csv", rdf_records)
    save_csv(HERE / "Figure13_round7_max_local_order.csv", max_order)
    plot_interface(interface_metrics)
    plot_rdf(rdf_records, len(rdf_core_indices))
    plot_max_order(max_order, continuity)

    mins = {}
    for t, name in ((1, "Ti"), (2, "Si"), (3, "C")):
        distance = periodic_tree(frames[0, types == 4]).query(frames[0, types == t], workers=2)[0]
        mins[name] = float(distance.min())
    report = {
        "status": "exploratory_round7_interface_and_MAX_local_structure_analysis",
        "times_ps": TIMES.tolist(),
        "trajectory": str(TRAJ),
        "trajectory_sha256": hashlib.file_digest(TRAJ.open("rb"), "sha256").hexdigest(),
        "corrected_source_core_csv_sha256": hashlib.file_digest(CORE_CSV.open("rb"), "sha256").hexdigest(),
        "atom_count_per_frame": counts,
        "type_counts": {N_TYPES[int(t)]: int(np.sum(types == t)) for t in np.unique(types)},
        "source_core_ag_count": int(len(core_indices)),
        "baseline_interface_ag_count_nearest_scaffold_within_3p5A_at_0p60ps": int(initial_core_interface.sum()),
        "interface_contact_cutoff_A": CORE_CONTACT_A,
        "rdf_definition": "Ag-core-centered partial g_AB(r), normalized by whole-material-slab element density and exact open-z shell accessibility; the bicontinuous local domain is inhomogeneous, so long-range convergence to unity is not assumed.",
        "rdf_center_population": f"same {len(rdf_core_indices)} corrected source-core Ag IDs that remain at z<150 A in all frames; atoms that leave the material slab are excluded from the RDF only",
        "rdf_center_count": int(len(rdf_core_indices)),
        "rdf_pair_types": ["Ag-Ag", "Ag-Ti", "Ag-Si", "Ag-C"],
        "rdf_material_z_slab_A": [0.0, MATERIAL_Z_A],
        "rdf_bin_width_A": RDF_DR_A,
        "rdf_max_A": RDF_MAX_A,
        "max_first_shell_cutoffs_A": {"Ti-C": TI_C_CUTOFF_A, "Ti-Si": TI_SI_CUTOFF_A},
        "max_shell_definition": "Fixed at first minimum after the first peak in 0.60 ps partial-distance histogram; geometric near-neighbor proxy, not bond-energy decomposition.",
        "max_skeleton_graph_cutoff_A": SCAFFOLD_CONNECTIVITY_A,
        "max_skeleton_connectivity": continuity["frames"],
        "baseline_MAX_pair_counts": continuity["baseline_pairs_by_type"],
        "min_Ag_to_scaffold_distance_by_type_at_0p60ps_A": mins,
        "potential_status": "The Round7 development MTP is not independently validated for liquid Ag or Ag/Ti3SiC2 interface transferability; the run summary reports a maximum force of 133.006 eV/A during the continuation.",
        "limitations": [
            "The QSGS phases form a bicontinuous, curved interface, so there is no unique global interface normal. Figure 11 uses source-core Ag distance and local cross-species contacts instead of an arbitrary planar concentration profile.",
            "The trajectory has only five stored frames 0.25 ps apart; sub-frame diffusion paths and fast bond rearrangements are not resolved.",
            "Geometric pair counts and distances do not measure bond energy or chemical bonding strength.",
            "No per-atom stress or force vectors are present in the trajectory dump; stress concentrations at scaffold branches cannot be inferred from these files.",
            "Connected components use a geometric 3.25 A criterion and do not prove mechanical load-bearing continuity.",
            "Results describe this exploratory trajectory and do not establish validated interfacial diffusion, reaction, melting, or welding."
        ],
        "outputs": [
            "Figure11_Round7_interface_contact_migration.png/.pdf/.tif",
            "Figure12_Round7_interfacial_partial_RDF.png/.pdf/.tif",
            "Figure13_Round7_MAX_skeleton_local_order.png/.pdf/.tif",
            "Figure11_round7_interface_metrics.csv",
            "Figure12_round7_interfacial_partial_RDF.csv",
            "Figure13_round7_max_local_order.csv"
        ]
    }
    (HERE / "round7_interface_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (HERE / "round7_interface_methods.md").write_text(methods_text(report), encoding="utf-8")
    print(json.dumps({"output_dir": str(HERE), "times_ps": TIMES.tolist(),
                      "baseline_interface_ag": report["baseline_interface_ag_count_nearest_scaffold_within_3p5A_at_0p60ps"],
                      "continuity": continuity["frames"]}, ensure_ascii=False, indent=2))


def methods_text(report: dict) -> str:
    return f"""# Round7 Ag/Ti₃SiC₂ 界面与 MAX 骨架分析方法

## 图11：界面接触与骨架迁移

固定跟踪更正源区筛选得到的 {report['source_core_ag_count']:,} 个 Ag ID，按 0.60、0.85、1.10、1.35、1.60 ps 五帧统计 Ag 到最近 Ti/Si/C 原子的距离，以及每个 Ag 在 3.5 Å 内的 Ti、Si、C 邻居数。右图按 0.60 ps 时最近 Ag 距离 ≤3.5 Å 将 Ti、Si、C 骨架原子标为界面邻域，比较该组与骨架内部的漂移修正 RMS 位移。漂移取每一帧间隔、z=40–150 Å 骨架原子的平均位移；x/y 采用周期最小镜像，z 不周期。

本体系为双连续 QSGS 形貌，界面弯曲且互相贯通，没有一个可代表全域的平直界面法向。因此不把某个任意 x/y/z 方向包装成“界面法向浓度”；图11改用局部距离和元素接触统计。

## 图12：热源核心 Ag 的条件部分 RDF

以同一组源区核心 Ag 为中心，统计 Ag–Ag、Ag–Ti、Ag–Si、Ag–C 距离分布；目标原子限制在材料层 z=0–150 Å。壳层按该元素在材料层的平均数密度和开放 z 边界下可达的球壳体积归一化，x/y 周期。由于 QSGS 两相在空间上强烈非均匀，长程基线不强制为 1；主要比较第一近邻峰的位置、宽度和时间变化。它是局部条件 RDF，不可当成均匀体相 RDF。

## 图13：MAX 骨架局部结构和几何连通

不使用 Ag 的 FCC/PTM 分类判据。Ti–C 与 Ti–Si 首壳截断值分别为 {TI_C_CUTOFF_A:.2f} Å 与 {TI_SI_CUTOFF_A:.2f} Å，依据 0.60 ps 距离直方图第一峰后的首个谷值选定，并固定用于所有帧。配位统计比较 Ti 周围 C 与 Si 的平均近邻数。键长代理跟踪 0.60 ps 首壳内同一对原子，并给出中位数及四分位数；保留率表示这些原子对在当前帧仍处于首壳截断内的比例，不是化学键能或反应判据。界面/内部按首帧骨架原子是否距最近 Ag ≤3.5 Å 分类。

骨架几何连续性以 Ti/Si/C 任意两原子间距 ≤{SCAFFOLD_CONNECTIVITY_A:.2f} Å 建图，报告最大连通分量比例和分量数。它是几何连通筛查，不等于力学承载或化学稳定性。

## 数据边界

输入是已有 Round7 全体系五帧 dump，时间间隔 0.25 ps。轨迹文件未保存逐原子力/应力向量，故无法识别分叉点应力集中。轨迹所属 Round7 开发 MTP 尚未独立验证 Ag 熔化与 Ag/Ti₃SiC₂ 界面可迁移性，且运行摘要报告了 133.006 eV/Å 的最大检查点力。因此这些图仅量化该探索性轨迹中的几何接触与局部结构变化，不证明界面扩散、化学反应、熔化或熔焊。
"""


if __name__ == "__main__":
    main()
