#!/usr/bin/env python3
"""在固定尺寸内同步加厚 Ti3SiC2 骨架和 Ag 网络，保持两相质量比 1:1。"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.data import atomic_masses, atomic_numbers
from ase.io import read, write

import build_step01 as common


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
BASE_MODEL = common.BASE_MODEL
BASE_PHASE = common.BASE_PHASE
BASE_REPORT = common.BASE_REPORT
BUILDER_PATH = common.BUILDER_PATH
OUTPUT = HERE / "output/step05_both_phases"
VISUALIZATION = HERE / "visualization"

# 0.80 Å/侧是在 2.2 Å 相间最小距离和 3.15 Å 银相连通阈值下，
# 从原银网络附近仍可达到目标 1:1 质量比的较大档位。
TI_OFFSET_A = 0.8
AG_GROWTH_RADIUS_A = 8.5
MIN_PHASE_SEPARATION_A = 2.2
AG_CONNECTIVITY_CUTOFF_A = 3.15
STEP = 5
MASS = {s: float(atomic_masses[atomic_numbers[s]]) for s in ("Ag", "Ti", "Si", "C")}
FORMULA_MASS = 3 * MASS["Ti"] + MASS["Si"] + 2 * MASS["C"]


def assert_one_component(builder, positions: np.ndarray, label: str, cutoff: float) -> dict:
    """检查相内最大连通分量；不允许把孤岛悄悄计入模型。"""
    keep, info = builder.component_filter(positions, cutoff=cutoff)
    if not keep.all():
        raise RuntimeError(f"{label} 相不是单一连通网络：{info}")
    return info


def select_connected_ag(builder, candidates: np.ndarray, baseline_ag: np.ndarray,
                        target_atoms: int, ti_positions: np.ndarray) -> tuple[np.ndarray, dict]:
    """在原银网络周围选出目标数量、仍连通且不与加粗骨架重叠的 Ag 位点。"""
    if not len(candidates):
        raise RuntimeError("加粗骨架后没有可用的 Ag 扩展位点")
    keep, component_info = builder.component_filter(candidates, cutoff=AG_CONNECTIVITY_CUTOFF_A)
    candidates = candidates[keep]
    if len(candidates) < target_atoms:
        raise RuntimeError(
            f"固定尺寸和最小相间距下可用连通 Ag 位点不足：{len(candidates)} < {target_atoms}"
        )

    # 优先保留原模型中仍可容纳且属于主网络的银原子，再从扩展位点外缘裁剪。
    baseline_distance = builder.tree(baseline_ag).query(candidates, workers=2)[0]
    is_baseline = baseline_distance < 1.0e-5
    to_remove = len(candidates) - target_atoms
    selected = np.ones(len(candidates), dtype=bool)
    if to_remove:
        ti_distance = builder.tree(ti_positions).query(candidates, workers=2)[0]
        # 新增的远端银位点优先移除；仍需裁剪时，再移除最贴近陶瓷界面的位点。
        priority = np.lexsort((ti_distance, is_baseline.astype(np.int8), -baseline_distance))
        selected[priority[:to_remove]] = False
    ag_positions = candidates[selected]

    connected, trim_info = builder.component_filter(ag_positions, cutoff=AG_CONNECTIVITY_CUTOFF_A)
    if not connected.all():
        # 若裁剪切断枝端，从原候选主网络沿 Ag 最近邻前沿补足原子数。
        ag_positions = ag_positions[connected]
        occupied = builder.tree(ag_positions).query(candidates, workers=2)[0] < 1.0e-7
        ti_distance = builder.tree(ti_positions).query(candidates, workers=2)[0]
        while int(occupied.sum()) < target_atoms:
            available = np.flatnonzero(~occupied)
            frontier_distance = builder.tree(candidates[occupied]).query(
                candidates[available], workers=2
            )[0]
            frontier = available[frontier_distance <= AG_CONNECTIVITY_CUTOFF_A]
            if not len(frontier):
                raise RuntimeError("银网络裁剪后无法沿连通前沿补回 1:1 质量比所需位点")
            # 先取原 Ag 附近、再远离 Ti 界面的候选位点。
            base_dist = builder.tree(baseline_ag).query(candidates[frontier], workers=2)[0]
            order = np.lexsort((ti_distance[frontier], base_dist))
            occupied[frontier[order[:target_atoms - int(occupied.sum())]]] = True
        ag_positions = candidates[occupied]
        connected, trim_info = builder.component_filter(ag_positions, cutoff=AG_CONNECTIVITY_CUTOFF_A)
        if not connected.all():
            raise RuntimeError(f"Ag 相补足目标原子数后仍不连通：{trim_info}")

    if len(ag_positions) != target_atoms:
        raise RuntimeError(f"Ag 原子数不等于 1:1 质量比目标：{len(ag_positions)} != {target_atoms}")
    info = {
        "candidate_main_component": component_info,
        "baseline_ag_atoms_retained_in_candidate_component": int(np.count_nonzero(is_baseline)),
        "trimmed_to_exact_mass_ratio": int(to_remove),
        "final_connectivity": trim_info,
        "candidate_shell_radius_A": AG_GROWTH_RADIUS_A,
        "growth_rule": (f"从原 Ag 晶格网络 {AG_GROWTH_RADIUS_A:g} Å 邻近候选位点扩展；"
                        f"排除距加粗 Ti3SiC2 原子小于 {MIN_PHASE_SEPARATION_A:g} Å 的位置"),
    }
    return ag_positions, info


def main() -> None:
    for path in (BASE_MODEL, BASE_PHASE, BASE_REPORT, BUILDER_PATH):
        if not path.is_file():
            raise FileNotFoundError(f"缺少基准输入：{path}")

    baseline = read(BASE_MODEL)
    phase = np.asarray(baseline.arrays["phase"], dtype=np.int32)
    baseline_ti = baseline.positions[phase == 1]
    baseline_ag = baseline.positions[phase == 2]
    if len(baseline) != 437235 or len(baseline_ti) != 335718 or len(baseline_ag) != 101517:
        raise ValueError("纳米基准模型原子数与已核验的基准值不一致")
    if not np.allclose(baseline.cell.lengths(), common.BOX, atol=1e-6):
        raise ValueError("基准模型模拟盒尺寸异常")

    builder, lattice_builder = common.load_builder()
    mask = np.load(BASE_PHASE)["skeleton"].astype(bool)
    if mask.shape != tuple(common.GRID):
        raise ValueError(f"QSGS 网格形状异常：{mask.shape}")
    field = common.signed_distance(mask)
    ti_pool, ti_symbols_pool, ag_pool, _, lattice_report = lattice_builder.build_lattice_pools(builder)

    # Ti3SiC2 按 signed distance 向原 Ag 侧拓宽，并重新保持 3:1:2 化学计量。
    inside = builder.at_mask(mask, ti_pool)
    distance_field = common.sample_field(field, ti_pool)
    ti_selected = inside | ((~inside) & (distance_field >= -TI_OFFSET_A))
    ti_positions = ti_pool[ti_selected]
    ti_symbols = ti_symbols_pool[ti_selected]
    ti_positions, ti_symbols, units, ti_cleanup = common.connected_phase(
        builder, ti_positions, ti_symbols, mask, ti_phase=True
    )
    ti_component = assert_one_component(builder, ti_positions, "Ti3SiC2", cutoff=3.25)

    # Ag 也由原 Ag 晶格网络向外扩展；允许它进入原 Ti 相掩膜中的原子级空位，
    # 但每个 Ag 位点都需避开实际 Ti/Si/C 原子并属于连通银网。
    baseline_distance = builder.tree(baseline_ag).query(ag_pool, workers=2)[0]
    ti_distance = builder.tree(ti_positions).query(ag_pool, workers=2)[0]
    ag_candidates = ag_pool[
        (baseline_distance <= AG_GROWTH_RADIUS_A)
        & (ti_distance >= MIN_PHASE_SEPARATION_A)
    ]
    target_ag_atoms = round(units * FORMULA_MASS / MASS["Ag"])
    ag_positions, ag_info = select_connected_ag(
        builder, ag_candidates, baseline_ag, target_ag_atoms, ti_positions
    )

    positions = np.vstack((ti_positions, ag_positions))
    symbols = ti_symbols.tolist() + ["Ag"] * len(ag_positions)
    atoms = Atoms(symbols=symbols, positions=positions, cell=common.BOX, pbc=[True, True, False])
    atoms.new_array("phase", np.r_[np.ones(len(ti_positions), np.int32),
                                   np.full(len(ag_positions), 2, np.int32)])
    atoms.new_array("fixed_layer", (positions[:, 2] < 20.0).astype(np.int32))
    atoms.new_array("thermostat_layer", ((positions[:, 2] >= 20.0) & (positions[:, 2] < 40.0)).astype(np.int32))
    atoms.new_array("newton_layer", ((positions[:, 2] >= 40.0) & (positions[:, 2] < 200.0)).astype(np.int32))
    atoms.new_array("sputtering_layer", ((positions[:, 2] >= 200.0) & (positions[:, 2] < 250.0)).astype(np.int32))
    atoms.new_array("evaporation_layer", ((positions[:, 2] >= 250.0) & (positions[:, 2] < 300.0)).astype(np.int32))
    atoms.info.update({
        "status": "UNRELAXED_GEOMETRY_CANDIDATE_STEP05_BOTH_PHASES_THICKENED",
        "thickening_step": STEP,
        "Ti3SiC2_offset_each_side_A": TI_OFFSET_A,
        "Ag_network_growth_radius_A": AG_GROWTH_RADIUS_A,
        "target_Ag_over_Ti3SiC2_mass_ratio": 1.0,
        "material_region_A": common.SIZE.tolist(),
        "phase_order": "phase 1=Ti3SiC2, phase 2=Ag",
    })

    OUTPUT.mkdir(parents=True, exist_ok=True)
    VISUALIZATION.mkdir(parents=True, exist_ok=True)
    extxyz = OUTPUT / "ag_ti3sic2_both_phases_thickened_step05.extxyz"
    datafile = OUTPUT / "ag_ti3sic2_both_phases_thickened_step05.data"
    write(extxyz, atoms, format="extxyz")
    write(datafile, atoms, format="lammps-data", atom_style="atomic",
          specorder=["Ti", "Si", "C", "Ag"], masses=True)
    np.savez_compressed(
        OUTPUT / "both_phases_step05_phase_field.npz",
        original_Ti3SiC2_mask=mask.astype(np.uint8),
        signed_distance_A=field.astype(np.float32),
        Ti3SiC2_offset_each_side_A=np.asarray(TI_OFFSET_A),
        Ag_growth_radius_A=np.asarray(AG_GROWTH_RADIUS_A),
        material_dimensions_A=common.SIZE,
        simulation_box_A=common.BOX,
    )

    # 重新从输出文件读取并核对相标签、盒子、质量比、相间距离与两相连通性。
    check = read(extxyz)
    check_phase = np.asarray(check.arrays["phase"], dtype=np.int32)
    ti_out = check.positions[check_phase == 1]
    ag_out = check.positions[check_phase == 2]
    check_symbols = np.asarray(check.get_chemical_symbols())
    counts = {s: int(np.count_nonzero(check_symbols == s)) for s in ("Ti", "Si", "C", "Ag")}
    check_ti = assert_one_component(builder, ti_out, "Ti3SiC2 输出", cutoff=3.25)
    check_ag = assert_one_component(builder, ag_out, "Ag 输出", cutoff=AG_CONNECTIVITY_CUTOFF_A)
    minimum_separation = float(builder.tree(ti_out).query(ag_out, workers=2)[0].min())
    mass_ratio = counts["Ag"] * MASS["Ag"] / (units * FORMULA_MASS)
    if abs(mass_ratio - 1.0) > 2e-4:
        raise RuntimeError(f"Ag/Ti3SiC2 质量比偏离 1:1：{mass_ratio:.8f}")
    if minimum_separation < MIN_PHASE_SEPARATION_A - 1e-5:
        raise RuntimeError(f"两相原子距离过近：{minimum_separation:.6f} Å")
    if not np.allclose(check.cell.lengths(), common.BOX, atol=1e-6):
        raise RuntimeError("输出模型的模拟盒尺寸发生变化")

    # 比较两相相对基准模型的净增量，确认 Ag 与陶瓷骨架都确实加粗。
    ti_retained = int(np.count_nonzero(builder.tree(ti_out).query(baseline_ti, workers=2)[0] < 1e-5))
    ag_retained = int(np.count_nonzero(builder.tree(ag_out).query(baseline_ag, workers=2)[0] < 1e-5))
    baseline_units = len(baseline_ti) // 6
    ti_added_net = len(ti_out) - len(baseline_ti)
    ag_added_net = len(ag_out) - len(baseline_ag)
    phase_growth_percent = {
        "Ti3SiC2": 100.0 * ti_added_net / len(baseline_ti),
        "Ag": 100.0 * ag_added_net / len(baseline_ag),
    }
    if ti_added_net <= 0 or ag_added_net <= 0:
        raise RuntimeError("目标模型没有同时增加 Ti3SiC2 与 Ag 原子数")

    report = {
        "status": atoms.info["status"],
        "baseline_model": str(BASE_MODEL.relative_to(REPO)),
        "baseline_atoms": len(baseline),
        "baseline_phase_atoms": {"Ti3SiC2": int(len(baseline_ti)), "Ag": int(len(baseline_ag))},
        "baseline_formula_units": int(baseline_units),
        "dimensions": {"material_region_A": common.SIZE.tolist(),
                       "simulation_box_A": common.BOX.tolist(), "pbc": [True, True, False]},
        "thickening": {"Ti3SiC2_offset_each_side_A": TI_OFFSET_A,
                       "Ag_growth_radius_A": AG_GROWTH_RADIUS_A,
                       "Ag_minimum_clearance_from_Ti3SiC2_A": MIN_PHASE_SEPARATION_A,
                       "Ag_connectivity_cutoff_A": AG_CONNECTIVITY_CUTOFF_A,
                       "source_QSGS_grid": common.GRID.tolist(),
                       "method": "Ti3SiC2 expands by signed-distance offset; Ag expands on original Ag lattice sites within 5.8 Å of its network, including atomic-scale voids inside the Ti mask"},
        "lattice": lattice_report,
        "final_atoms": len(check),
        "final_element_counts": counts,
        "final_Ti3SiC2_atoms": int(len(ti_out)),
        "final_Ag_atoms": int(len(ag_out)),
        "Ti3SiC2_formula_units": int(units),
        "Ag_over_Ti3SiC2_mass_ratio": float(mass_ratio),
        "target_Ag_over_Ti3SiC2_mass_ratio": 1.0,
        "phase_net_atom_growth_vs_baseline": {"Ti3SiC2": ti_added_net, "Ag": ag_added_net},
        "phase_atom_growth_percent_vs_baseline": phase_growth_percent,
        "baseline_atoms_retained": {"Ti3SiC2": ti_retained, "Ag": ag_retained},
        "minimum_cross_phase_distance_A": minimum_separation,
        "connectivity": {"Ti3SiC2": check_ti, "Ag": check_ag,
                         "preselection_Ti3SiC2": ti_component, "Ag_selection": ag_info},
        "outputs": {"extxyz": extxyz.name, "lammps_data": datafile.name,
                    "phase_field": "both_phases_step05_phase_field.npz"},
        "limitations": [
            "固定材料区 200 × 200 × 200 Å 和模拟盒 200 × 200 × 300 Å；此模型仍是纳米尺度。",
            "两相在几何上同时扩展且保持约 1:1 质量比；为留出 Ag 扩展空间，Ag 可进入原 Ti 相掩膜内的原子级空位。",
            "这是未弛豫的几何候选，尚未完成势函数力学稳定性、冷态准入或 TTM/MD 加热验证。",
        ],
    }
    report_path = OUTPUT / "both_phases_step05_report.json"
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"status": report["status"], "atoms": len(check),
                      "element_counts": counts, "mass_ratio": mass_ratio,
                      "dimensions_A": report["dimensions"],
                      "net_growth": report["phase_net_atom_growth_vs_baseline"],
                      "growth_percent": phase_growth_percent,
                      "minimum_cross_phase_distance_A": minimum_separation,
                      "connectivity": report["connectivity"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
