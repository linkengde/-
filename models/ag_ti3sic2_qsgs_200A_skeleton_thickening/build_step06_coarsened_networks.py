#!/usr/bin/env python3
"""固定相质量与盒子尺寸，将两相网络粗化为较粗、较平滑的连通枝干。"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.data import atomic_masses, atomic_numbers
from ase.io import read, write
from scipy import ndimage as ndi

import build_step01 as common


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
BASELINE = HERE / "output/step05_both_phases/ag_ti3sic2_both_phases_thickened_step05.extxyz"
OUTPUT = HERE / "output/step06_coarsened_networks"
GRID = np.asarray(common.GRID)
VOXEL = common.VOXEL
SIGMA_VOXELS = 1.5
MIN_CROSS_PHASE_DISTANCE_A = 2.2
TI_CONNECTIVITY_CUTOFF_A = 3.25
AG_CONNECTIVITY_CUTOFF_A = 3.15
TI_CANDIDATE_MULTIPLIER = 1.02
AG_CANDIDATE_MULTIPLIER = 1.15
MASS = {s: float(atomic_masses[atomic_numbers[s]]) for s in ("Ag", "Ti", "Si", "C")}
FORMULA_MASS = 3 * MASS["Ti"] + MASS["Si"] + 2 * MASS["C"]


def occupancy_grid(positions: np.ndarray) -> np.ndarray:
    """把相原子按位置计数到 QSGS 网格。"""
    idx = np.floor(positions / VOXEL).astype(int)
    idx = np.clip(idx, 0, GRID - 1)
    flat = np.ravel_multi_index(idx.T, tuple(GRID))
    return np.bincount(flat, minlength=int(np.prod(GRID))).reshape(tuple(GRID)).astype(np.float32)


def smooth_periodic_xy(field: np.ndarray) -> np.ndarray:
    """在周期 x/y、开放 z 的边界条件下平滑相占位密度。"""
    tiled = np.tile(field, (3, 3, 1))
    smoothed = ndi.gaussian_filter(
        tiled, sigma=(SIGMA_VOXELS, SIGMA_VOXELS, SIGMA_VOXELS), mode="nearest"
    )
    return smoothed[GRID[0]:2 * GRID[0], GRID[1]:2 * GRID[1], :]


def sample_scores(positions: np.ndarray, field: np.ndarray) -> np.ndarray:
    """将平滑相密度采样到候选晶格点。"""
    idx = np.floor(positions / VOXEL).astype(int)
    idx = np.clip(idx, 0, GRID - 1)
    return field[idx[:, 0], idx[:, 1], idx[:, 2]]


def top_indices(ids: np.ndarray, scores: np.ndarray, count: int) -> np.ndarray:
    """按相场密度稳定选取指定数量的候选位点。"""
    if len(ids) < count:
        raise RuntimeError(f"候选晶格位点不足：{len(ids)} < {count}")
    local = np.argpartition(scores[ids], -count)[-count:]
    chosen = ids[local]
    return chosen[np.argsort(scores[chosen], kind="stable")[::-1]]


def connected_exact_selection(builder, pool: np.ndarray, symbols: np.ndarray,
                              scores: np.ndarray, target: dict[str, int],
                              cutoff: float, label: str) -> tuple[np.ndarray, np.ndarray, dict]:
    """从平滑后的高分位点中选出精确配额，并沿连通前沿补全主网。"""
    candidate_ids = []
    for species, quota in target.items():
        species_ids = np.flatnonzero(symbols == species)
        n_candidates = min(len(species_ids), int(np.ceil(quota * TI_CANDIDATE_MULTIPLIER)))
        if species == "Ag":
            n_candidates = min(len(species_ids), int(np.ceil(quota * AG_CANDIDATE_MULTIPLIER)))
        candidate_ids.append(top_indices(species_ids, scores, n_candidates))
    candidate_ids = np.concatenate(candidate_ids)
    candidate_ids = np.unique(candidate_ids)
    candidate_pos = pool[candidate_ids]
    candidate_symbol = symbols[candidate_ids]
    candidate_score = scores[candidate_ids]

    candidate_keep, candidate_info = builder.component_filter(candidate_pos, cutoff=cutoff)
    candidate_ids = candidate_ids[candidate_keep]
    candidate_pos = pool[candidate_ids]
    candidate_symbol = symbols[candidate_ids]
    candidate_score = scores[candidate_ids]
    if not len(candidate_ids):
        raise RuntimeError(f"{label} 平滑候选没有连通分量")

    selected_ids = []
    for species, quota in target.items():
        ids = np.flatnonzero(candidate_symbol == species)
        if len(ids) < quota:
            raise RuntimeError(f"{label} 主连通候选中的 {species} 位点不足：{len(ids)} < {quota}")
        local = np.argpartition(candidate_score[ids], -quota)[-quota:]
        selected_ids.append(candidate_ids[ids[local]])
    selected_ids = np.unique(np.concatenate(selected_ids))

    # 若按平滑密度截取后留下了小孤岛，先保留最大网络，再沿候选网络前沿补足各元素配额。
    selected_pos = pool[selected_ids]
    selected_keep, selected_info = builder.component_filter(selected_pos, cutoff=cutoff)
    if not selected_keep.all():
        selected_ids = selected_ids[selected_keep]
    for _ in range(12):
        current_symbols = symbols[selected_ids]
        deficit = {s: n - int(np.count_nonzero(current_symbols == s)) for s, n in target.items()}
        if all(n == 0 for n in deficit.values()):
            break
        if any(n < 0 for n in deficit.values()):
            raise RuntimeError(f"{label} 连通修整后元素数量超过配额：{deficit}")
        current_tree = builder.tree(pool[selected_ids])
        remaining = np.setdiff1d(candidate_ids, selected_ids, assume_unique=True)
        frontier_dist = current_tree.query(pool[remaining], workers=2)[0]
        frontier_ids = remaining[frontier_dist <= cutoff]
        additions = []
        for species, missing in deficit.items():
            if missing <= 0:
                continue
            ids = frontier_ids[symbols[frontier_ids] == species]
            if len(ids) < missing:
                raise RuntimeError(
                    f"{label} 无法在连通前沿补足 {species} 配额：可补 {len(ids)}，需 {missing}"
                )
            additions.append(top_indices(ids, scores, missing))
        if not additions:
            break
        selected_ids = np.unique(np.concatenate((selected_ids, *additions)))
    else:
        raise RuntimeError(f"{label} 连通前沿修整超过迭代上限")

    final_symbols = symbols[selected_ids]
    final_counts = {s: int(np.count_nonzero(final_symbols == s)) for s in target}
    if final_counts != target:
        raise RuntimeError(f"{label} 元素配额不符：{final_counts} != {target}")
    final_keep, final_info = builder.component_filter(pool[selected_ids], cutoff=cutoff)
    if not final_keep.all():
        raise RuntimeError(f"{label} 最终相仍不连通：{final_info}")
    return pool[selected_ids], final_symbols, {
        "candidate_main_component": candidate_info,
        "initial_selected_connectivity": selected_info,
        "final_connectivity": final_info,
        "candidate_pool_counts": {s: int(np.count_nonzero(candidate_symbol == s)) for s in target},
    }


def equal_volume_smooth_mask(field: np.ndarray, original_occupancy: np.ndarray) -> np.ndarray:
    """生成与基准相同占位体素数的粗化掩膜，用于估计枝干厚度变化。"""
    count = int(np.count_nonzero(original_occupancy))
    flat = field.ravel()
    chosen = np.argpartition(flat, -count)[-count:]
    out = np.zeros(flat.size, dtype=bool)
    out[chosen] = True
    return out.reshape(tuple(GRID))


def thickness_proxy(mask: np.ndarray) -> dict:
    """按网格距离变换计算等体积分辨掩膜的局部半厚度统计。"""
    tiled = np.tile(mask, (3, 3, 1))
    dist = ndi.distance_transform_edt(tiled, sampling=VOXEL)
    local = dist[GRID[0]:2 * GRID[0], GRID[1]:2 * GRID[1], :][mask]
    return {"mean_local_half_thickness_A": float(local.mean()),
            "median_local_half_thickness_A": float(np.median(local)),
            "p90_local_half_thickness_A": float(np.percentile(local, 90)),
            "occupied_voxels": int(mask.sum())}


def main() -> None:
    if not BASELINE.is_file():
        raise FileNotFoundError(f"缺少 step05 基准原子模型：{BASELINE}")
    atoms0 = read(BASELINE)
    phase0 = np.asarray(atoms0.arrays["phase"], dtype=np.int32)
    ti0 = atoms0.positions[phase0 == 1]
    ag0 = atoms0.positions[phase0 == 2]
    base_symbols = np.asarray(atoms0.get_chemical_symbols())
    base_counts = {s: int(np.count_nonzero(base_symbols == s)) for s in ("Ti", "Si", "C", "Ag")}
    if base_counts != {"Ti": 175257, "Si": 58419, "C": 116838, "Ag": 105991}:
        raise ValueError(f"step05 相原子数与质量比基准不符：{base_counts}")
    if not np.allclose(atoms0.cell.lengths(), common.BOX, atol=1e-6):
        raise ValueError("step05 模拟盒尺寸异常")

    builder, lattice_builder = common.load_builder()
    ti_pool, ti_pool_symbols, ag_pool, _, lattice_report = lattice_builder.build_lattice_pools(builder)
    ti_occ = occupancy_grid(ti0)
    ag_occ = occupancy_grid(ag0)
    ti_smooth = smooth_periodic_xy(ti_occ)
    ag_smooth = smooth_periodic_xy(ag_occ)
    ti_score = sample_scores(ti_pool, ti_smooth)
    ag_score = sample_scores(ag_pool, ag_smooth)

    units = base_counts["Si"]
    ti_target = {"Ti": 3 * units, "Si": units, "C": 2 * units}
    ti_positions, ti_symbols, ti_info = connected_exact_selection(
        builder, ti_pool, ti_pool_symbols, ti_score, ti_target,
        TI_CONNECTIVITY_CUTOFF_A, "Ti3SiC2"
    )
    ti_tree = builder.tree(ti_positions)

    # 按 Ag 自身平滑密度筛选银晶格位点；优先保留高占位密度区域，并避开骨架原子。
    ag_symbols = np.full(len(ag_pool), "Ag")
    eligible_ag = ti_tree.query(ag_pool, workers=2)[0] >= MIN_CROSS_PHASE_DISTANCE_A
    ag_pool_clear = ag_pool[eligible_ag]
    ag_score_clear = ag_score[eligible_ag]
    ag_target = round(units * FORMULA_MASS / MASS["Ag"])
    ag_positions, _, ag_info = connected_exact_selection(
        builder, ag_pool_clear, ag_symbols[eligible_ag], ag_score_clear,
        {"Ag": ag_target}, AG_CONNECTIVITY_CUTOFF_A, "Ag"
    )

    positions = np.vstack((ti_positions, ag_positions))
    symbols = ti_symbols.tolist() + ["Ag"] * len(ag_positions)
    atoms = Atoms(symbols=symbols, positions=positions, cell=common.BOX,
                  pbc=[True, True, False])
    atoms.new_array("phase", np.r_[np.ones(len(ti_positions), np.int32),
                                   np.full(len(ag_positions), 2, np.int32)])
    atoms.new_array("fixed_layer", (positions[:, 2] < 20.0).astype(np.int32))
    atoms.new_array("thermostat_layer", ((positions[:, 2] >= 20.0) & (positions[:, 2] < 40.0)).astype(np.int32))
    atoms.new_array("newton_layer", ((positions[:, 2] >= 40.0) & (positions[:, 2] < 200.0)).astype(np.int32))
    atoms.new_array("sputtering_layer", ((positions[:, 2] >= 200.0) & (positions[:, 2] < 250.0)).astype(np.int32))
    atoms.new_array("evaporation_layer", ((positions[:, 2] >= 250.0) & (positions[:, 2] < 300.0)).astype(np.int32))
    atoms.info.update({
        "status": "UNRELAXED_GEOMETRY_CANDIDATE_STEP06_COARSENED_NETWORKS",
        "thickening_step": 6,
        "coarsening_gaussian_sigma_voxels": SIGMA_VOXELS,
        "coarsening_length_scale_A": SIGMA_VOXELS * float(VOXEL.mean()),
        "phase_atom_counts_preserved_from_step05": True,
        "target_Ag_over_Ti3SiC2_mass_ratio": 1.0,
        "material_region_A": common.SIZE.tolist(),
        "phase_order": "phase 1=Ti3SiC2, phase 2=Ag",
    })

    OUTPUT.mkdir(parents=True, exist_ok=True)
    extxyz = OUTPUT / "ag_ti3sic2_coarsened_networks_step06.extxyz"
    datafile = OUTPUT / "ag_ti3sic2_coarsened_networks_step06.data"
    write(extxyz, atoms, format="extxyz")
    write(datafile, atoms, format="lammps-data", atom_style="atomic",
          specorder=["Ti", "Si", "C", "Ag"], masses=True)
    ti_mask0 = ti_occ > 0
    ag_mask0 = ag_occ > 0
    ti_mask6 = equal_volume_smooth_mask(ti_smooth, ti_occ)
    ag_mask6 = equal_volume_smooth_mask(ag_smooth, ag_occ)
    np.savez_compressed(
        OUTPUT / "coarsened_step06_phase_fields.npz",
        step05_Ti3SiC2_occupancy=ti_occ.astype(np.uint16),
        step05_Ag_occupancy=ag_occ.astype(np.uint16),
        smooth_Ti3SiC2_density=ti_smooth.astype(np.float32),
        smooth_Ag_density=ag_smooth.astype(np.float32),
        step06_equal_volume_Ti3SiC2_mask=ti_mask6.astype(np.uint8),
        step06_equal_volume_Ag_mask=ag_mask6.astype(np.uint8),
        gaussian_sigma_voxels=np.asarray(SIGMA_VOXELS),
        material_dimensions_A=common.SIZE,
        simulation_box_A=common.BOX,
    )

    check = read(extxyz)
    phases = np.asarray(check.arrays["phase"], dtype=np.int32)
    ti_out = check.positions[phases == 1]
    ag_out = check.positions[phases == 2]
    out_symbols = np.asarray(check.get_chemical_symbols())
    counts = {s: int(np.count_nonzero(out_symbols == s)) for s in ("Ti", "Si", "C", "Ag")}
    formula_units = counts["Si"]
    ratio = counts["Ag"] * MASS["Ag"] / (formula_units * FORMULA_MASS)
    ti_connectivity = builder.component_filter(ti_out, cutoff=TI_CONNECTIVITY_CUTOFF_A)[1]
    ag_connectivity = builder.component_filter(ag_out, cutoff=AG_CONNECTIVITY_CUTOFF_A)[1]
    min_cross = float(builder.tree(ti_out).query(ag_out, workers=2)[0].min())
    if counts != base_counts:
        raise RuntimeError(f"粗化后相元素原子数必须保持 step05 数量：{counts} != {base_counts}")
    if ti_connectivity["discarded_atoms"] or ag_connectivity["discarded_atoms"]:
        raise RuntimeError(f"粗化后有断开的相网络：Ti={ti_connectivity}, Ag={ag_connectivity}")
    if abs(ratio - 1.0) > 2e-4:
        raise RuntimeError(f"粗化后 Ag/Ti3SiC2 质量比偏离 1:1：{ratio:.8f}")
    if min_cross < MIN_CROSS_PHASE_DISTANCE_A - 1e-5:
        raise RuntimeError(f"粗化后相间原子距离过近：{min_cross:.6f} Å")

    report = {
        "status": atoms.info["status"],
        "baseline_model": str(BASELINE.relative_to(REPO)),
        "dimensions": {"material_region_A": common.SIZE.tolist(),
                       "simulation_box_A": common.BOX.tolist(), "pbc": [True, True, False]},
        "coarsening": {"method": "对 step05 两相原子占位密度施加周期 x/y 的高斯平滑，再在原 Ti3SiC2 与 Ag 晶格上按平滑密度选回相同的元素原子数",
                       "gaussian_sigma_voxels": SIGMA_VOXELS,
                       "approximate_smoothing_length_A": SIGMA_VOXELS * float(VOXEL.mean()),
                       "Ti3SiC2_candidate_oversampling": TI_CANDIDATE_MULTIPLIER,
                       "Ag_candidate_oversampling": AG_CANDIDATE_MULTIPLIER,
                       "minimum_cross_phase_distance_A": MIN_CROSS_PHASE_DISTANCE_A,
                       "limitations": "此步骤将原子从细枝区域重新分配到较粗的平滑网络；相总原子数和质量保持不变，不是增加材料质量。"},
        "lattice": lattice_report,
        "atoms": len(check),
        "element_counts": counts,
        "Ag_over_Ti3SiC2_mass_ratio": float(ratio),
        "minimum_cross_phase_distance_A": min_cross,
        "connectivity": {"Ti3SiC2": ti_connectivity, "Ag": ag_connectivity,
                         "Ti3SiC2_selection": ti_info, "Ag_selection": ag_info},
        "equal_volume_smooth_mask_thickness_proxy": {
            "Ti3SiC2_step05": thickness_proxy(ti_mask0),
            "Ti3SiC2_step06": thickness_proxy(ti_mask6),
            "Ag_step05": thickness_proxy(ag_mask0),
            "Ag_step06": thickness_proxy(ag_mask6),
        },
        "outputs": {"extxyz": extxyz.name, "lammps_data": datafile.name,
                    "phase_fields": "coarsened_step06_phase_fields.npz"},
        "caveats": [
            "材料区与模拟盒物理尺寸未变化，模型仍为纳米尺寸。",
            "该粗化通过原子位置重分配构造，属于未弛豫几何候选。",
            "没有运行能量最小化、冷态稳定性检查或 TTM/MD 加热。",
            "相的原子总数保持 step05，因此 1:1 质量比精确保留；局部枝干增粗伴随细枝平滑与重新分配。",
        ],
    }
    report_path = OUTPUT / "coarsened_step06_report.json"
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"status": report["status"], "atoms": len(check), "counts": counts,
                      "mass_ratio": ratio, "minimum_cross_phase_distance_A": min_cross,
                      "connectivity": {"Ti3SiC2": ti_connectivity, "Ag": ag_connectivity},
                      "thickness_proxy": report["equal_volume_smooth_mask_thickness_proxy"]},
                     ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
