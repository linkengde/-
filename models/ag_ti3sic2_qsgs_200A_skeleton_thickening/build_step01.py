#!/usr/bin/env python3
"""在原始 20 nm Ag/Ti3SiC2 原子模型上做一档 1 Å/侧的骨架加粗。"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
from ase import Atoms
from ase.data import atomic_masses, atomic_numbers
from ase.io import read, write
from scipy import ndimage as ndi

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
BASE = REPO / "models/ag_ti3sic2_qsgs_200A_ttm_finer"
BASE_MODEL = BASE / "ag_ti3sic2_qsgs.extxyz"
BASE_PHASE = BASE / "qsgs_morphology.npz"
BASE_REPORT = BASE / "geometry_report.json"
BUILDER_PATH = REPO / "models/ag_ti3sic2_qsgs_200A_spacing10x_atomic/build_atomic_model.py"
OUT = HERE / "output"

SIZE = np.array([200.0, 200.0, 200.0])
BOX = np.array([200.0, 200.0, 300.0])
GRID = np.array([80, 80, 80])
VOXEL = SIZE / GRID
THICKEN_A = 1.0
CLEARANCE_A = 2.2
MASS = {s: float(atomic_masses[atomic_numbers[s]]) for s in ("Ag", "Ti", "Si", "C")}


def load_builder():
    """重用原模型的晶格、化学计量和周期连通性实现。"""
    spec = importlib.util.spec_from_file_location("qsgs_atomic_builder", BUILDER_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"无法加载原模型建模工具：{BUILDER_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module.load_original_builder(), module


def signed_distance(mask: np.ndarray) -> np.ndarray:
    """得到 x/y 周期、z 开放的近似带符号距离场，单位 Å。"""
    tiled = np.tile(mask, (3, 3, 1))
    inside = ndi.distance_transform_edt(tiled, sampling=VOXEL)
    outside = ndi.distance_transform_edt(~tiled, sampling=VOXEL)
    # 距离变换位于体素中心；减半个体素边长，将零面移至两相界面。
    field = inside - outside - float(VOXEL.min() / 2.0)
    return field


def sample_field(field: np.ndarray, positions: np.ndarray) -> np.ndarray:
    """在原子坐标采样距离场；索引按体素中心和周期 x/y 平移校正。"""
    grid = positions / VOXEL
    coordinates = np.vstack((grid[:, 0] + 79.5, grid[:, 1] + 79.5, grid[:, 2] - 0.5))
    return ndi.map_coordinates(field, coordinates, order=1, mode="nearest", prefilter=False)


def connected_phase(base, positions: np.ndarray, symbols: np.ndarray,
                    mask: np.ndarray, ti_phase: bool):
    """保留最大连通相，并对 Ti3SiC2 反复修整为 3:1:2。"""
    cleanup = []
    if ti_phase:
        units = 0
        for _ in range(8):
            keep, info = base.component_filter(positions)
            cleanup.append(info)
            positions, symbols = positions[keep], symbols[keep]
            positions, symbols, units, trimmed = base.stoichiometry(positions, symbols, mask)
            if info["discarded_atoms"] == 0 and trimmed == 0:
                break
        else:
            raise RuntimeError("加粗后的 Ti3SiC2 连通性和化学计量未稳定")
        return positions, symbols, units, cleanup

    keep, info = base.component_filter(positions, cutoff=3.15)
    cleanup.append(info)
    positions = positions[keep]
    if not keep.all():
        keep, info = base.component_filter(positions, cutoff=3.15)
        cleanup.append(info)
        positions = positions[keep]
    if not len(positions):
        raise RuntimeError("加粗后 Ag 连通相为空")
    return positions, symbols, 0, cleanup


def main() -> None:
    for path in (BASE_MODEL, BASE_PHASE, BASE_REPORT, BUILDER_PATH):
        if not path.exists():
            raise FileNotFoundError(f"缺少输入文件：{path}")
    original = read(BASE_MODEL)
    original_report = json.loads(BASE_REPORT.read_text())
    original_counts = {s: int(np.count_nonzero(np.asarray(original.get_chemical_symbols()) == s))
                       for s in ("Ti", "Si", "C", "Ag")}
    if len(original) != 437235 or not np.allclose(original.cell.lengths(), BOX):
        raise ValueError("原始基准原子数或模拟盒尺寸不符合记录")
    if original_counts != {"Ti": 167859, "Si": 55953, "C": 111906, "Ag": 101517}:
        raise ValueError(f"原始元素计数不符合基准模型：{original_counts}")

    base, lattice_builder = load_builder()
    mask = np.load(BASE_PHASE)["skeleton"].astype(bool)
    if mask.shape != tuple(GRID):
        raise ValueError(f"原始 QSGS 网格形状异常：{mask.shape}")
    field = signed_distance(mask)
    ti_pool, ti_pool_symbols, ag_pool, _, lattice = lattice_builder.build_lattice_pools(base)

    # 原相场内的 Ti 位点全保留；仅将相界面外 1 Å 内的晶格位点并入骨架。
    ti_inside = base.at_mask(mask, ti_pool)
    ti_distance = sample_field(field, ti_pool)
    ti_selected = ti_inside | ((~ti_inside) & (ti_distance >= -THICKEN_A))
    ti_positions = ti_pool[ti_selected]
    ti_symbols = ti_pool_symbols[ti_selected]
    ti_positions, ti_symbols, units, ti_cleanup = connected_phase(
        base, ti_positions, ti_symbols, mask, ti_phase=True
    )

    # Ag 只从未扩张相区选位点，并移除界面过近的候选原子。
    ag_outside = ~base.at_mask(mask, ag_pool)
    ag_distance = sample_field(field, ag_pool)
    ag_selected = ag_outside & (ag_distance < -THICKEN_A)
    ag_positions = ag_pool[ag_selected]
    if len(ag_positions):
        separation = base.tree(ti_positions).query(ag_positions, workers=2)[0]
        ag_positions = ag_positions[separation >= CLEARANCE_A]
    ag_available_before_connectivity = int(len(ag_positions))
    ag_positions, _, _, ag_cleanup = connected_phase(
        base, ag_positions, np.full(len(ag_positions), "Ag"), mask, ti_phase=False
    )

    # 写出实际可用的连通 Ag 位点；加粗后若不足以保持严格等质量比，如实记录。
    formula_mass = 3 * MASS["Ti"] + MASS["Si"] + 2 * MASS["C"]
    ag_mass_ratio = len(ag_positions) * MASS["Ag"] / (units * formula_mass)
    positions = np.vstack((ti_positions, ag_positions))
    symbols = ti_symbols.tolist() + ["Ag"] * len(ag_positions)
    atoms = Atoms(symbols=symbols, positions=positions, cell=BOX, pbc=[True, True, False])
    atoms.new_array("phase", np.r_[np.ones(len(ti_positions), np.int32),
                                   np.full(len(ag_positions), 2, np.int32)])
    atoms.new_array("fixed_layer", (positions[:, 2] < 20.0).astype(np.int32))
    atoms.new_array("thermostat_layer", ((positions[:, 2] >= 20.0) & (positions[:, 2] < 40.0)).astype(np.int32))
    atoms.new_array("newton_layer", ((positions[:, 2] >= 40.0) & (positions[:, 2] < 200.0)).astype(np.int32))
    atoms.new_array("sputtering_layer", ((positions[:, 2] >= 200.0) & (positions[:, 2] < 250.0)).astype(np.int32))
    atoms.new_array("evaporation_layer", ((positions[:, 2] >= 250.0) & (positions[:, 2] < 300.0)).astype(np.int32))
    atoms.info.update({
        "status": "UNRELAXED_GEOMETRY_CANDIDATE_STEP01_SKELETON_THICKENED",
        "thickening_each_side_A": THICKEN_A,
        "material_region_A": SIZE.tolist(),
        "phase_order": "phase 1=Ti3SiC2, phase 2=Ag",
    })

    OUT.mkdir(parents=True, exist_ok=True)
    extxyz = OUT / "ag_ti3sic2_skeleton_thickened_step01.extxyz"
    datafile = OUT / "ag_ti3sic2_skeleton_thickened_step01.data"
    write(extxyz, atoms, format="extxyz")
    write(datafile, atoms, format="lammps-data", atom_style="atomic",
          specorder=["Ti", "Si", "C", "Ag"], masses=True)
    np.savez_compressed(
        OUT / "thickened_step01_phase_field.npz",
        original_Ti3SiC2_mask=mask.astype(np.uint8),
        signed_distance_A=field.astype(np.float32),
        thickness_increase_each_side_A=np.asarray(THICKEN_A),
        material_dimensions_A=SIZE,
        simulation_box_A=BOX,
    )

    final_counts = {s: int(np.count_nonzero(np.asarray(symbols) == s)) for s in ("Ti", "Si", "C", "Ag")}
    phase_atom_changes = {}
    for phase_id, label in ((1, "Ti3SiC2"), (2, "Ag")):
        old_pos = original.positions[original.arrays["phase"] == phase_id]
        new_pos = positions[atoms.arrays["phase"] == phase_id]
        old_to_new = base.tree(new_pos).query(old_pos, workers=2)[0] < 1.0e-5
        new_to_old = base.tree(old_pos).query(new_pos, workers=2)[0] < 1.0e-5
        phase_atom_changes[label] = {
            "baseline_atoms_retained": int(old_to_new.sum()),
            "baseline_atoms_removed": int((~old_to_new).sum()),
            "new_atoms_added": int((~new_to_old).sum()),
        }
    report = {
        "status": "UNRELAXED_GEOMETRY_CANDIDATE_STEP01_SKELETON_THICKENED",
        "description": "在原纳米原子模型同尺寸材料区内，将 Ti3SiC2 相界面向 Ag 侧外扩 1.0 Å；Ag 只保留扩张后空隙中的连通位点。",
        "baseline_model": str(BASE_MODEL.relative_to(REPO)),
        "baseline_atoms": len(original),
        "baseline_element_counts": original_counts,
        "dimensions": {"material_region_A": SIZE.tolist(), "simulation_box_A": BOX.tolist(),
                       "pbc": [True, True, False]},
        "thickening": {"offset_each_side_A": THICKEN_A,
                       "method": "QSGS voxel-center signed-distance interpolation sampled at Ti3SiC2 and Ag crystal sites",
                       "qsgs_grid": GRID.tolist(), "qsgs_voxel_size_A": VOXEL.tolist()},
        "final_atoms": len(atoms),
        "final_element_counts": final_counts,
        "Ti3SiC2_formula_units": int(units),
        "Ag_over_Ti3SiC2_mass_ratio": float(ag_mass_ratio),
        "Ag_available_before_connectivity_filter": ag_available_before_connectivity,
        "Ti3SiC2_connectivity_cleanup": ti_cleanup,
        "Ag_connectivity_cleanup": ag_cleanup,
        "phase_atom_changes_vs_baseline": phase_atom_changes,
        "Ti3SiC2_atoms_added_vs_baseline_count": int(sum(final_counts[s] for s in ("Ti", "Si", "C")) - sum(original_counts[s] for s in ("Ti", "Si", "C"))),
        "limitations": [
            "此版本保持 200×200×200 Å 材料区及 200×200×300 Å 模拟盒，物理尺寸仍为纳米级。",
            "这是未弛豫的几何原子模型；还没有进行势函数力学稳定性或 TTM/MD 加热验证。",
            "Ag 按扩张后的可用晶格位点填充；若骨架继续加粗，Ag/Ti3SiC2 质量比会继续变化。",
            "形貌可作为后续粗化/微米风格结构的比较版本，但不表示实体模型尺寸已经达到微米。",
        ],
        "outputs": {"extxyz": extxyz.name, "lammps_data": datafile.name,
                    "phase_field": "thickened_step01_phase_field.npz", "report": "thickened_step01_report.json"},
    }
    (OUT / "thickened_step01_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"status": report["status"], "atoms": len(atoms),
                      "counts": final_counts, "mass_ratio": ag_mass_ratio,
                      "material_A": SIZE.tolist(), "box_A": BOX.tolist(),
                      "thicken_each_side_A": THICKEN_A,
                      "Ti3SiC2_cleanup": ti_cleanup, "Ag_cleanup": ag_cleanup},
                     ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
