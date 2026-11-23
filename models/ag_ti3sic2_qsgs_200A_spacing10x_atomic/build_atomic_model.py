#!/usr/bin/env python3
"""从已建纳米原子模型的 QSGS 形貌构建枝间距约放大十倍的原子级双相模型。"""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
from ase import Atoms
from ase.build import bulk, make_supercell
from ase.data import atomic_masses, atomic_numbers
from ase.io import read, write
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
BASE = REPO / "models/ag_ti3sic2_qsgs_200A_ttm_finer"
BASE_SCRIPT = BASE / "build_model.py"
BASE_ATOM_DATA = BASE / "ag_ti3sic2_qsgs.data"
BASE_EXTXYZ = BASE / "ag_ti3sic2_qsgs.extxyz"
BASE_REPORT = BASE / "geometry_report.json"
MORPHOLOGY = BASE / "qsgs_morphology.npz"
CRYSTAL_CIF = REPO / "research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures/source/Ti3SiC2_COD_9009647.cif"
OUT = HERE / "output"

SIZE = np.array([200.0, 200.0, 200.0])
BOX = np.array([200.0, 200.0, 300.0])
GRID = np.array([80, 80, 80])
VOXEL_A = SIZE / GRID
VACUUM = 100.0
FIXED_A = 20.0
THERMOSTAT_TOP_A = 40.0
SPUTTER_TOP_A = 250.0
EVAPORATION_TOP_A = 300.0
CLEARANCE_A = 2.2
SOURCE_CROP_START = (24, 53, 41)
SOURCE_CROP_SIZE = 8
SPATIAL_SCALE = 10
MASK_THRESHOLD = 0.5

MASS = {s: float(atomic_masses[atomic_numbers[s]]) for s in ("Ag", "Ti", "Si", "C")}
FORMULA_MASS = 3 * MASS["Ti"] + MASS["Si"] + 2 * MASS["C"]


def load_original_builder():
    """载入原始纳米模型生成器的晶格和连通性函数。"""
    spec = importlib.util.spec_from_file_location("baseline_qsgs_builder", BASE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"无法载入原始建模脚本：{BASE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def make_coarse_mask(base) -> tuple[np.ndarray, np.ndarray, list[int]]:
    """读取原纳米原子模型所用的 Ti 相标签场，再将 20 Å 子域放大十倍。"""
    baseline = read(BASE_EXTXYZ)
    if len(baseline) != 437235 or "phase" not in baseline.arrays:
        raise ValueError("纳米基准原子文件的原子数或 phase 标签不符合预期")
    if not np.allclose(baseline.cell.lengths(), BOX, atol=1e-6):
        raise ValueError(f"纳米基准模拟盒尺寸异常：{baseline.cell.lengths()}")
    if not np.count_nonzero(baseline.arrays["phase"] == 1):
        raise ValueError("纳米基准原子模型中没有 Ti₃SiC₂ 相标签")
    # 该相标签场正是生成已核验纳米原子模型时使用的几何拓扑；直接取其局部窗，
    # 避免晶格采样噪声在 8³ 源网格上改写分支连通关系。
    source = np.load(MORPHOLOGY)["skeleton"].astype(np.float32)
    if source.shape != tuple(GRID):
        raise ValueError(f"原始模型 QSGS 网格异常：{source.shape}")
    x, y, z = SOURCE_CROP_START
    n = SOURCE_CROP_SIZE
    crop = source[x:x+n, y:y+n, z:z+n].astype(np.float32)
    if crop.shape != (n, n, n):
        raise ValueError(f"QSGS 子域尺寸异常：{crop.shape}")
    scaled = ndi.zoom(crop, SPATIAL_SCALE, order=1, mode="nearest", prefilter=False)
    if scaled.shape != tuple(GRID):
        raise ValueError(f"十倍映射后的网格异常：{scaled.shape}")
    mask = scaled >= MASK_THRESHOLD
    # 原模型同样用周期 x/y、开放 z 的连通性清理，保留 Ti 与 Ag 双连续通路。
    mask, changes = base.clean_mask(mask)
    return mask, crop, changes


def build_lattice_pools(base) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict]:
    """用原模型相同的 Ti3SiC2 晶向和 Ag 晶格准备候选原子位点。"""
    cell = make_supercell(read(CRYSTAL_CIF), [[1, 0, 0], [1, 2, 0], [0, 0, 1]])
    ti_repeats = tuple(int(v) for v in np.ceil(SIZE / cell.cell.lengths()))
    ti = cell.repeat(ti_repeats)
    ti_lengths = ti.cell.lengths().copy()
    ti_target = ti_lengths.copy()
    ti_target[:2] = SIZE[:2]
    ti.set_cell(ti_target, scale_atoms=True)
    ti_strain = (ti_target / ti_lengths - 1.0).tolist()
    ti_valid = ti.positions[:, 2] < SIZE[2]
    ti_pos = ti.positions[ti_valid].copy()
    ti_symbols = np.asarray(ti.get_chemical_symbols())[ti_valid]
    del ti

    ag = bulk("Ag", "fcc", a=4.09, cubic=True)
    ag_repeats = tuple(int(v) for v in np.ceil(SIZE / ag.cell.lengths()))
    ag = ag.repeat(ag_repeats)
    ag_lengths = ag.cell.lengths().copy()
    ag_target = ag_lengths.copy()
    ag_target[:2] = SIZE[:2]
    ag.set_cell(ag_target, scale_atoms=True)
    ag_strain = (ag_target / ag_lengths - 1.0).tolist()
    ag_pos = ag.positions[ag.positions[:, 2] < SIZE[2]].copy()
    del ag

    for positions in (ti_pos, ag_pos):
        positions[:, :2] %= SIZE[:2]
    lattice_report = {
        "Ti3SiC2_source_cif": str(CRYSTAL_CIF.relative_to(REPO)),
        "Ti3SiC2_repeats": list(ti_repeats),
        "Ti3SiC2_strain_xyz": ti_strain,
        "Ti3SiC2_orientation": "c axis along z, matching the original nanomodel",
        "Ag_lattice_constant_A": 4.09,
        "Ag_repeats": list(ag_repeats),
        "Ag_strain_xyz": ag_strain,
    }
    return ti_pos, ti_symbols, ag_pos, ag_target, lattice_report


def make_atomic_model(base, mask: np.ndarray):
    """将放大后的完整相场填入 Ti3SiC2 与 Ag 原子，并清理两相连通性。"""
    ti_pool, ti_pool_symbols, ag_pool, _, lattice_report = build_lattice_pools(base)
    ti_pos = ti_pool[base.at_mask(mask, ti_pool)]
    ti_symbols = ti_pool_symbols[base.at_mask(mask, ti_pool)]
    ti_pos, ti_symbols, units, stoich_trim = base.stoichiometry(ti_pos, ti_symbols, mask)

    # 移除孤立 Ti 晶粒并重整边界化学计量，直到 Ti3SiC2 相成为单一连通网络。
    ti_cleanup = []
    for _ in range(8):
        keep, info = base.component_filter(ti_pos)
        ti_cleanup.append(info)
        ti_pos, ti_symbols = ti_pos[keep], ti_symbols[keep]
        ti_pos, ti_symbols, units, trimmed = base.stoichiometry(ti_pos, ti_symbols, mask)
        if info["discarded_atoms"] == 0 and trimmed == 0:
            break
    else:
        raise RuntimeError("Ti3SiC2 原子连通性/化学计量未能稳定")

    # Ag 填入互补相区，保持最小界面间距并使 Ag 相连通。
    ag_pos = ag_pool[~base.at_mask(mask, ag_pool)]
    distance_to_ti, _ = base.tree(ti_pos).query(ag_pos, workers=2)
    ag_pos = ag_pos[distance_to_ti >= CLEARANCE_A]
    available_mass_ratio = len(ag_pos) * MASS["Ag"] / (units * FORMULA_MASS)
    desired_ag = round(units * FORMULA_MASS / MASS["Ag"])
    if len(ag_pos) < desired_ag:
        raise RuntimeError(f"可用连通 Ag 位点不足：{len(ag_pos)} < {desired_ag}")

    ag_original = ag_pos.copy()
    keep, ag_info = base.component_filter(ag_pos, cutoff=3.15)
    ag_cleanup = [ag_info]
    ag_pos = ag_pos[keep]
    if len(ag_pos) < desired_ag:
        raise RuntimeError("Ag 连通主相小于等质量目标")
    ti_distance, _ = base.tree(ti_pos).query(ag_pos, workers=2)
    trim = len(ag_pos) - desired_ag
    if trim:
        keep = np.ones(len(ag_pos), dtype=bool)
        keep[np.argsort(ti_distance, kind="stable")[:trim]] = False
        ag_pos = ag_pos[keep]
    keep, ag_info = base.component_filter(ag_pos, cutoff=3.15)
    if not keep.all():
        # 若界面修剪断开银相，则从剩余连通 Ag 的前沿向外补回内部位点。
        ag_pos = ag_pos[keep]
        current_tree = base.tree(ag_pos)
        occupied = current_tree.query(ag_original, workers=2)[0] < 1.0e-7
        depth = base.tree(ti_pos).query(ag_original, workers=2)[0]
        while occupied.sum() < desired_ag:
            candidates = np.flatnonzero(~occupied)
            neighbor_d = base.tree(ag_original[occupied]).query(
                ag_original[candidates], workers=2
            )[0]
            candidates = candidates[neighbor_d <= 3.15]
            if not len(candidates):
                raise RuntimeError("Ag 连通生长前沿为空")
            candidates = candidates[np.argsort(-depth[candidates], kind="stable")]
            occupied[candidates[:desired_ag - int(occupied.sum())]] = True
        ag_pos = ag_original[occupied]
        keep, ag_info = base.component_filter(ag_pos, cutoff=3.15)
        if not keep.all():
            raise RuntimeError("Ag 相最终连通性失败")
    ag_cleanup.append(ag_info)

    keep, ti_info = base.component_filter(ti_pos)
    if not keep.all():
        raise RuntimeError("Ti3SiC2 相最终连通性失败")

    positions = np.vstack((ti_pos, ag_pos))
    symbols = ti_symbols.tolist() + ["Ag"] * len(ag_pos)
    atoms = Atoms(
        symbols=symbols,
        positions=positions,
        cell=[SIZE[0], SIZE[1], BOX[2]],
        pbc=[True, True, False],
    )
    atoms.new_array("phase", np.r_[np.ones(len(ti_pos), np.int32), np.full(len(ag_pos), 2, np.int32)])
    atoms.new_array("fixed_layer", (positions[:, 2] < FIXED_A).astype(np.int32))
    atoms.new_array("thermostat_layer", ((positions[:, 2] >= FIXED_A) & (positions[:, 2] < THERMOSTAT_TOP_A)).astype(np.int32))
    atoms.new_array("newton_layer", ((positions[:, 2] >= THERMOSTAT_TOP_A) & (positions[:, 2] < SIZE[2])).astype(np.int32))
    atoms.new_array("sputtering_layer", ((positions[:, 2] >= SIZE[2]) & (positions[:, 2] < SPUTTER_TOP_A)).astype(np.int32))
    atoms.new_array("evaporation_layer", ((positions[:, 2] >= SPUTTER_TOP_A) & (positions[:, 2] < EVAPORATION_TOP_A)).astype(np.int32))
    atoms.info.update({
        "status": "UNRELAXED_ATOMISTIC_10X_SPACING_GEOMETRY",
        "material_thickness_A": float(SIZE[2]),
        "vacuum_A": VACUUM,
        "nominal_spacing_scale": SPATIAL_SCALE,
        "phase_order": "phase 1=Ti3SiC2, phase 2=Ag",
    })

    mass_ratio = len(ag_pos) * MASS["Ag"] / (units * FORMULA_MASS)
    counts = {s: int(np.sum(np.asarray(symbols) == s)) for s in ("Ti", "Si", "C", "Ag")}
    # 计算相内点密度，用原子坐标重建 Ti3SiC2 的可视化相边界。
    ti_pool_density = len(ti_pool) / float(np.prod(SIZE))
    ag_pool_density = len(ag_pool) / float(np.prod(SIZE))
    ti_volume = len(ti_pos) / ti_pool_density
    ag_volume = len(ag_pos) / ag_pool_density
    material_volume = float(np.prod(SIZE))
    ti_fraction_box = float(ti_volume / material_volume)
    ag_fraction_box = float(ag_volume / material_volume)
    atomic_volume_fraction_ti = float(ti_volume / (ti_volume + ag_volume))
    details = {
        "atoms": int(len(atoms)),
        "counts": counts,
        "Ti3SiC2_formula_units": int(units),
        "Ti3SiC2_stoichiometry_trimmed_atoms": int(stoich_trim),
        "mass_ratio_Ag_over_Ti3SiC2": float(mass_ratio),
        "available_Ag_mass_ratio_before_trim": float(available_mass_ratio),
        "Ti3SiC2_volume_fraction_of_box_from_lattice_site_densities": ti_fraction_box,
        "Ag_volume_fraction_of_box_from_lattice_site_densities": ag_fraction_box,
        "void_volume_fraction_of_box_from_lattice_site_densities": float(max(0.0, 1.0 - ti_fraction_box - ag_fraction_box)),
        "Ti3SiC2_fraction_of_occupied_solid_volume_from_lattice_site_densities": atomic_volume_fraction_ti,
        "Ti3SiC2_connectivity": ti_info,
        "Ag_connectivity": ag_info,
        "Ti3SiC2_cleanup": ti_cleanup,
        "Ag_cleanup": ag_cleanup,
        "lattice": lattice_report,
    }
    density = {
        "Ti3SiC2_atoms_per_A3": float(ti_pool_density),
        "Ag_atoms_per_A3": float(ag_pool_density),
        "Ti3SiC2_volume_A3": float(ti_volume),
        "Ag_volume_A3": float(ag_volume),
    }
    return atoms, mask, details, density


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for path in (BASE_ATOM_DATA, BASE_EXTXYZ, BASE_REPORT, MORPHOLOGY, CRYSTAL_CIF):
        if not path.exists():
            raise FileNotFoundError(f"缺少纳米基准模型输入：{path}")
    baseline_report = json.loads(BASE_REPORT.read_text())
    if baseline_report.get("atoms") != 437235:
        raise ValueError(f"纳米基准原子数不是预期值：{baseline_report.get('atoms')}")
    base = load_original_builder()
    mask, crop, changes = make_coarse_mask(base)
    atoms, mask, atomic_report, density_report = make_atomic_model(base, mask)

    extxyz = OUT / "ag_ti3sic2_spacing10x.extxyz"
    datafile = OUT / "ag_ti3sic2_spacing10x.data"
    write(extxyz, atoms, format="extxyz")
    write(datafile, atoms, format="lammps-data", atom_style="atomic",
          specorder=["Ti", "Si", "C", "Ag"], masses=True)
    np.savez_compressed(
        OUT / "spacing10x_atomic_phase_field.npz",
        Ti3SiC2_mask=mask.astype(np.uint8),
        Ag_mask=(~mask).astype(np.uint8),
        source_crop=crop.astype(np.uint8),
        grid_xyz=GRID,
        voxel_spacing_A=VOXEL_A,
        material_dimensions_A=SIZE,
        simulation_box_A=BOX,
    )
    ratio_gap = 9.3719
    report = {
        "status": "UNRELAXED_ATOMISTIC_SPACING10X_GEOMETRY",
        "construction_order": [
            "inspect_existing_437235_atom_nanomodel",
            "scale_20A_QSGS_skeleton_window_to_200A_material_region",
            "fill_Ti3SiC2_and_Ag_crystal_lattice_sites",
            "extract_skeleton_from_final_Ti3SiC2_atoms_in_separate_step",
        ],
        "baseline_model": {
            "file": str(BASE_EXTXYZ.relative_to(REPO)),
            "atoms": int(baseline_report["atoms"]),
            "material_dimensions_A": baseline_report["material_dimensions_A"],
            "simulation_box_A": baseline_report["simulation_box_A"],
            "sha256_extxyz": hashlib.sha256(BASE_EXTXYZ.read_bytes()).hexdigest(),
        },
        "source_qsgs": {
            "file": str(MORPHOLOGY.relative_to(REPO)),
            "sha256": hashlib.sha256(MORPHOLOGY.read_bytes()).hexdigest(),
            "shape_source": "QSGS Ti3SiC2 phase labels underlying the baseline 437235-atom model; baseline extxyz phase labels, cell and atom count are checked before use",
            "baseline_Ti3SiC2_atom_count": int(sum(baseline_report["counts"][s] for s in ("Ti", "Si", "C"))),
            "crop_start_voxel_xyz": list(SOURCE_CROP_START),
            "crop_size_voxels_xyz": [SOURCE_CROP_SIZE] * 3,
            "crop_extent_A": [20.0] * 3,
            "spatial_scale_factor": SPATIAL_SCALE,
            "threshold": MASK_THRESHOLD,
            "mask_volume_fraction_after_two_phase_connectivity_cleanup": float(mask.mean()),
            "cleanup_changes_by_pass": changes,
        },
        "dimensions": {
            "material_region_A": SIZE.tolist(),
            "simulation_box_A": BOX.tolist(),
            "pbc": [True, True, False],
            "fixed_layer_A": [0.0, FIXED_A],
            "thermostat_layer_A": [FIXED_A, THERMOSTAT_TOP_A],
            "newton_layer_A": [THERMOSTAT_TOP_A, SIZE[2]],
            "sputtering_layer_A": [SIZE[2], SPUTTER_TOP_A],
            "evaporation_layer_A": [SPUTTER_TOP_A, EVAPORATION_TOP_A],
        },
        "atomic_model": atomic_report,
        "phase_lattice_number_densities": density_report,
        "nominal_target_source_spacing_A": ratio_gap,
        "limitations": [
            "This is an unrelaxed atomistic geometry; no MD or TTM was run.",
            "The 10x morphology is made by mapping a 20 A region of the existing QSGS field onto the unchanged 200 A material region.",
            "The Ag/Ti3SiC2 interface and the current candidate potentials have not passed mechanical stability validation.",
            "The extracted skeleton is a visualization product generated from the final Ti3SiC2 atom coordinates, not an additional atom species.",
        ],
        "outputs": {
            "lammps_data": datafile.name,
            "ovito_extxyz": extxyz.name,
            "phase_field": "spacing10x_atomic_phase_field.npz",
            "report": "atomic_model_report.json",
        },
    }
    (OUT / "atomic_model_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({
        "status": report["status"],
        "atoms": atomic_report["atoms"],
        "counts": atomic_report["counts"],
        "material_A": SIZE.tolist(),
        "box_A": BOX.tolist(),
        "Ti_fraction_voxel": round(float(mask.mean()), 6),
        "mass_ratio_Ag_over_Ti3SiC2": round(atomic_report["mass_ratio_Ag_over_Ti3SiC2"], 10),
        "available_Ag_mass_ratio": round(atomic_report["available_Ag_mass_ratio_before_trim"], 4),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
