#!/usr/bin/env python3
"""整理经过审查的候选标签，生成仅供探索性拟合的数据集。"""
from __future__ import annotations
import hashlib, json
from pathlib import Path
import numpy as np
from ase.io import read, write

BASE = Path('/tmp/mace-fit-data-supplement-20261123-repo/archive/mace-fit-ace-exploratory-20261119')
SUPP = Path('/tmp/mace-fit-data-supplement-20261123-repo/archive/mace-fit-data-supplement-20261123')
OUT = Path(__file__).resolve().parent
ELEMENTS = {'Ag', 'C', 'Si', 'Ti'}

# 将晶胞、边界和按元素排序后的坐标纳入哈希，检测原子顺序变化后的精确重复。
def geometry_hash(atoms) -> str:
    h = hashlib.sha256()
    cell = np.asarray(atoms.cell, dtype=float).round(8)
    pbc = np.asarray(atoms.pbc, dtype=np.uint8)
    if np.linalg.det(cell) > 1e-12:
        q = atoms.get_scaled_positions(wrap=False)
        q[:, pbc.astype(bool)] %= 1.0
    else:
        q = np.asarray(atoms.positions, dtype=float)
    records = sorted((int(z), *np.asarray(pos, dtype=float).round(8))
                     for z, pos in zip(atoms.numbers, q))
    h.update(cell.astype('<f8').tobytes())
    h.update(pbc.tobytes())
    h.update(np.asarray(records, dtype='<f8').tobytes())
    return h.hexdigest()

# 为统一标签字段，将来源输出中的能量和力映射到训练程序使用的 REF 键。
def normalize_label(atoms, source: str, family: str):
    a = atoms.copy()
    energy = a.info.get('REF_energy', a.info.get('PW_PBE_energy_eV'))
    forces = a.arrays.get('REF_forces', a.arrays.get('PW_PBE_forces'))
    if energy is None or forces is None:
        raise ValueError(f'{source}: 缺少 DFT 能量或力')
    forces = np.asarray(forces, dtype=float)
    if forces.shape != (len(a), 3) or not np.isfinite(forces).all():
        raise ValueError(f'{source}: 力数组不合法')
    energy = float(energy)
    if not np.isfinite(energy):
        raise ValueError(f'{source}: 能量不是有限值')
    if set(a.get_chemical_symbols()) - ELEMENTS:
        raise ValueError(f'{source}: 出现非 Ag/C/Si/Ti 元素')
    if np.linalg.det(np.asarray(a.cell)) <= 0:
        raise ValueError(f'{source}: 晶胞体积非正')
    a.info['REF_energy'] = energy
    a.arrays['REF_forces'] = forces.copy()
    a.info['exploratory_source_family'] = family
    a.info['exploratory_data_role'] = 'training_candidate_unassigned'
    a.info['source_method'] = a.info.get(
        'source_method', 'GPAW 26.7.0 PBE/PW500 eV, Fermi-Dirac sigma=0.10 eV; source mesh per manifest')
    return a

# 原始数据不修改；本目录生成独立的派生 train/dev。
train = read(BASE / 'train.extxyz', index=':')
dev = read(BASE / 'dev.extxyz', index=':')
assert (len(train), len(dev)) == (71, 7)

supplement_specs = [
    (SUPP / 'ag_solid_controls/Ag_solid_controls_train_candidate_unassigned.extxyz',
     'Ag32_COD9009647_pure_solid_controls'),
    (SUPP / 'ag_liquid_first_label/runs/attempt_001/training_candidate_unassigned.extxyz',
     'Ag_EAM_WMH_seed20261123_liquid_frame'),
]
# 同一 Ti3SiC2 晶体来源的四个响应构型作为一个训练家族整体处理。
ti_dir = SUPP / 'ti3sic2_force_response_controls'
ti_paths = sorted(ti_dir.glob('*/*_PW_PBE.extxyz'))
if len(ti_paths) != 4:
    raise RuntimeError(f'预期4个Ti3SiC2响应标签，实际{len(ti_paths)}')
new_frames = []
source_records = []
for path, family in supplement_specs:
    frames = read(path, index=':')
    for i, atoms in enumerate(frames):
        label = normalize_label(atoms, str(path), family)
        new_frames.append(label)
        source_records.append({'path': str(path), 'frame': i, 'family': family})
for path in ti_paths:
    atoms = read(path)
    family = 'Ti3SiC2_COD9009647_force_response_controls'
    label = normalize_label(atoms, str(path), family)
    label.info['config_type'] = path.stem
    new_frames.append(label)
    source_records.append({'path': str(path), 'frame': 0, 'family': family})

train_hashes = {geometry_hash(a): i for i, a in enumerate(train)}
dev_hashes = {geometry_hash(a): i for i, a in enumerate(dev)}
accepted, skipped = [], []
seen = set(train_hashes) | set(dev_hashes)
for atoms, rec in zip(new_frames, source_records):
    gh = geometry_hash(atoms)
    if gh in seen:
        skipped.append({**rec, 'geometry_sha256': gh,
                        'reason': 'exact geometry already exists in base train/dev or earlier candidate'})
        continue
    seen.add(gh)
    accepted.append(atoms)

# 开发集保持原样。新增标签来自若干与开发样本来源相关的晶体家族，故结果仅作探索性诊断。
train_out = train + accepted
write(OUT / 'train.extxyz', train_out, format='extxyz')
write(OUT / 'dev.extxyz', dev, format='extxyz')

train_final_hashes = [geometry_hash(a) for a in train_out]
dev_final_hashes = [geometry_hash(a) for a in dev]
intersections = sorted(set(train_final_hashes) & set(dev_final_hashes))
if intersections:
    raise RuntimeError(f'派生 train/dev 存在精确几何重复: {len(intersections)}')

new_summary = []
for a in accepted:
    f = np.asarray(a.arrays['REF_forces'])
    new_summary.append({
        'config_type': a.info.get('config_type'),
        'family': a.info['exploratory_source_family'],
        'natoms': len(a),
        'formula': a.get_chemical_formula(),
        'energy_eV_cell': float(a.info['REF_energy']),
        'fmax_eV_A': float(np.linalg.norm(f, axis=1).max()),
        'geometry_sha256': geometry_hash(a),
    })
manifest = {
    'status': 'EXPLORATORY_ACE_REFIT_DATASET_ONLY',
    'base_split': {'train_frames': len(train), 'development_frames': len(dev)},
    'supplement_candidates': len(new_frames),
    'supplement_added_after_exact_deduplication': len(accepted),
    'supplement_skipped': skipped,
    'derived_split': {'train_frames': len(train_out), 'development_frames': len(dev),
                      'test_frames': 0, 'blind_test_available': False,
                      'exact_train_dev_geometry_intersection_count': 0},
    'added_source_families': sorted({a.info['exploratory_source_family'] for a in accepted}),
    'validation_limit': 'The unchanged 7-frame development set is source-family correlated with some added COD 9009647 controls. Treat its metrics only as exploratory diagnostics, not independent generalization evidence.',
    'energy_convention': 'REF_energy / native GPAW extrapolated total energy; free energies remain separately in source metadata where available.',
    'force_convention': 'REF_forces / GPAW finite-smearing forces at sigma 0.10 eV.',
    'new_frames': new_summary,
    'source_records': source_records,
    'files_sha256': {},
}
for name in ('train.extxyz', 'dev.extxyz'):
    manifest['files_sha256'][name] = hashlib.sha256((OUT / name).read_bytes()).hexdigest()
(OUT / 'dataset_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
print(json.dumps({
    'status': manifest['status'],
    'base_train': len(train), 'added_candidates': len(new_frames),
    'added_unique': len(accepted), 'skipped_duplicates': skipped,
    'derived_train': len(train_out), 'development': len(dev),
    'exact_train_dev_intersections': len(intersections),
    'source_families': manifest['added_source_families'],
    'dev_is_independent': False,
}, indent=2))
