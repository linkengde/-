#!/usr/bin/env python3
"""用派生的探索性数据集拟合 ACE 候选并评估训练/开发误差。"""
from __future__ import annotations
import hashlib, json, os, subprocess, time
from pathlib import Path
import numpy as np
import pandas as pd
from ase.io import read
from pyace import create_multispecies_basis_config, PyACECalculator
from pyace.linearacefit import LinearACEFit, LinearACEDataset
from sklearn.linear_model import Ridge

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'ace_models'
OUT.mkdir(exist_ok=True)
RCUT = 6.0

# 将 ASE 构型转换为 PyACE 需要的标签表；非周期方向扩胞时保持邻居距离不变。
def make_dataset(path: Path) -> pd.DataFrame:
    rows = []
    for atoms in read(path, index=':'):
        energy = float(atoms.info['REF_energy'])
        forces = np.asarray(atoms.arrays['REF_forces'], dtype=float).copy()
        a = atoms.copy()
        if a.pbc.any() and not a.pbc.all():
            old_dist = a.get_all_distances(mic=True)
            cell = np.asarray(a.cell).copy()
            scaled = a.get_scaled_positions(wrap=False)
            for axis in range(3):
                if not a.pbc[axis]:
                    length = float(np.linalg.norm(cell[axis]))
                    span = float(np.ptp(scaled[:, axis]) * length)
                    if length - span < RCUT + 2.0:
                        cell[axis] *= (span + RCUT + 2.0) / length
            a.set_cell(cell, scale_atoms=False)
            a.pbc = True
            new_dist = a.get_all_distances(mic=True)
            mask = (old_dist < RCUT + 1e-8) | (new_dist < RCUT + 1e-8)
            if not np.allclose(old_dist[mask], new_dist[mask], atol=1e-8):
                raise RuntimeError(f'{path}: 非周期扩胞改变了截断内距离')
        rows.append({'ase_atoms': a, 'energy': energy,
                     'energy_corrected': energy, 'forces': forces,
                     'NUMBER_OF_ATOMS': len(a)})
    return pd.DataFrame(rows)

# 固定与上一轮 ACE 基线相同的 6 Å 紧凑基组，便于比较。
potential_config = {
    'deltaSplineBins': 0.001,
    'elements': ['Ag', 'C', 'Si', 'Ti'],
    'embeddings': {'ALL': {'npot': 'FinnisSinclairShiftedScaled',
                           'fs_parameters': [1, 1], 'ndensity': 1}},
    'bonds': {'ALL': {'radbase': 'SBessel', 'radparameters': [6.25],
                      'rcut': RCUT, 'dcut': 0.01}},
    'functions': {'ALL': {'nradmax_by_orders': [8, 3, 1],
                          'lmax_by_orders': [0, 3, 1]}},
}
(OUT / 'basis_config.json').write_text(json.dumps(potential_config, indent=2) + '\n')
basis = create_multispecies_basis_config(potential_config=potential_config)
train = LinearACEDataset(basis, make_dataset(ROOT / 'train.extxyz'))
valid = LinearACEDataset(basis, make_dataset(ROOT / 'dev.extxyz'))
start = time.monotonic()
train.get_design_matrix(max_workers=2, verbose=True)
valid.get_design_matrix(max_workers=2, verbose=True)
y = train.get_target_vector()
n_energy_rows = len(train.df)

def evaluate(fit, frame_df):
    calc = PyACECalculator(fit.get_bbasis())
    de, force_errors, si_errors, per_config = [], [], [], []
    for _, row in frame_df.iterrows():
        atoms = row['ase_atoms'].copy()
        atoms.calc = calc
        e_pred = float(atoms.get_potential_energy())
        f_pred = atoms.get_forces()
        f_ref = np.asarray(row['forces'], dtype=float)
        df = f_pred - f_ref
        e_err = (e_pred - float(row['energy'])) / len(atoms)
        de.append(e_err)
        force_errors.append(df.reshape(-1, 3))
        site_error = np.linalg.norm(df, axis=1)
        symbols = np.asarray(atoms.get_chemical_symbols())
        if np.any(symbols == 'Si'):
            si_errors.extend(site_error[symbols == 'Si'])
        per_config.append({
            'config_type': atoms.info.get('config_type'),
            'natoms': len(atoms),
            'energy_error_meV_atom': 1000 * e_err,
            'force_component_rmse_eV_A': float(np.sqrt(np.mean(df**2))),
            'force_vector_rmse_eV_A': float(np.sqrt(np.mean(np.sum(df**2, axis=1)))),
            'force_vector_error_max_eV_A': float(site_error.max()),
        })
    force_errors = np.concatenate(force_errors, axis=0)
    return {
        'energy_rmse_meV_atom': float(1000 * np.sqrt(np.mean(np.square(de)))),
        'force_component_rmse_eV_A': float(np.sqrt(np.mean(np.square(force_errors)))),
        'force_vector_rmse_eV_A': float(np.sqrt(np.mean(np.sum(np.square(force_errors), axis=1)))),
        'force_vector_error_max_eV_A': float(max(x['force_vector_error_max_eV_A'] for x in per_config)),
        'si_site_vector_error_max_proxy_eV_A': float(max(si_errors)) if si_errors else None,
        'per_configuration': per_config,
    }

results = []
for alpha in [1e-4, 1e-3, 1e-2, 1e-1, 1.0]:
    # 显式使用能量权重10000、力权重1，与前一轮 ACE 基线保持一致。
    weights = np.ones(len(y), dtype=float)
    weights[:n_energy_rows] = 10000.0
    fit = LinearACEFit(train_dataset=train, model=Ridge(alpha=alpha))
    fit.fit(sample_weight=weights)
    stem = f'ACE_exploratory_B1616_E10000_F1_alpha{alpha:g}'
    yaml_path = OUT / f'{stem}.ace.yaml'
    fit.get_bbasis().save(str(yaml_path))
    yace_path = OUT / f'{stem}.yace'
    subprocess.run(['/workspace/.venvs/ace/bin/pace_yaml2yace',
                    str(yaml_path), '-o', str(yace_path)], check=True,
                   env={**os.environ, 'MPLCONFIGDIR': '/tmp/mplconfig-ace'})
    result = {
        'basis_functions': int(basis.total_number_of_functions),
        'cutoff_A': RCUT,
        'alpha': alpha,
        'weights': {'energy': 10000.0, 'force': 1.0},
        'train': evaluate(fit, train.df),
        'development': evaluate(fit, valid.df),
        'model_yaml': yaml_path.name,
        'model_yace': yace_path.name,
        'model_yaml_sha256': hashlib.sha256(yaml_path.read_bytes()).hexdigest(),
        'model_yace_sha256': hashlib.sha256(yace_path.read_bytes()).hexdigest(),
    }
    (OUT / f'{stem}.metrics.json').write_text(json.dumps(result, indent=2) + '\n')
    results.append(result)
    print(json.dumps({
        'alpha': alpha,
        'train': {k: result['train'][k] for k in ('energy_rmse_meV_atom', 'force_component_rmse_eV_A')},
        'development': {k: result['development'][k] for k in ('energy_rmse_meV_atom', 'force_component_rmse_eV_A', 'force_vector_error_max_eV_A')},
        'yaml': yaml_path.name,
        'yace': yace_path.name,
    }, separators=(',', ':')), flush=True)

summary = {
    'status': 'EXPLORATORY_FIT_COMPLETE_SCREENING_PENDING',
    'train_frames': len(train.df),
    'development_frames': len(valid.df),
    'development_independence': 'source-family correlated; diagnostic only',
    'wall_seconds': time.monotonic() - start,
    'results': results,
    'selected_candidate': min(results, key=lambda x: x['development']['force_component_rmse_eV_A'])['model_yace'],
}
(OUT / 'fit_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps({'selected_candidate': summary['selected_candidate'],
                  'wall_seconds': summary['wall_seconds']}, indent=2))
