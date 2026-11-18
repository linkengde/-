#!/workspace/venvs/dft-labels/bin/python
"""Validate and checksum the eight PW350/Gamma local force labels."""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import numpy as np
from ase.io import read

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = ROOT / 'target_active_geometries.extxyz'
EXPECTED_CUTOFF = 350


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    source_frames = read(SOURCE, index=':', format='extxyz')
    rows = []
    issues = []
    if len(source_frames) != 8:
        issues.append(f'Expected 8 source crops, found {len(source_frames)}')
    for idx in range(8):
        stem = f'target_local_{idx:02d}_pw350_gamma'
        extxyz = HERE / f'{stem}.extxyz'
        record_path = HERE / f'{stem}.json'
        log_path = HERE / f'{stem}.txt'
        if not all(p.is_file() for p in (extxyz, record_path, log_path)):
            issues.append(f'Frame {idx}: one or more label files are missing')
            continue
        labeled = read(extxyz, format='extxyz')
        record = json.loads(record_path.read_text())
        source = source_frames[idx]
        log = log_path.read_text(errors='replace')
        match = re.search(r'Converged in (\d+) steps', log)
        scf_converged = match is not None
        if not scf_converged:
            issues.append(f'Frame {idx}: no explicit SCF convergence marker')
        if len(source) != len(labeled):
            issues.append(f'Frame {idx}: source/output atom-count mismatch')
        if source.get_chemical_formula() != labeled.get_chemical_formula():
            issues.append(f'Frame {idx}: source/output formula mismatch')
        pos_delta = float(np.max(np.abs(source.positions - labeled.positions)))
        cell_delta = float(np.max(np.abs(np.asarray(source.cell) - np.asarray(labeled.cell))))
        if pos_delta > 1e-10 or cell_delta > 1e-10 or not np.array_equal(source.pbc, labeled.pbc):
            issues.append(f'Frame {idx}: source geometry/cell/PBC changed during labelling')
        required_arrays = ('REF_forces', 'MTP_predicted_forces_eV_A', 'DFT_force_error_eV_A')
        missing_arrays = [name for name in required_arrays if name not in labeled.arrays]
        if missing_arrays:
            issues.append(f'Frame {idx}: missing force arrays {missing_arrays}')
            dft_forces = np.asarray(labeled.arrays.get('REF_forces', np.empty((0, 3))))
        else:
            dft_forces = np.asarray(labeled.arrays['REF_forces'], dtype=float)
            if dft_forces.shape != (len(labeled), 3) or not np.isfinite(dft_forces).all():
                issues.append(f'Frame {idx}: malformed/nonfinite force targets')
        energy = float(labeled.info.get('REF_energy', np.nan))
        free_energy = float(labeled.info.get('free_energy_eV', np.nan))
        if not (np.isfinite(energy) and np.isfinite(free_energy)):
            issues.append(f'Frame {idx}: energy/free-energy target is nonfinite')
        if int(labeled.info.get('DFT_cutoff_eV', -1)) != EXPECTED_CUTOFF:
            issues.append(f'Frame {idx}: cutoff metadata mismatch')
        if int(record['iterations']) != int(match.group(1)) if match else True:
            issues.append(f'Frame {idx}: SCF iteration count mismatch')
        rows.append({
            'frame_index': idx,
            'source_atom_id': int(record['source_atom_id']),
            'formula': labeled.get_chemical_formula(),
            'natoms': len(labeled),
            'scf_converged': scf_converged,
            'scf_iterations': int(match.group(1)) if match else None,
            'max_dft_force_norm_eV_A': float(np.linalg.norm(dft_forces, axis=1).max()) if len(dft_forces) else None,
            'center_dft_force_norm_eV_A': float(record['center_DFT_force_norm_eV_A']),
            'center_original_mtp_force_norm_eV_A': float(record['center_MTP_force_norm_eV_A']),
            'energy_eV': energy,
            'free_energy_eV': free_energy,
            'energy_minus_free_energy_eV': energy - free_energy,
            'source_position_max_delta_A': pos_delta,
            'source_cell_max_delta_A': cell_delta,
            'pbc_match': bool(np.array_equal(source.pbc, labeled.pbc)),
            'label_files_sha256': {p.name: sha256(p) for p in (extxyz, record_path, log_path)},
        })

    result = {
        'status': 'PASS' if not issues and len(rows) == 8 else 'FAIL',
        'method': 'GPAW PBE/PW350 Gamma, Fermi-Dirac sigma 0.10 eV',
        'labels_completed': len(rows),
        'source_geometry_unchanged': all(r['source_position_max_delta_A'] <= 1e-10 and r['source_cell_max_delta_A'] <= 1e-10 and r['pbc_match'] for r in rows),
        'all_scf_converged': all(r['scf_converged'] for r in rows),
        'all_force_targets_finite': not any('malformed/nonfinite force targets' in x for x in issues),
        'issues': issues,
        'frames': rows,
        'fit_compatibility': {
            'direct_energy_merge_into_v10': 'NO',
            'reason': 'These are finite local crops at PW350/Gamma with finite-crop boundary effects and explicitly raw mixed-cutoff energy references; v10 uses full periodic training configurations and a different validated energy convention.',
            'force_use': 'Potentially useful for a separately weighted center-force fit after confirming the center has a complete 5 A neighborhood; edge forces must be excluded/downweighted.'
        },
    }
    (HERE / 'local_dft_label_manifest.json').write_text(json.dumps(result, indent=2) + '\n')
    md = [
        '# Round7 局部 DFT 标签验收', '',
        f"状态：**{result['status']}**（{len(rows)}/8）", '',
        '8 个 PW350/Γ 单点任务均检查 SCF 收敛标记、输入/输出坐标与晶胞、PBC、元素组成、力阵列有限性，并计算 extxyz、JSON 和 GPAW 日志的 SHA-256。', '',
        '|帧|来源原子|式量|SCF 步数|中心 DFT 力 (eV/Å)|原 MTP 力 (eV/Å)|E−F (eV)|',
        '|---:|---:|---|---:|---:|---:|---:|'
    ]
    for r in rows:
        md.append(f"|{r['frame_index']}|{r['source_atom_id']}|{r['formula']}|{r['scf_iterations']}|{r['center_dft_force_norm_eV_A']:.4f}|{r['center_original_mtp_force_norm_eV_A']:.3f}|{r['energy_minus_free_energy_eV']:.6f}|")
    md += ['', '## 拟合兼容性', '',
           '这些局部裁剪标签不直接并入 v10。它们是 PW350/Γ 有限裁剪的能量标签，不能与全周期 v10 能量不加处理地混合；裁剪边缘原子的力还受边界影响。若要用于新拟合，应优先只用中心区域的力标签，并在新版本中单独记录能量基准和原子权重。', '',
           '逐文件 SHA-256 和完整逐帧信息见 `local_dft_label_manifest.json`。']
    if issues:
        md += ['', '## 发现的问题', ''] + [f'- {issue}' for issue in issues]
    (HERE / 'local_dft_label_audit.md').write_text('\n'.join(md) + '\n')
    print(json.dumps({k: result[k] for k in ('status', 'labels_completed', 'source_geometry_unchanged', 'all_scf_converged', 'all_force_targets_finite', 'issues')}, indent=2))


if __name__ == '__main__':
    main()
