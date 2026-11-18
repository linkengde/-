#!/workspace/venvs/dft-labels/bin/python
"""Audit the four new targeted PBE/PW350 Gamma hotspot labels."""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import numpy as np
from ase.io import read

HERE = Path(__file__).resolve().parent
SOURCE = HERE / 'candidate_geometries.extxyz'
LABELS = HERE / 'labels'


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def audit_one(idx: int, source) -> dict:
    stem = f'si_hotspot_{idx:02d}_pw350_gamma'
    xyz, meta, log = (LABELS / f'{stem}.extxyz', LABELS / f'{stem}.json', LABELS / f'{stem}.txt')
    missing = [str(p) for p in (xyz, meta, log) if not p.is_file() or p.stat().st_size == 0]
    if missing:
        return {'index': idx, 'center_atom_id': int(source.info['center_atom_id']),
                'status': 'PENDING', 'missing': missing}
    labeled = read(xyz)
    record = json.loads(meta.read_text())
    issues = []
    for what, a, b in (('element order', source.numbers, labeled.numbers),
                       ('coordinates', source.positions, labeled.positions),
                       ('cell', source.cell.array, labeled.cell.array),
                       ('PBC', source.pbc, labeled.pbc),
                       ('source IDs', source.arrays['source_atom_id'], labeled.arrays.get('source_atom_id', []))):
        if not np.array_equal(a, b):
            issues.append(f'{what} changed')
    force = np.asarray(labeled.arrays.get('REF_forces', []), dtype=float)
    if force.shape != (len(labeled), 3) or not np.isfinite(force).all():
        issues.append('DFT force array is missing, malformed, or non-finite')
    center_id = int(source.info['center_atom_id'])
    center_indices = np.flatnonzero(np.asarray(labeled.arrays.get('source_atom_id', []), dtype=int) == center_id)
    if len(center_indices) != 1:
        issues.append('center atom ID does not map uniquely')
        center_force = np.full(3, np.nan)
    elif force.shape != (len(labeled), 3):
        center_force = np.full(3, np.nan)
    else:
        center_force = force[int(center_indices[0])]
    rows = [line for line in log.read_text(errors='replace').splitlines() if '|iter:' in line]
    markers = re.findall(r'[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?c(?=\s*\|)', rows[-1]) if rows else []
    converged = len(markers) >= 3
    if not converged:
        issues.append('last GPAW SCF row lacks three convergence markers')
    if int(record.get('center_atom_id', -1)) != center_id:
        issues.append('JSON center ID disagrees with source')
    if int(record.get('natoms', -1)) != len(labeled):
        issues.append('JSON atom count disagrees')
    dft_norm = float(np.linalg.norm(center_force))
    pod_force = np.asarray(source.info['center_POD_v11_force_eV_A'], dtype=float)
    error = float(np.linalg.norm(center_force - pod_force))
    if abs(float(record.get('center_force_DFT_norm_eV_A', np.nan)) - dft_norm) > 1e-9:
        issues.append('JSON center-force norm disagrees with extxyz')
    return {'index': idx, 'center_atom_id': center_id,
            'role': str(source.info['training_split']),
            'status': 'PASS' if not issues else 'FAIL', 'issues': issues,
            'natoms': len(labeled), 'SCF_iterations': len(rows), 'SCF_converged': converged,
            'POD_v11_center_force_eV_A': pod_force.tolist(),
            'DFT_center_force_eV_A': center_force.tolist(),
            'POD_v11_force_norm_eV_A': float(np.linalg.norm(pod_force)),
            'DFT_force_norm_eV_A': dft_norm,
            'DFT_minus_POD_v11_force_error_eV_A': error,
            'sha256': {'extxyz': sha256(xyz), 'json': sha256(meta), 'log': sha256(log)}}


def main() -> None:
    sources = read(SOURCE, index=':')
    frames = [audit_one(i, atom) for i, atom in enumerate(sources)]
    complete = all(f['status'] == 'PASS' for f in frames)
    report = {'status': 'PASS' if complete else 'INCOMPLETE',
              'method': 'GPAW PBE/PW350 Gamma; Fermi-Dirac sigma 0.10 eV',
              'expected_labels': len(sources),
              'passed_labels': sum(f['status'] == 'PASS' for f in frames),
              'role': 'targeted development data; center force only; ID 104359 reserved as v12 development holdout',
              'frames': frames}
    out = HERE / 'hotspot_label_audit.json'
    out.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    raise SystemExit(0 if complete else 2)


if __name__ == '__main__':
    main()
