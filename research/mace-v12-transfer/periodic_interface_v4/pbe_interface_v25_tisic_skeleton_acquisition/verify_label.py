#!/usr/bin/env python3
"""Verify one completed v25 skeleton DFT archive."""
from __future__ import annotations
import hashlib, json, sys
from pathlib import Path
import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parent
if len(sys.argv) != 2:
    raise SystemExit('usage: verify_label.py LABEL')
label = sys.argv[1]
manifest = json.loads((ROOT/'input_manifest.json').read_text())
job = next(r for r in manifest['preflight_jobs'] if r['label'] == label)
folder = ROOT/'calculations'/label
inventory = {line.split('  ',1)[1]:line.split('  ',1)[0] for line in (folder/'SHA256SUMS.txt').read_text().splitlines()}
for name, expected in inventory.items():
    actual = hashlib.sha256((folder/name).read_bytes()).hexdigest()
    if actual != expected: raise SystemExit(f'hash mismatch: {label}/{name}')
summary = json.loads((folder/'summary.json').read_text())
atoms = read(folder/f'{label}_PW_PBE.extxyz')
source = read(ROOT/job['input'])
checks = {
    'scf_converged': summary.get('scf_converged') is True,
    'registered_input_hash': summary.get('input_sha256') == job['input_sha256'],
    'composition': atoms.get_chemical_formula() == 'C4Si2Ti6' and set(atoms.get_chemical_symbols()) == {'Ti','Si','C'},
    'geometry_matches_input': len(atoms)==len(source) and np.allclose(atoms.positions,source.positions,rtol=0,atol=1e-10) and np.allclose(atoms.cell.array,source.cell.array,rtol=0,atol=1e-10),
    'finite_force_consistent_targets': np.isfinite(float(atoms.info['REF_energy'])) and np.isfinite(atoms.arrays['REF_forces']).all(),
    'finite_native_and_free_energy': np.isfinite(float(atoms.info['PW_PBE_energy_eV'])) and np.isfinite(float(atoms.info['PW_PBE_free_energy_eV'])),
    'finite_stress': all(np.isfinite(float(atoms.info[f'PW_PBE_stress_{x}_eV_A3'])) for x in ('xx','yy','zz','yz','xz','xy')),
    'no_checkpoint_in_archive': not any(p.suffix=='.gpw' for p in folder.iterdir()),
}
if not all(checks.values()): raise SystemExit(json.dumps({'status':'FAIL','checks':checks},indent=2))
print(json.dumps({'status':'PASS','label':label,'checks':checks,'archive_sha256':inventory},indent=2))
