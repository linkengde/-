import csv
import json
import hashlib
from collections import Counter
from pathlib import Path

import numpy as np
import torch
from ase.data import chemical_symbols
from ase.io import read
from mace.calculators import MACECalculator

base = Path('/workspace/.local/lammps-ttm/run_1ps_diagnostic_20261001/potential_screening/periodic_interface_v4')
model_path = base / 'mace_periodic_v9_forceonly/models/Ag_Ti_Si_C_MACE_periodic_v9_forceonly.model'
data_dir = base / 'mace_periodic_v9_forceonly/data'
out = base / 'mace_periodic_v9_forceonly/energy_relative_screen_20261003'
torch.set_default_dtype(torch.float64)
model = torch.load(model_path, map_location='cpu', weights_only=False).double().eval()
calc = MACECalculator(models=[model], device='cpu', default_dtype='float64')
elements = [chemical_symbols[int(z)] for z in model.atomic_numbers.tolist()]
results = {}
for split in ('train', 'valid', 'test'):
    structures = read(data_dir / f'{split}.extxyz', index=':', format='extxyz')
    rows = []
    for idx, atoms in enumerate(structures):
        if 'REF_energy' not in atoms.info:
            continue
        y = float(atoms.info['REF_energy'])
        atoms.calc = calc
        pred = float(atoms.get_potential_energy())
        counts = Counter(atoms.get_chemical_symbols())
        rows.append({
            'split': split,
            'frame': idx,
            'formula': atoms.get_chemical_formula(),
            'n_atoms': len(atoms),
            'composition': dict(sorted(counts.items())),
            'reference_energy_eV': y,
            'mace_energy_eV': pred,
            'residual_DFT_minus_MACE_eV': y-pred,
        })
    results[split] = rows

train = results['train']
if not train:
    raise RuntimeError('No energy-labeled training frames found')
A = np.array([[r['composition'].get(el, 0) for el in elements] for r in train], dtype=float)
y = np.array([r['residual_DFT_minus_MACE_eV'] for r in train], dtype=float)
mu, _, rank, svals = np.linalg.lstsq(A, y, rcond=None)
for split, rows in results.items():
    for r in rows:
        n = np.array([r['composition'].get(el, 0) for el in elements], dtype=float)
        r['training_only_element_reference_correction_eV'] = float(n @ mu)
        r['corrected_residual_eV'] = r['residual_DFT_minus_MACE_eV'] - float(n @ mu)

summary = {
    'model': str(model_path),
    'model_sha256': hashlib.sha256(model_path.read_bytes()).hexdigest(),
    'elements': elements,
    'train_energy_frame_count': len(train),
    'composition_matrix_rank': int(rank),
    'composition_matrix_columns': len(elements),
    'fitted_element_reference_corrections_eV_per_atom': {el: float(v) for el, v in zip(elements, mu)},
    'singular_values': [float(v) for v in svals],
    'train_raw_residual_rmse_meV_atom': float(1000*np.sqrt(np.mean([r['residual_DFT_minus_MACE_eV']**2/r['n_atoms']**2 for r in train]))),
    'train_raw_energy_residual_rmse_meV_atom': float(1000*np.sqrt(np.mean([(r['residual_DFT_minus_MACE_eV']/r['n_atoms'])**2 for r in train]))),
    'interpretation': 'Least-squares species reference terms fit exclusively on training energy frames. This is a diagnostic only; no model file was modified. Rank deficiency means individual elemental shifts are not identifiable from the training compositions.',
    'records': results,
}
# Save a tidy table for each energy-labeled split.
with (out / 'training_reference_fit_metrics.csv').open('w', newline='') as f:
    fields = ['split','frame','formula','n_atoms','composition','reference_energy_eV','mace_energy_eV','residual_DFT_minus_MACE_eV','training_only_element_reference_correction_eV','corrected_residual_eV']
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    for split in ('train','valid','test'):
        for r in results[split]:
            q = r.copy()
            q['composition'] = json.dumps(q['composition'], sort_keys=True)
            writer.writerow({k:q.get(k) for k in fields})
(out / 'training_reference_fit.json').write_text(json.dumps(summary, indent=2) + '\n')
lines = [
    '# v9 training-only energy reference diagnostic',
    '',
    'Only the v9 training energy labels were used to fit reference-energy correction coefficients. Validation/test and external PW-PBE frames were not used in the fit. No model artifact was changed.',
    '',
    f"- Energy-labeled training frames: {len(train)}",
    f"- Species composition matrix rank: {rank}/{len(elements)}",
    f"- Species terms (eV/atom): {', '.join(f'{el}={v:+.8f}' for el,v in zip(elements,mu))}",
    f"- Training energy RMSE after the species-reference fit: {summary['train_raw_energy_residual_rmse_meV_atom']:.3f} meV/atom",
    '',
    '| Split | Energy frames | Raw RMSE (meV/atom) | After train-only ref fit (meV/atom) |',
    '|---|---:|---:|---:|',
]
for split in ('train','valid','test'):
    rs=results[split]
    if rs:
        raw=np.array([r['residual_DFT_minus_MACE_eV']/r['n_atoms'] for r in rs])
        corr=np.array([r['corrected_residual_eV']/r['n_atoms'] for r in rs])
        lines.append(f"| {split} | {len(rs)} | {1000*np.sqrt(np.mean(raw**2)):.3f} | {1000*np.sqrt(np.mean(corr**2)):.3f} |")
    else:
        lines.append(f'| {split} | 0 | n/a | n/a |')
lines += ['', 'This diagnostic tests whether the observed absolute-energy bias is mostly a reference term. If the design-matrix rank is below four, individual element references cannot be uniquely inferred from the current training set; add energy-labeled structures with independent compositions before using per-element corrections. A small corrected RMSE on train alone is not validation; valid/test and untouched external structures control that claim.', '']
(out / 'training_reference_fit_report.md').write_text('\n'.join(lines))
print(json.dumps({k:v for k,v in summary.items() if k!='records'}, indent=2))
for split in ('train','valid','test'):
    for r in results[split]:
        print(split, r['frame'], r['formula'], r['n_atoms'], f"raw={r['residual_DFT_minus_MACE_eV']/r['n_atoms']*1000:+.2f} meV/atom", f"corr={r['corrected_residual_eV']/r['n_atoms']*1000:+.2f} meV/atom")
