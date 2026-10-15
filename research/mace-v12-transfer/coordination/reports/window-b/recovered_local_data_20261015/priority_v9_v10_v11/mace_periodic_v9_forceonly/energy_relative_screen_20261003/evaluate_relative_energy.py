import csv
import json
from pathlib import Path

import numpy as np
import torch
from ase.io import read
from mace.calculators import MACECalculator

root = Path('/workspace/.local/lammps-ttm/run_1ps_diagnostic_20261001/potential_screening/periodic_interface_v4')
out = root / 'mace_periodic_v9_forceonly' / 'energy_relative_screen_20261003'
model_path = root / 'mace_periodic_v9_forceonly' / 'models' / 'Ag_Ti_Si_C_MACE_periodic_v9_forceonly.model'
frames = [
    ('AgTi_2p3871', root / 'pbe_external_v8_motif_holdouts' / 'AgTi' / 'retry1' / 'AgTi_cluster_PW_PBE.extxyz', 2.3870911918314524, 'force-structure was included in v9 training; its DFT energy was excluded'),
    ('AgTi_2p55', root / 'pbe_external_v9_holdouts' / 'AgTi_d2p55' / 'retry1' / 'AgTi_d2p55_cluster_PW_PBE.extxyz', 2.55, 'new separation within the same interface motif'),
    ('AgTi_2p7', root / 'pbe_external_v9_holdouts' / 'AgTi_d2p7' / 'retry1' / 'AgTi_d2p7_cluster_PW_PBE.extxyz', 2.6999999924868017, 'geometry represented in LCAO force data; PW-PBE method-transfer check'),
]
torch.set_default_dtype(torch.float64)
model = torch.load(model_path, map_location='cpu', weights_only=False).double().eval()
calc = MACECalculator(models=[model], device='cpu', default_dtype='float64')
results = []
reference = None
for label, path, pair_distance, note in frames:
    atoms = read(path, index=0, format='extxyz')
    symbols = atoms.get_chemical_symbols()
    from collections import Counter
    composition = dict(sorted(Counter(symbols).items()))
    dft_energy = float(atoms.get_potential_energy())
    atoms.calc = calc
    mace_energy = float(atoms.get_potential_energy())
    if reference is None:
        reference = (dft_energy, mace_energy)
    results.append({
        'label': label,
        'source': str(path),
        'n_atoms': len(atoms),
        'composition': composition,
        'formula': atoms.get_chemical_formula(),
        'cell_A': [float(x) for x in atoms.cell.lengths()],
        'pbc': [bool(x) for x in atoms.pbc],
        'central_pair_distance_A': pair_distance,
        'dft_energy_eV': dft_energy,
        'mace_v9_energy_eV': mace_energy,
        'absolute_error_eV': mace_energy - dft_energy,
        'note': note,
    })
base_dft, base_mace = next((r['dft_energy_eV'], r['mace_v9_energy_eV']) for r in results if r['label'] == 'AgTi_2p7')
for r in results:
    r['delta_dft_vs_2p7_eV'] = r['dft_energy_eV'] - base_dft
    r['delta_mace_vs_2p7_eV'] = r['mace_v9_energy_eV'] - base_mace
    r['relative_energy_error_vs_2p7_eV'] = r['delta_mace_vs_2p7_eV'] - r['delta_dft_vs_2p7_eV']

metadata = {
    'model': str(model_path),
    'model_sha256': __import__('hashlib').sha256(model_path.read_bytes()).hexdigest(),
    'reference': 'AgTi_2p7',
    'calibration_applied': False,
    'interpretation': 'Same-composition relative-energy diagnostic only; not a fully independent energy holdout and not an energy calibration. DFT setup is PW-PBE, 500 eV, Gamma-only, 28 A cubic cell for all three structures.',
    'records': results,
}
(out / 'relative_energy_results.json').write_text(json.dumps(metadata, indent=2) + '\n')
with (out / 'relative_energy_results.csv').open('w', newline='') as f:
    cols = list(results[0])
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader()
    for r in results:
        row = r.copy()
        row['composition'] = json.dumps(row['composition'], sort_keys=True)
        row['cell_A'] = json.dumps(row['cell_A'])
        row['pbc'] = json.dumps(row['pbc'])
        w.writerow(row)

lines = [
    '# v9 Ag–Ti existing-point relative-energy diagnostic',
    '',
    'This is an uncalibrated, post-hoc diagnostic on three existing periodic PW-PBE structures. No model parameters or energy zero were fitted. All structures have the same 49 atoms, composition C11Ag28SiTi9, 28 Å cubic cell, PBC, and DFT method (PBE, 500 eV, 1×1×1 k-points).',
    '',
    '| Separation (Å) | DFT E (eV) | v9 E (eV) | DFT ΔE from 2.70 Å (eV) | v9 ΔE (eV) | ΔE error (meV/cell) |',
    '|---:|---:|---:|---:|---:|---:|',
]
for r in results:
    lines.append(f"| {r['central_pair_distance_A']:.4f} | {r['dft_energy_eV']:.9f} | {r['mace_v9_energy_eV']:.9f} | {r['delta_dft_vs_2p7_eV']:.9f} | {r['delta_mace_vs_2p7_eV']:.9f} | {1000*r['relative_energy_error_vs_2p7_eV']:.3f} |")
lines += [
    '',
    'Caveats:',
    '- 2.3871 Å was included in v9 force training, while its energy label was excluded; this is not an independent structural holdout.',
    '- 2.70 Å is a PW-PBE method-transfer point on a geometry represented in LCAO force data.',
    '- 2.55 Å is a new distance in the same motif; one interval is insufficient to validate a whole energy curve.',
    '- These relative energies only probe the sampled Ag–Ti coordinate in one motif. They do not validate Ag–Si/Ag–C binding energies, thermal disorder, energy conservation under MD, or formation energies across compositions.',
    '- No calibration was applied. Element-reference calibration would require a training-only fit; validation energies must stay blind.',
    '',
    f"Model SHA256: `{metadata['model_sha256']}`",
    '',
]
(out / 'relative_energy_report.md').write_text('\n'.join(lines))
print(json.dumps(metadata, indent=2))
