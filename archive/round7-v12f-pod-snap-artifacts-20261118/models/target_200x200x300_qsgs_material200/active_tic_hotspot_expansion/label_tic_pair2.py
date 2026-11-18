#!/workspace/venvs/dft-labels/bin/python
"""GPAW PBE single-point forces for the independent second Ti-C crop."""
from __future__ import annotations
import json
import time
from pathlib import Path
import numpy as np
from ase.io import read, write
from gpaw import FermiDirac, GPAW, Mixer, PW
from gpaw.mpi import world

HERE = Path(__file__).resolve().parent
SOURCE = HERE/'candidates/tic_pair2_cold_10fs.extxyz'
OUTDIR = HERE/'labels'
STEM = 'tic_pair2_cold_10fs_pw350_gamma'
OUTDIR.mkdir(exist_ok=True)
atoms = read(SOURCE)
ids = np.asarray(atoms.arrays['source_atom_id'], dtype=int)
ci, pi = int(atoms.info['center_atom_id']), int(atoms.info['partner_atom_id'])
index = {int(atom_id): i for i, atom_id in enumerate(ids)}
if len(atoms) < 2 or ci not in index or pi not in index:
    raise RuntimeError('Unexpected crop size or center IDs')
log = OUTDIR/f'{STEM}.txt'
start = time.time()
atoms.calc = GPAW(mode=PW(350), xc='PBE', kpts=(1, 1, 1), symmetry='off',
    occupations=FermiDirac(0.10), mixer=Mixer(beta=0.05, nmaxold=8, weight=100),
    convergence={'density': 1e-5, 'eigenstates': 1e-5, 'energy': 1e-5}, maxiter=160,
    parallel={'domain': 1, 'band': max(1, world.size)}, txt=str(log))
energy = float(atoms.get_potential_energy())
free_energy = float(atoms.get_potential_energy(force_consistent=True))
forces = np.asarray(atoms.get_forces(), dtype=float)
iterations = int(atoms.calc.get_number_of_iterations())
if not np.isfinite(forces).all() or not np.isfinite(energy):
    raise RuntimeError('DFT produced non-finite results')
center_forces = {str(i): forces[index[i]].tolist() for i in (ci, pi)}
pod = np.asarray(atoms.arrays['POD_v12g_forces'], dtype=float)
errors = {str(i): float(np.linalg.norm(forces[index[i]]-pod[index[i]])) for i in (ci, pi)}
target = np.asarray(atoms.arrays['is_DFT_force_target'], dtype=int)
atoms.arrays['REF_forces'] = forces
atoms.info.update({'REF_energy': energy, 'free_energy_eV': free_energy,
    'DFT_method': 'GPAW PBE/PW350 Gamma, Fermi-Dirac sigma 0.10 eV',
    'DFT_center_forces_eV_A': json.dumps(center_forces),
    'DFT_POD_v12g_center_force_errors_eV_A': json.dumps(errors),
    'DFT_target_ids': json.dumps([ci, pi]),
    'label_status': 'GPAW single-point; target force rows are Ti-C center pair; edge rows excluded',
    'model_data_role': 'training_candidate; Si01 remains development holdout'})
atoms.calc = None
out = OUTDIR/f'{STEM}.extxyz'
write(out, atoms, format='extxyz')
logtext = log.read_text(errors='replace')
record = {'status': 'CONVERGED' if 'Converged in ' in logtext else 'SCF_CONVERGENCE_REVIEW_REQUIRED',
    'state': 'pair2_cold_10fs', 'center_atom_ids': [ci, pi], 'elements': ['Ti', 'C'],
    'pair_distance_A': float(atoms.info['center_partner_distance_A']), 'natoms': len(atoms),
    'energy_eV': energy, 'free_energy_eV': free_energy,
    'center_DFT_forces_eV_A': center_forces, 'center_POD_v12g_force_errors_eV_A': errors,
    'max_DFT_force_norm_eV_A': float(np.linalg.norm(forces, axis=1).max()),
    'DFT_target_rows': int(target.sum()), 'iterations': iterations, 'elapsed_s': time.time()-start,
    'output': str(out), 'log': str(log), 'method': 'GPAW PBE/PW350 Gamma, Fermi-Dirac sigma 0.10 eV'}
(OUTDIR/f'{STEM}.json').write_text(json.dumps(record, indent=2)+'\n')
if world.rank == 0:
    print(json.dumps(record), flush=True)
