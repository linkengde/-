#!/workspace/venvs/dft-labels/bin/python
"""Single-point GPAW PBE/PW350 Gamma label for a selected Si hotspot crop."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from ase.io import read, write
from gpaw import FermiDirac, GPAW, Mixer, PW
from gpaw.mpi import world

HERE = Path(__file__).resolve().parent
INDEX = int(sys.argv[1])
SOURCE = HERE / 'candidate_geometries.extxyz'
OUT = HERE / 'labels'
OUT.mkdir(exist_ok=True)
atoms = read(SOURCE, index=INDEX)
center_id = int(atoms.info['center_atom_id'])
ids = np.asarray(atoms.arrays['source_atom_id'], dtype=int)
center_matches = np.flatnonzero(ids == center_id)
if len(center_matches) != 1:
    raise RuntimeError(f'Center ID {center_id} does not map uniquely')
ci = int(center_matches[0])
start = time.time()
log = OUT / f'si_hotspot_{INDEX:02d}_pw350_gamma.txt'
atoms.calc = GPAW(mode=PW(350), xc='PBE', kpts=(1, 1, 1), symmetry='off',
                  occupations=FermiDirac(0.10),
                  mixer=Mixer(beta=0.05, nmaxold=8, weight=100),
                  convergence={'density': 1e-5, 'eigenstates': 1e-5, 'energy': 1e-5},
                  maxiter=160, parallel={'domain': max(1, world.size)}, txt=str(log))
energy = float(atoms.get_potential_energy())
free_energy = float(atoms.get_potential_energy(force_consistent=True))
forces = np.asarray(atoms.get_forces(), dtype=float)
iterations = int(atoms.calc.get_number_of_iterations())
if not np.isfinite(forces).all() or not np.isfinite(energy):
    raise RuntimeError('DFT returned non-finite energy or forces')
atoms.info.update({
    'REF_energy': energy,
    'REF_energy_unshifted_PBE_PW350_Gamma_eV': energy,
    'free_energy_eV': free_energy,
    'label_status': 'GPAW PBE/PW350 Gamma single-point label for a POD v11 Si hotspot',
    'DFT_method': 'GPAW PBE/PW350 Gamma, Fermi-Dirac sigma 0.10 eV',
    'DFT_cutoff_eV': 350,
    'DFT_kpts': [1, 1, 1],
    'DFT_center_force_eV_A': forces[ci].tolist(),
    'DFT_center_force_norm_eV_A': float(np.linalg.norm(forces[ci])),
    'DFT_minus_POD_v11_center_force_eV_A': (forces[ci] - np.asarray(atoms.info['center_POD_v11_force_eV_A'])).tolist(),
    'DFT_minus_POD_v11_center_force_error_eV_A': float(np.linalg.norm(forces[ci] - np.asarray(atoms.info['center_POD_v11_force_eV_A']))),
    'DFT_edge_force_note': 'Finite crop; only the center force is a DFT target.',
    'DFT_energy_note': 'Finite-crop energy has zero fit weight.',
    'model_data_role': str(atoms.info['training_split']),
})
atoms.arrays['REF_forces'] = forces
atoms.calc = None
stem = f'si_hotspot_{INDEX:02d}_pw350_gamma'
outxyz = OUT / f'{stem}.extxyz'
write(outxyz, atoms, format='extxyz')
record = {'index': INDEX, 'center_atom_id': center_id,
          'center_species': 'Si', 'natoms': len(atoms),
          'energy_eV': energy, 'free_energy_eV': free_energy,
          'center_force_DFT_eV_A': forces[ci].tolist(),
          'center_force_DFT_norm_eV_A': float(np.linalg.norm(forces[ci])),
          'center_force_POD_v11_eV_A': np.asarray(
              atoms.info['center_POD_v11_force_eV_A'], dtype=float).tolist(),
          'center_force_error_eV_A': atoms.info['DFT_minus_POD_v11_center_force_error_eV_A'],
          'max_force_norm_eV_A': float(np.linalg.norm(forces, axis=1).max()),
          'iterations': iterations, 'elapsed_s': time.time() - start,
          'role': str(atoms.info['training_split']),
          'output': str(outxyz), 'log': str(log), 'status': 'computed; audit pending'}
(OUT / f'{stem}.json').write_text(json.dumps(record, indent=2) + '\n')
if world.rank == 0:
    print(json.dumps(record), flush=True)
