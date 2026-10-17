#!/usr/bin/env python3
"""Run one frozen Ti3SiC2 diagnostic single point. Called under MPI."""
import argparse
import hashlib
import json
import os
import re
import time
from pathlib import Path

import ase
import gpaw as gpaw_module
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write
from gpaw import FermiDirac, GPAW, PW
from gpaw.dft import MSR1
from gpaw.mpi import world


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def parse_scf_iterations(log_text):
    matches = re.findall(r'Converged in (\d+) steps', log_text)
    if not matches:
        raise RuntimeError('GPAW log has no converged-SCF marker; refusing to write labels')
    return int(matches[-1])


parser = argparse.ArgumentParser()
parser.add_argument('--task-id', required=True)
parser.add_argument('--input', required=True)
parser.add_argument('--run-dir', required=True)
parser.add_argument('--input-sha256', required=True)
parser.add_argument('--execution-lock', required=True)
parser.add_argument('--execution-lock-sha256', required=True)
args = parser.parse_args()

input_path = Path(args.input).resolve()
run_dir = Path(args.run_dir).resolve()
lock_path = Path(args.execution_lock).resolve()
run_dir.mkdir(parents=True, exist_ok=True)
assert sha256(input_path) == args.input_sha256, 'input SHA mismatch'
assert sha256(lock_path) == args.execution_lock_sha256, 'execution lock SHA mismatch'
lock = json.loads(lock_path.read_text())
assert lock.get('execution_authorized') is True
assert lock.get('role', '').startswith('DIAGNOSTIC_ONLY')
method = lock['scientific_method']
scf = method['scf_convergence']
mixer_spec = scf['mixer']
assert mixer_spec['type'] == 'MSR1' and method['eigensolver']['type'] == 'GPAW default Davidson'
assert gpaw_module.__version__ == method['version']
assert ase.__version__ == method['ASE_version']
assert int(os.environ.get('OMP_NUM_THREADS', '0')) == method['parallel']['OMP_NUM_THREADS']
assert os.environ.get('GPAW_MPI_BACKEND') == 'cgpaw', 'GPAW 26.7 MPI backend must be explicitly set to cgpaw'
assert os.environ.get('GPAW_SETUP_PATH') and Path(os.environ['GPAW_SETUP_PATH']).is_dir()
atoms = read(input_path)
symbols = atoms.get_chemical_symbols()
assert len(atoms) == 48
assert symbols.count('Ti') == 24 and symbols.count('Si') == 8 and symbols.count('C') == 16
assert atoms.pbc.tolist() == [True, True, True]
assert 'energy' not in atoms.info and 'free_energy' not in atoms.info
assert 'forces' not in atoms.arrays

t0 = time.monotonic()
calc = GPAW(
    mode=PW(method['cutoff_eV']),
    xc=method['xc'],
    kpts=tuple(method['kpoints']['size']),  # tuple uses GPAW's Monkhorst-Pack convention.
    occupations=FermiDirac(method['occupations_width_eV']),
    convergence={k: scf[k] for k in ('energy', 'density', 'eigenstates')},
    mixer=MSR1(nmaxold=mixer_spec['nmaxold'], beta=mixer_spec['beta'],
              reg=mixer_spec['reg'], gb_scale=mixer_spec['gb_scale'],
              max_A=mixer_spec['max_A'], trust_scale=mixer_spec['trust_scale'],
              soft_lim=mixer_spec['soft_lim'], hard_lim=mixer_spec['hard_lim'],
              weight=mixer_spec['weight'], sigma=mixer_spec['sigma']),
    maxiter=scf['maxiter'],
    parallel={'sl_auto': method['parallel']['sl_auto']},
    txt=str(run_dir / 'gpaw.log'))
atoms.calc = calc

energy = float(atoms.get_potential_energy(force_consistent=False))
free_energy = float(atoms.get_potential_energy(force_consistent=True))
forces = np.asarray(atoms.get_forces(), dtype=float)
stress = np.asarray(atoms.get_stress(voigt=False), dtype=float)
log_text = (run_dir / 'gpaw.log').read_text(errors='replace')
iterations = parse_scf_iterations(log_text)
reported_iterations = int(calc.get_number_of_iterations())
if reported_iterations != iterations:
    raise RuntimeError(f'SCF iteration mismatch: log={iterations}, calculator={reported_iterations}')
# GPAW 26.7's new ASECalculator raises on an unconverged SCF and reports a
# "Converged in N steps" marker only after all SCF criteria and post-SCF hooks pass.
converged = True
elapsed = time.monotonic() - t0

if not converged:
    raise RuntimeError(f'SCF did not converge after {iterations} iterations; refusing to write labels')
if forces.shape != (48, 3) or stress.shape != (3, 3):
    raise RuntimeError(f'Unexpected force/stress shapes: {forces.shape}, {stress.shape}')
if not (np.isfinite(energy) and np.isfinite(free_energy) and np.all(np.isfinite(forces)) and np.all(np.isfinite(stress))):
    raise RuntimeError('Non-finite result; refusing to write labels')

if world.rank == 0:
    labelled = atoms.copy()
    labelled.calc = SinglePointCalculator(labelled, energy=energy,
                                         free_energy=free_energy,
                                         forces=forces, stress=stress)
    output_path = run_dir / 'labelled.extxyz'
    write(output_path, labelled, format='extxyz')
    result = {
        'task_id': args.task_id,
        'status': 'SUCCESS',
        'scf_converged': converged,
        'scf_iterations': iterations,
        'energy_eV': energy,
        'free_energy_eV': free_energy,
        'forces_eV_per_A': forces.tolist(),
        'stress_eV_per_A3': stress.tolist(),
        'stress_GPa': (stress * 160.21766208).tolist(),
        'atom_count': len(atoms),
        'composition': {'Ti': 24, 'Si': 8, 'C': 16},
        'input_sha256': args.input_sha256,
        'parameter_lock_sha256': args.execution_lock_sha256,
        'kpoint_type': method['kpoints']['type'],
        'kpoints': method['kpoints']['size'],
        'eigensolver': method['eigensolver'],
        'mixer': mixer_spec,
        'GPAW_version': gpaw_module.__version__,
        'ASE_version': ase.__version__,
        'output_extxyz_sha256': sha256(output_path),
        'elapsed_seconds': elapsed,
        'units': {'energy': 'eV', 'forces': 'eV/Angstrom', 'stress': 'eV/Angstrom^3 and GPa'},
        'role': 'DIAGNOSTIC_ONLY_NOT_TRAIN_DEV_TEST'
    }
    (run_dir / 'RESULT.json').write_text(json.dumps(result, indent=2) + '\n')
