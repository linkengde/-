#!/usr/bin/env python3
"""Run one registered periodic Ti3SiC2 PBE/PW label with durable checkpoints."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time

import numpy as np
from gpaw.mpi import world

ROOT = Path(__file__).resolve().parent
if len(sys.argv) != 3:
    raise SystemExit('usage: run_preflight_label.py LABEL SCRATCH_DIR')
LABEL, SCRATCH_ARG = sys.argv[1:]
manifest = json.loads((ROOT / 'input_manifest.json').read_text())
job = next(row for row in manifest['preflight_jobs'] if row['label'] == LABEL)
assert job['launch_enabled'] and world.size == 4
source = ROOT / job['input']
if hashlib.sha256(source.read_bytes()).hexdigest() != job['input_sha256']:
    raise SystemExit('registered input SHA256 changed')
scratch = Path(SCRATCH_ARG).resolve()
if world.rank == 0:
    scratch.mkdir(parents=True, exist_ok=False)
world.barrier()

import gpaw_data
os.environ['GPAW_SETUP_PATH'] = str(gpaw_data.datapath())
from ase.io import read, write
from gpaw import GPAW, PW, restart
from gpaw.mixer import Mixer
from gpaw.occupations import FermiDirac

atoms = read(source)
assert len(atoms) == 12 and atoms.get_chemical_formula() == 'C4Si2Ti6'
assert atoms.pbc.tolist() == [True, True, True]
assert not any(k in atoms.info for k in ('energy','REF_energy','PW_PBE_energy_eV'))
assert not any(k in atoms.arrays for k in ('forces','REF_forces','PW_PBE_forces'))
method = manifest['method']
started = time.time()
progress_path = scratch / 'progress.json'
checkpoint_meta_path = ROOT / 'calculations/latest_checkpoint/latest_checkpoint.json'
checkpoint_dir = checkpoint_meta_path.parent
checkpoint_dir.mkdir(parents=True, exist_ok=True)
restart_path = os.environ.get('V25_RESTART_CHECKPOINT')
iteration_offset = int(os.environ.get('V25_RESTART_ITERATION_OFFSET', '0'))
if restart_path:
    restarted_atoms, calc = restart(restart_path, txt=str(scratch / 'gpaw.log'), parallel={'sl_auto': True})
    if (restarted_atoms.get_chemical_symbols() != atoms.get_chemical_symbols()
            or restarted_atoms.pbc.tolist() != atoms.pbc.tolist()
            or not np.array_equal(restarted_atoms.cell.array, atoms.cell.array)
            or not np.allclose(restarted_atoms.positions, atoms.positions, rtol=0, atol=1e-10)):
        raise SystemExit('checkpoint geometry does not match this registered input')
else:
    calc = GPAW(
        mode=PW(method['cutoff_eV']), xc=method['xc'], kpts=tuple(method['kpts']),
        occupations=FermiDirac(method['smearing_eV']),
        mixer=Mixer(method['mixer']['beta'], method['mixer']['nmaxold'], method['mixer']['weight']),
        parallel={'sl_auto': True}, convergence=method['convergence'], maxiter=method['maxiter'],
        txt=str(scratch / 'gpaw.log'))
atoms.calc = calc
if world.rank == 0:
    progress_path.write_text(json.dumps({
        'label': LABEL, 'status': 'running', 'iteration': iteration_offset,
        'source_sha256': job['input_sha256'], 'mpi_ranks': world.size,
        'started_epoch_s': started, 'resumed_from_checkpoint': bool(restart_path),
        'checkpoint_path': restart_path, 'iteration_offset': iteration_offset,
    }, indent=2) + '\n')
world.barrier()

pause_path = ROOT / 'pause.request'
checkpoint_interval = int(method['checkpoint_interval_scf_steps'])

def scf_step(context):
    iteration = iteration_offset + int(context.niter)
    if world.rank == 0 and iteration % 5 == 0:
        progress_path.write_text(json.dumps({
            'label': LABEL, 'status': 'running', 'iteration': iteration,
            'elapsed_s': time.time() - started, 'source_sha256': job['input_sha256'],
            'mpi_ranks': world.size, 'resumed_from_checkpoint': bool(restart_path),
        }, indent=2) + '\n')
    if iteration > 0 and iteration % checkpoint_interval == 0:
        tmp = checkpoint_dir / 'latest_state.gpw.tmp'
        stable = checkpoint_dir / 'latest_state.gpw'
        if world.rank == 0:
            tmp.unlink(missing_ok=True)
        world.barrier()
        calc.write(str(tmp), mode='all')
        world.barrier()
        if world.rank == 0:
            digest = hashlib.sha256(tmp.read_bytes()).hexdigest()
            os.replace(tmp, stable)
            meta = {
                'status': 'stable_checkpoint_snapshot', 'label': LABEL,
                'input_sha256': job['input_sha256'], 'checkpoint_sha256': digest,
                'checkpoint_bytes': stable.stat().st_size,
                'last_reported_scf_iteration': iteration,
                'created_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
            }
            meta_tmp = checkpoint_meta_path.with_suffix('.json.tmp')
            meta_tmp.write_text(json.dumps(meta, indent=2) + '\n')
            meta_tmp.replace(checkpoint_meta_path)
    pause_flag = np.array([int(pause_path.exists()) if world.rank == 0 else 0], dtype=np.int32)
    world.broadcast(pause_flag, 0)
    pause_requested = bool(pause_flag[0])
    if pause_requested and iteration > 0 and iteration % checkpoint_interval == 0:
        if world.rank == 0:
            progress_path.write_text(json.dumps({
                'label': LABEL, 'status': 'paused_checkpointed', 'iteration': iteration,
                'source_sha256': job['input_sha256'], 'mpi_ranks': world.size,
                'checkpoint_sha256': json.loads(checkpoint_meta_path.read_text())['checkpoint_sha256'],
            }, indent=2) + '\n')
        world.barrier()
        raise RuntimeError('user pause requested after durable SCF checkpoint')

calc.hooks['scf_step'] = scf_step
converged = {'value': False}
calc.hooks['converged'] = lambda: converged.__setitem__('value', True)
native_energy = float(atoms.get_potential_energy(force_consistent=False))
free_energy = float(calc.results['free_energy'])
forces = np.asarray(atoms.get_forces(), dtype=float)
stress = np.asarray(atoms.get_stress(voigt=True), dtype=float)
iterations = iteration_offset + int(calc.get_number_of_iterations())
if not converged['value']:
    raise SystemExit(f'SCF failed to converge within {iterations} iterations')
if not np.isfinite(native_energy) or not np.isfinite(free_energy) or not np.isfinite(forces).all() or not np.isfinite(stress).all():
    raise SystemExit('nonfinite DFT result')

if world.rank == 0:
    out = atoms.copy()
    out.calc = None
    out.info.update({
        'config_type': LABEL, 'PW_PBE_energy_eV': native_energy,
        'PW_PBE_free_energy_eV': free_energy, 'REF_energy': free_energy,
        'energy_target_convention': 'force-consistent finite-width free energy at Fermi-Dirac sigma=0.10 eV',
        'dataset_role': job['role'], 'dataset_split': job['split'],
        'source_doi': '10.1007/BF00899949',
        'source_family': 'Jeitschko_Nowotny_1967_single_crystal',
    })
    for name, value in zip(('xx','yy','zz','yz','xz','xy'), stress):
        out.info[f'PW_PBE_stress_{name}_eV_A3'] = float(value)
    out.arrays['PW_PBE_forces'] = forces.copy()
    out.arrays['REF_forces'] = forces.copy()
    write(scratch / f'{LABEL}_PW_PBE.extxyz', out, format='extxyz')
    summary = {
        'label': LABEL, 'input': job['input'], 'input_sha256': job['input_sha256'],
        'source_family': 'Jeitschko_Nowotny_1967_single_crystal',
        'source_doi': '10.1007/BF00899949', 'source_cif_sha256': manifest['source']['source_cif_sha256'],
        'atom_count': len(atoms), 'formula': atoms.get_chemical_formula(),
        'cell_angstrom': atoms.cell.array.tolist(), 'pbc': atoms.pbc.tolist(),
        'native_extrapolated_energy_eV_cell': native_energy,
        'free_energy_eV_cell': free_energy,
        'training_energy_target_eV_cell': free_energy,
        'energy_target_convention': 'force-consistent finite-width free energy at sigma 0.10 eV',
        'forces_eV_A': forces.tolist(), 'stress_voigt_eV_A3': stress.tolist(),
        'scf_converged': True, 'scf_iterations': iterations,
        'resumed_from_checkpoint': bool(restart_path), 'checkpoint_iteration_offset': iteration_offset,
        'max_force_eV_A': float(np.linalg.norm(forces, axis=1).max()),
        'rms_force_vector_eV_A': float(np.sqrt(np.mean(np.sum(forces**2, axis=1)))),
        'method': method, 'dataset_role': job['role'], 'dataset_split': job['split'],
        'training_integration': 'not performed; candidate acquisition only',
        'elapsed_s': time.time() - started,
    }
    (scratch / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    progress_path.write_text(json.dumps({
        'label': LABEL, 'status': 'complete', 'iteration': iterations,
        'elapsed_s': summary['elapsed_s'], 'scf_converged': True,
        'source_sha256': job['input_sha256'], 'mpi_ranks': world.size,
    }, indent=2) + '\n')
    print(json.dumps(summary, indent=2), flush=True)
world.barrier()
