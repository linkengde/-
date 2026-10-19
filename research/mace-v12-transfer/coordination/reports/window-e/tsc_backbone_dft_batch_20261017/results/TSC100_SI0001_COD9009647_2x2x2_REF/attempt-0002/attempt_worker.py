#!/usr/bin/env python3
"""One fixed-geometry Si-terminated reference SCF attempt; no automatic retry."""
import argparse
import datetime as dt
import hashlib
import json
import os
import re
import shutil
import sys
import time
from pathlib import Path

import ase
import gpaw
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write
from gpaw import FermiDirac, GPAW, PW
from gpaw.dft import MSR1
from gpaw.mpi import world

ROOT = Path(__file__).resolve().parent
TASK = 'TSC100_SI0001_COD9009647_2x2x2_REF'
EXPECTED_INPUT_SHA = 'ac8b804516424215694c284d3bcc36d23da590f80390f47a0805c87cedb764c4'
GUARD = 28 * 1024**3
EXPECTED_CHECKPOINT_BYTES = 3_239_353_332
CHECKPOINT_HEADROOM = 4 * 1024**3
SAVE_EVERY = 10
last_checkpoint = None
last_checkpoint_step = None
next_checkpoint_due = 1
started = time.monotonic()


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')


def sha256(path):
    h = hashlib.sha256()
    fd = os.open(path, os.O_RDONLY)
    try:
        offset = 0
        while True:
            block = os.read(fd, 1024 * 1024)
            if not block:
                break
            h.update(block)
            os.posix_fadvise(fd, offset, len(block), os.POSIX_FADV_DONTNEED)
            offset += len(block)
    finally:
        os.close(fd)
    return h.hexdigest()


def atomic_json(path, value):
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    os.replace(temp, path)


def geometry_sha(atoms):
    h = hashlib.sha256()
    h.update(' '.join(atoms.get_chemical_symbols()).encode())
    h.update(np.asarray(atoms.cell, dtype='<f8').tobytes(order='C'))
    h.update(np.asarray(np.round(atoms.positions, 8), dtype='<f8').tobytes(order='C'))
    h.update(np.asarray(atoms.pbc, dtype=np.uint8).tobytes())
    h.update(np.asarray(atoms.arrays['source_atom_id'], dtype='<i8').tobytes())
    return h.hexdigest()


def memory_current():
    return int(Path('/sys/fs/cgroup/memory.current').read_text())


class SafeBudgetStop(Exception):
    pass


parser = argparse.ArgumentParser()
parser.add_argument('--config', required=True)
parser.add_argument('--config-sha256', required=True)
parser.add_argument('--execution-lock', required=True)
parser.add_argument('--execution-lock-sha256', required=True)
args = parser.parse_args()

config_path = Path(args.config).resolve()
lock_path = Path(args.execution_lock).resolve()
if sha256(config_path) != args.config_sha256:
    raise RuntimeError('retry config SHA mismatch')
if sha256(lock_path) != args.execution_lock_sha256:
    raise RuntimeError('retry execution-lock SHA mismatch')
cfg = json.loads(config_path.read_text())
lock = json.loads(lock_path.read_text())
assert cfg['structure_id'] == lock['task_id'] == TASK
assert cfg['input_sha256'] == lock['input_sha256'] == EXPECTED_INPUT_SHA
assert gpaw.__version__ == cfg['software']['GPAW'] == '26.7.0'
assert ase.__version__ == cfg['software']['ASE'] == '3.29.0'
assert os.environ.get('GPAW_MPI_BACKEND') == 'cgpaw'
assert os.environ.get('OMP_NUM_THREADS') == '1'
assert world.size == 4
assert os.environ.get('ATTEMPT_ID') == 'attempt-0002'

input_path = ROOT / 'input.extxyz'
assert sha256(input_path) == EXPECTED_INPUT_SHA
atoms = read(input_path, format='extxyz')
symbols = atoms.get_chemical_symbols()
assert len(atoms) == 100
assert {s: symbols.count(s) for s in sorted(set(symbols))} == {'C': 32, 'Si': 20, 'Ti': 48}
assert atoms.info.get('surface_termination') == 'Si'
assert atoms.pbc.tolist() == [True, True, True]
assert np.array_equal(atoms.arrays.get('source_atom_id'), np.arange(1, 101))
assert not any(k in atoms.info for k in ('energy', 'free_energy', 'stress'))
assert not any(k in atoms.arrays for k in ('forces', 'stress'))
geometry_before = geometry_sha(atoms)
assert geometry_before == cfg['geometry_sha256']

method = cfg['method']
mixer = method['mixer']
calc = GPAW(
    mode=PW(method['cutoff_eV']),
    xc=method['xc'],
    kpts=tuple(method['kpoints_size']),
    occupations=FermiDirac(method['fermi_dirac_width_eV']),
    spinpol=method['spin_polarized'],
    convergence=method['scf_convergence'],
    maxiter=method['maxiter'],
    mixer=MSR1(
        nmaxold=mixer['nmaxold'], beta=mixer['beta'], reg=mixer['reg'],
        gb_scale=mixer['gb_scale'], max_A=mixer['max_A'],
        trust_scale=mixer['trust_scale'], soft_lim=mixer['soft_lim'],
        hard_lim=mixer['hard_lim'], weight=mixer['weight'], sigma=mixer['sigma']),
    parallel={'sl_auto': True},
    txt=str(ROOT / 'gpaw.log'))
atoms.calc = calc
checkpoint_dir = ROOT / 'checkpoints'
checkpoint_dir.mkdir(exist_ok=True)
temp_checkpoint = checkpoint_dir / 'latest.partial.gpw'
final_checkpoint = checkpoint_dir / 'latest.gpw'
stop_request = ROOT / 'STOP_REQUEST'


def write_checkpoint(step, forced=False):
    global last_checkpoint, last_checkpoint_step, next_checkpoint_due
    decision = np.array([1, 0, 0], dtype=np.int64)
    if world.rank == 0:
        free = shutil.disk_usage(ROOT).free
        expected = EXPECTED_CHECKPOINT_BYTES
        disk_required = expected * (2 if final_checkpoint.exists() else 1) + 6 * 1024**3
        mem = memory_current()
        if free < disk_required:
            decision[0] = -1
            decision[1] = free
            decision[2] = disk_required
        elif temp_checkpoint.exists():
            decision[0] = -2
        elif mem > GUARD - CHECKPOINT_HEADROOM:
            decision[0] = 0
            decision[1] = mem
            decision[2] = GUARD - CHECKPOINT_HEADROOM
    world.broadcast(decision, 0)
    if decision[0] == -1:
        raise RuntimeError(f'checkpoint disk guard: {decision[1]} free < {decision[2]} required')
    if decision[0] == -2:
        raise RuntimeError('partial checkpoint path exists; refusing overwrite')
    if decision[0] == 0:
        if world.rank == 0:
            atomic_json(ROOT / 'CHECKPOINT_DEFERRED.json', {
                'status': 'DEFERRED_LOW_MEMORY_HEADROOM', 'step': int(step),
                'memory_current_bytes': int(decision[1]), 'required_headroom_bytes': CHECKPOINT_HEADROOM,
                'guard_bytes': GUARD, 'latest_checkpoint_step': last_checkpoint_step,
                'utc': utc()})
            print(f'CHECKPOINT_DEFERRED step={step} memory={int(decision[1])}', flush=True)
            next_checkpoint_due = step + 1
        due_step = np.array([next_checkpoint_due if world.rank == 0 else 0], dtype=np.int64)
        world.broadcast(due_step, 0)
        next_checkpoint_due = int(due_step[0])
        return False
    t0 = time.monotonic()
    calc.write(str(temp_checkpoint), mode='all', precision='single')
    world.barrier()
    outcome = np.array([0, -1], dtype=np.int64)
    failed = 0
    message = ''
    if world.rank == 0:
        try:
            if not temp_checkpoint.is_file() or temp_checkpoint.stat().st_size == 0:
                raise RuntimeError('GPAW returned without a nonempty checkpoint')
            fd = os.open(temp_checkpoint, os.O_RDWR)
            try:
                os.fsync(fd)
            finally:
                os.close(fd)
            digest = sha256(temp_checkpoint)
            size = temp_checkpoint.stat().st_size
            os.replace(temp_checkpoint, final_checkpoint)
            last_checkpoint = final_checkpoint
            last_checkpoint_step = int(step)
            next_checkpoint_due = int(step) + SAVE_EVERY
            atomic_json(ROOT / 'CHECKPOINT_STATUS.json', {
                'status': 'VALID_CHECKPOINT_WRITTEN', 'step': int(step),
                'path': 'checkpoints/latest.gpw', 'bytes': size,
                'sha256': digest, 'precision': 'single', 'mode': 'all',
                'wavefunctions_included': True,
                'write_elapsed_seconds': round(time.monotonic() - t0, 3),
                'older_checkpoint_replaced_only_after_new_checkpoint_hashed': True,
                'utc': utc()})
            outcome[:] = [1, int(step)]
            print(f'TARGET_CHECKPOINT step={step} bytes={size} sha256={digest}', flush=True)
        except Exception as exc:
            failed = 1
            message = f'{type(exc).__name__}: {exc}'
            try:
                temp_checkpoint.unlink()
            except FileNotFoundError:
                pass
    failed = int(world.sum(failed))
    if failed:
        raise RuntimeError(message or 'rank zero checkpoint finalization failed')
    world.broadcast(outcome, 0)
    if outcome[0] != 1:
        raise RuntimeError('checkpoint finalization status was not broadcast')
    last_checkpoint_step = int(outcome[1])
    last_checkpoint = final_checkpoint
    next_checkpoint_due = last_checkpoint_step + SAVE_EVERY
    world.barrier()
    return True


def scf_step(ctx):
    step = int(ctx.niter)
    if world.rank == 0:
        calc.log.fd.flush()
        atomic_json(ROOT / 'SCF_PROGRESS.json', {
            'status': 'SCF_RUNNING', 'iteration': step, 'iteration_limit': method['maxiter'],
            'elapsed_seconds': round(time.monotonic() - started, 3),
            'latest_checkpoint_step': last_checkpoint_step,
            'no_intermediate_energy_force_stress_labels': True, 'utc': utc()})
    world.barrier()
    if stop_request.exists():
        if last_checkpoint_step != step:
            write_checkpoint(step, forced=True)
        raise SafeBudgetStop('walltime shutdown reserve reached; latest checkpoint preserved')
    if step >= next_checkpoint_due:
        write_checkpoint(step)


calc.hooks['scf_step'] = scf_step
try:
    energy = float(atoms.get_potential_energy(force_consistent=False))
except SafeBudgetStop as exc:
    if world.rank == 0:
        atomic_json(ROOT / 'SCF_STOPPED.json', {
            'status': 'STOPPED_BEFORE_CONVERGENCE', 'reason': str(exc),
            'iterations_completed': int(calc.get_number_of_iterations()),
            'latest_checkpoint_step': last_checkpoint_step,
            'labels_written': False, 'utc': utc()})
    world.barrier()
    sys.exit(75)
except Exception as exc:
    if world.rank == 0:
        atomic_json(ROOT / 'SCF_FAILURE.json', {
            'status': 'SCF_FAILED_OR_NOT_CONVERGED', 'exception_type': type(exc).__name__,
            'message': str(exc), 'iterations_completed': int(calc.get_number_of_iterations()),
            'latest_checkpoint_step': last_checkpoint_step,
            'labels_written': False, 'utc': utc()})
    world.barrier()
    raise

calc.log.fd.flush()
world.barrier()
markers = []
if world.rank == 0:
    log_text = (ROOT / 'gpaw.log').read_text(errors='replace')
    markers = re.findall(r'Converged in (\d+) steps', log_text)
marker_ok = np.array([bool(markers), int(markers[-1]) if markers else -1], dtype=np.int64)
world.broadcast(marker_ok, 0)
steps = int(marker_ok[1])
assert marker_ok[0] and steps == int(calc.get_number_of_iterations())
assert geometry_sha(atoms) == geometry_before
assert np.isfinite(energy)

# Extract official labels only after GPAW convergence is proven.
free_energy = float(atoms.get_potential_energy(force_consistent=True))
forces = np.asarray(atoms.get_forces(), dtype=float)
stress = np.asarray(atoms.get_stress(voigt=False), dtype=float)
assert forces.shape == (100, 3) and np.isfinite(forces).all()
assert stress.shape == (3, 3) and np.isfinite(stress).all()
assert np.isfinite(free_energy)
if world.rank == 0:
    labelled = atoms.copy()
    labelled.calc = SinglePointCalculator(labelled, energy=energy,
                                          free_energy=free_energy,
                                          forces=forces, stress=stress)
    tmp = ROOT / 'labelled.extxyz.tmp'
    write(tmp, labelled, format='extxyz')
    os.replace(tmp, ROOT / 'labelled.extxyz')
    cfg_sha = hashlib.sha256(config_path.read_bytes()).hexdigest()
    lock_sha = hashlib.sha256(lock_path.read_bytes()).hexdigest()
    result = {
        'structure_id': TASK, 'attempt_id': 'attempt-0002', 'status': 'SUCCESS',
        'scf_converged': True, 'scf_iterations': steps,
        'energy_eV': energy, 'force_consistent_free_energy_eV': free_energy,
        'forces_eV_per_angstrom': forces.tolist(), 'stress_eV_per_angstrom3': stress.tolist(),
        'max_force_eV_per_angstrom': float(np.linalg.norm(forces, axis=1).max()),
        'atom_count': len(atoms), 'composition': {'Ti': 48, 'Si': 20, 'C': 32},
        'input_sha256': EXPECTED_INPUT_SHA, 'geometry_sha256': geometry_before,
        'config_sha256': cfg_sha, 'execution_lock_sha256': lock_sha,
        'data_role': 'TRAIN_FORCE_LABEL_CANDIDATE',
        'labels_only_after_scf_convergence': True,
    }
    atomic_json(ROOT / 'SCF_CONVERGENCE.json', {
        'status': 'SCF_CONVERGED_LABELS_WRITTEN', 'scf_converged': True,
        'scf_iterations': steps, 'convergence': method['scf_convergence'],
        'gpaw_log_sha256': sha256(ROOT / 'gpaw.log'), 'utc': utc()})
    atomic_json(ROOT / 'ENERGIES.json', {
        'status': 'WRITTEN', 'energy_eV': energy,
        'force_consistent_free_energy_eV': free_energy,
        'labels_only_after_scf_convergence': True})
    atomic_json(ROOT / 'FORCES.json', {
        'status': 'WRITTEN', 'units': 'eV/Angstrom',
        'forces': forces.tolist(), 'labels_only_after_scf_convergence': True})
    atomic_json(ROOT / 'STRESS.json', {
        'status': 'WRITTEN', 'units': 'eV/Angstrom^3',
        'stress': stress.tolist(), 'labels_only_after_scf_convergence': True})
    atomic_json(ROOT / 'RESULT.json', result)
    print(f'TARGET_SCF_CONVERGED_LABELS_WRITTEN steps={steps}', flush=True)
world.barrier()
