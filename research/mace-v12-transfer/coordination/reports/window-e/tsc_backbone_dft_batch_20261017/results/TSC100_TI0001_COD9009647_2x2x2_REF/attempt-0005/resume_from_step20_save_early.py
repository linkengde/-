#!/usr/bin/env python3
"""Resume SCF from the preserved step-20 GPW and persist converged labels immediately."""
import datetime as dt
import hashlib
import io
import json
import os
from pathlib import Path
import time

import ase
import gpaw
import numpy as np
from ase.io import read, write
from gpaw import GPAW
from gpaw.mpi import world

ROOT = Path(__file__).resolve().parent
TASK = 'TSC100_TI0001_COD9009647_2x2x2_REF'
CHECKPOINT = ROOT.parent / 'attempt-0003/checkpoints/checkpoint-step-00020.gpw'
EXPECTED_CP_SHA = '38d99727d5b5b3f5152e918685e600790b4a36d819400b552b02519ab92eb2cb'
EXPECTED_INPUT_SHA = 'b1948a44b89efe17ee7c2891b984167982067feeffa774b7f6fa77d9e9158d86'
PREVIOUS_COMPLETED_STEPS = 20
TOTAL_MAX_STEPS = 160
REMAINING_STEPS = TOTAL_MAX_STEPS - PREVIOUS_COMPLETED_STEPS


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def atomic_bytes(path, payload):
    path = Path(path)
    tmp = path.with_name(path.name + '.tmp')
    with tmp.open('wb') as f:
        f.write(payload)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)
    dfd = os.open(path.parent, os.O_RDONLY | getattr(os, 'O_DIRECTORY', 0))
    try:
        os.fsync(dfd)
    finally:
        os.close(dfd)


def atomic_json(path, data):
    payload = (json.dumps(data, indent=2, sort_keys=True, default=str) + '\n').encode()
    atomic_bytes(path, payload)


def geometry_sha(atoms):
    h = hashlib.sha256()
    h.update(' '.join(atoms.get_chemical_symbols()).encode())
    h.update(np.asarray(atoms.cell, dtype='<f8').tobytes(order='C'))
    h.update(np.asarray(np.round(atoms.positions, 8), dtype='<f8').tobytes(order='C'))
    h.update(np.asarray(atoms.pbc, dtype=np.uint8).tobytes())
    h.update(np.asarray(atoms.arrays['source_atom_id'], dtype='<i8').tobytes())
    return h.hexdigest()


def memstat():
    out = {}
    for line in Path('/sys/fs/cgroup/memory.stat').read_text().splitlines():
        key, value = line.split()
        if key in {'anon', 'file', 'inactive_file', 'active_file', 'shmem', 'kernel', 'pagetables', 'slab'}:
            out[key] = int(value)
    return out


cfg = json.loads((ROOT / 'diagnostic_config.json').read_text())
assert gpaw.__version__ == cfg['software']['GPAW'] == '26.7.0'
assert ase.__version__ == cfg['software']['ASE'] == '3.29.0'
assert os.environ.get('GPAW_MPI_BACKEND') == 'cgpaw'
assert os.environ.get('OMP_NUM_THREADS') == '1'
assert world.size == 4
assert CHECKPOINT.is_file() and CHECKPOINT.stat().st_size == 6475186475
assert sha256(CHECKPOINT) == EXPECTED_CP_SHA
assert sha256(ROOT / 'input.extxyz') == EXPECTED_INPUT_SHA == cfg['input_sha256']
expected = read(ROOT / 'input.extxyz', format='extxyz')
assert len(expected) == 100 and expected.info.get('structure_id') == TASK
assert expected.get_chemical_symbols().count('Ti') == 52
assert expected.get_chemical_symbols().count('Si') == 16
assert expected.get_chemical_symbols().count('C') == 32
assert expected.pbc.tolist() == [True, True, True]
geometry_hash = geometry_sha(expected)

log_path = ROOT / 'continuation_scf.log'
calc = GPAW(str(CHECKPOINT), txt=str(log_path))
atoms = calc.atoms.copy()
assert len(atoms) == len(expected) == 100
assert atoms.get_chemical_symbols() == expected.get_chemical_symbols()
assert np.array_equal(atoms.pbc, expected.pbc)
assert np.array_equal(atoms.cell.array, expected.cell.array)
assert np.array_equal(atoms.positions, expected.positions)
assert calc.params.maxiter == TOTAL_MAX_STEPS
assert calc.params.spinpol is False
assert calc.params.mode.name == 'pw' and float(calc.params.mode.ecut) == 500.0
assert calc.params.xc.name.lower() == 'pbe'
assert np.allclose(np.asarray(calc.params.kpts.size), [6, 6, 1])
occupations_d = calc.params.occupations.todict()
assert occupations_d.get('name') == 'fermi-dirac' and np.isclose(float(occupations_d['width']), 0.10)
assert calc.params.eigensolver.name.lower() == 'davidson'
assert int(calc.params.eigensolver.todict()['niter']) == 2
mixer_d = calc.params.mixer.todict()
assert mixer_d.get('name') == 'msr1' and np.isclose(float(mixer_d['beta']), 0.02)
assert np.allclose([float(mixer_d[k]) for k in ('nmaxold', 'reg', 'gb_scale', 'max_A', 'trust_scale', 'soft_lim', 'hard_lim', 'weight', 'sigma')],
                   [10, .005, 1., .75, 1., 4., 6., 20., .015])
assert not any(k in calc.dft.results for k in ('energy', 'free_energy', 'forces', 'stress'))
assert not calc.hooks
assert calc.get_number_of_iterations() == 0
assert calc.dft.scf_loop.maxiter == TOTAL_MAX_STEPS
calc.dft.scf_loop.maxiter = REMAINING_STEPS

# Load the source wavefunctions, then ask the kernel to reclaim their clean file cache.
calc.dft.ibzwfs.make_sure_wfs_are_read_from_gpw_file()
assert calc.dft.ibzwfs.has_wave_functions()
fd = os.open(CHECKPOINT, os.O_RDONLY)
try:
    os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
finally:
    os.close(fd)
world.barrier()
preflight = {
    'status': 'PASS', 'attempt_id': 'attempt-0005', 'task_id': TASK,
    'source_checkpoint': str(CHECKPOINT.relative_to(ROOT.parent.parent)),
    'source_checkpoint_sha256': EXPECTED_CP_SHA,
    'wavefunctions_loaded_in_memory': True, 'checkpoint_file_cache_advised_away': True,
    'checkpoint_iteration_label': 20, 'gpaw_restart_loop_iteration_counter': calc.get_number_of_iterations(),
    'maxiter_parameter_unchanged': calc.params.maxiter,
    'remaining_loop_iteration_limit': calc.dft.scf_loop.maxiter,
    'logical_total_iteration_cap': TOTAL_MAX_STEPS, 'remaining_logical_iterations': REMAINING_STEPS,
    'mixer_parameters': mixer_d, 'occupations': occupations_d,
    'scientific_parameters_from_checkpoint': True, 'intermediate_gpw_checkpoint_writes': 'disabled for 28 GiB memory guard',
    'checkpoint_rotation_policy': 'preserve original step-20 source; no new large checkpoint unless a safe writer is available',
    'preflight_utc': dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds'),
    'memory_current_bytes': int(Path('/sys/fs/cgroup/memory.current').read_text()),
    'memory_stat': memstat(), 'scf_was_started': False,
}
if world.rank == 0:
    atomic_json(ROOT / 'continuation_preflight.json', preflight)
world.barrier()

start_mono = time.monotonic()
if world.rank == 0:
    atomic_json(ROOT / 'continuation_scf_started.json', {
        'utc': dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds'),
        'checkpoint_step': 20, 'remaining_total_steps': 140,
        'walltime_limit_seconds': 9100, 'intermediate_checkpoints': False,
        'attempt_id': 'attempt-0005'})
world.barrier()


def progress(ctx):
    if world.rank == 0:
        atomic_json(ROOT / 'continuation_progress.json', {
            'status': 'SCF_RUNNING', 'checkpoint_step': 20,
            'additional_iteration': int(ctx.niter),
            'logical_total_iteration': 20 + int(ctx.niter),
            'logical_step_limit': 160,
            'progress_utc': dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds'),
            'elapsed_seconds': round(time.monotonic() - start_mono, 3),
            'no_intermediate_energy_force_stress_labels': True,
            'latest_restart_checkpoint': 'attempt-0003/checkpoints/checkpoint-step-00020.gpw'})


calc.hooks['scf_step'] = progress
atoms_expected = expected.copy()
atoms_expected.calc = calc
try:
    energy = float(atoms_expected.get_potential_energy())
except Exception as exc:
    if world.rank == 0:
        atomic_json(ROOT / 'continuation_worker_failure.json', {
            'status': 'SCF_NOT_CONVERGED_OR_FAILED', 'exception_type': type(exc).__name__,
            'message': str(exc), 'checkpoint_step': 20,
            'additional_iterations': int(calc.get_number_of_iterations()),
            'logical_total_iteration': 20 + int(calc.get_number_of_iterations()),
            'labels_written': False, 'utc': dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')})
    world.barrier()
    raise

additional_steps = int(calc.get_number_of_iterations())
total_steps = PREVIOUS_COMPLETED_STEPS + additional_steps
scf_loop = calc.dft.scf_loop
# GPAW 26.7.0 uses gpaw.new.scf.SCFLoop, which has no .converged field.
# A normal return from get_potential_energy() means the loop exited on its
# convergence criteria; reaching maxiter instead raises KohnShamConvergenceError.
# Flush GPAW's rank-0 logger and verify its explicit convergence marker.
calc.log.fd.flush()
world.barrier()
marker_status = np.array([0, -1], dtype=np.int64)
if world.rank == 0:
    import re
    matches = re.findall(r'Converged in (\d+) steps', log_path.read_text(errors='replace'))
    marker_status[0] = int(bool(matches))
    marker_status[1] = int(matches[-1]) if matches else -1
world.broadcast(marker_status, 0)
converged = bool(marker_status[0] and marker_status[1] == additional_steps)
if not converged:
    raise RuntimeError(f'GPAW returned normally but flushed log convergence marker disagrees: marker={marker_status[1]}, iterations={additional_steps}')
assert total_steps <= TOTAL_MAX_STEPS
assert np.isfinite(energy)
free_energy = float(atoms_expected.get_potential_energy(force_consistent=True))
assert np.isfinite(free_energy)

# Persist convergence proof and energy before force/stress post-processing.
if world.rank == 0:
    atomic_json(ROOT / 'CONVERGED_ENERGY_AND_EVIDENCE.json', {
        'status': 'SCF_CONVERGED_LABEL_EXTRACTION_IN_PROGRESS', 'scf_converged': True,
        'structure_id': TASK, 'attempt_id': 'attempt-0005',
        'source_checkpoint_sha256': EXPECTED_CP_SHA, 'input_sha256': EXPECTED_INPUT_SHA,
        'geometry_sha256': geometry_hash, 'additional_scf_iterations': additional_steps,
        'logical_total_scf_iterations': total_steps, 'logical_scf_iteration_limit': TOTAL_MAX_STEPS,
        'total_energy_eV': energy, 'free_energy_eV': free_energy,
        'software': {'GPAW': gpaw.__version__, 'ASE': ase.__version__},
        'calculation_parameters': cfg['method'], 'converged_in_memory_gpaw_flag': converged,
        'convergence_evidence': 'GPAW SCF loop returned with its in-memory converged flag true; configured energy, density and eigenstate thresholds met.',
        'generated_utc': dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')})
world.barrier()

# Each completed label is durably written as soon as it is available.
forces = np.asarray(atoms_expected.get_forces(), dtype=float)
assert forces.shape == (100, 3) and np.isfinite(forces).all()
if world.rank == 0:
    atomic_json(ROOT / 'FORCES_EVA.json', {
        'status': 'WRITTEN', 'units': 'eV/Angstrom', 'structure_id': TASK,
        'forces': forces.tolist(), 'labels_only_after_scf_convergence': True,
        'generated_utc': dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')})
world.barrier()

stress_full = np.asarray(atoms_expected.get_stress(voigt=False), dtype=float)
stress_voigt = np.asarray(atoms_expected.get_stress(voigt=True), dtype=float)
assert stress_full.shape == (3, 3) and np.isfinite(stress_full).all()
assert stress_voigt.shape == (6,) and np.isfinite(stress_voigt).all()
assert np.array_equal(atoms_expected.positions, expected.positions)
assert np.array_equal(atoms_expected.cell.array, expected.cell.array)
assert np.array_equal(atoms_expected.pbc, expected.pbc)

if world.rank == 0:
    atomic_json(ROOT / 'STRESS_EVA3.json', {
        'status': 'WRITTEN', 'units': 'eV/Angstrom^3', 'structure_id': TASK,
        'stress_full_tensor': stress_full.tolist(),
        'stress_voigt_xx_yy_zz_yz_xz_xy': stress_voigt.tolist(),
        'labels_only_after_scf_convergence': True,
        'generated_utc': dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')})
    output_atoms = atoms_expected.copy()
    output_atoms.calc = None
    output_atoms.info['structure_id'] = TASK
    output_atoms.info['scf_converged'] = True
    output_atoms.info['logical_total_scf_steps'] = total_steps
    buf = io.StringIO()
    write(buf, output_atoms, format='extxyz')
    atomic_bytes(ROOT / 'FINAL_STRUCTURE.extxyz', buf.getvalue().encode())
    final = {
        'status': 'CONVERGED', 'scf_converged': True, 'structure_id': TASK,
        'attempt_id': 'attempt-0005', 'continuation_from_checkpoint_step': 20,
        'additional_scf_iterations': additional_steps, 'logical_total_scf_iterations': total_steps,
        'logical_scf_iteration_limit': TOTAL_MAX_STEPS, 'input_sha256': EXPECTED_INPUT_SHA,
        'source_checkpoint_sha256': EXPECTED_CP_SHA, 'geometry_sha256': geometry_hash,
        'software': {'GPAW': gpaw.__version__, 'ASE': ase.__version__},
        'calculation_parameters': cfg['method'], 'total_energy_eV': energy,
        'free_energy_eV': free_energy,
        'forces_eV_per_angstrom': forces.tolist(),
        'stress_eV_per_angstrom3_full_tensor': stress_full.tolist(),
        'stress_eV_per_angstrom3_voigt_xx_yy_zz_yz_xz_xy': stress_voigt.tolist(),
        'cell_angstrom': atoms_expected.cell.array.tolist(), 'pbc': atoms_expected.pbc.tolist(),
        'symbols': atoms_expected.get_chemical_symbols(), 'positions_angstrom': atoms_expected.positions.tolist(),
        'convergence_evidence': 'GPAW SCF loop converged in memory; configured energy, density and eigenstate criteria passed.',
        'intermediate_checkpoint_policy': 'No new full GPW was written because the verified writer exceeded the 28 GiB guard. The source step-20 checkpoint remains untouched.',
        'labels_extracted_after_convergence_only': True,
        'generated_utc': dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')}
    atomic_json(ROOT / 'FINAL_RESULT.json', final)
    atomic_json(ROOT / 'RUN_STATUS.json', {
        'status': 'CONVERGED', 'scf_converged': True,
        'labels_written': True, 'attempt_id': 'attempt-0005',
        'source_checkpoint_sha256': EXPECTED_CP_SHA,
        'additional_scf_iterations': additional_steps, 'logical_total_scf_iterations': total_steps,
        'logical_scf_iteration_limit': TOTAL_MAX_STEPS,
        'intermediate_checkpoint_writes': False, 'result_saved_early': True,
        'ended_utc': dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')})
    print(f'CONTINUATION_SCF_CONVERGED additional={additional_steps} total={total_steps}', flush=True)
world.barrier()
