#!/usr/bin/env python3
"""One-shot, guarded GPAW retry for the Si-terminated reference structure."""
import argparse
import datetime as dt
import hashlib
import importlib.metadata
import json
import os
import shutil
import signal
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

import ase
import gpaw
import gpaw_data
import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parent
BASE = ROOT.parents[2]
REPO = ROOT.parents[7]
TASK = 'TSC100_SI0001_COD9009647_2x2x2_REF'
BRANCH = 'cloud-e-b1-dft'
LIMIT = 18000
STOP_RESERVE = 600
GUARD = 28 * 1024**3
MIN_START_FREE = 16 * 1024**3
MIN_RUNTIME_FREE = 6 * 1024**3
EXPECTED_INPUT_SHA = 'ac8b804516424215694c284d3bcc36d23da590f80390f47a0805c87cedb764c4'
VENV_PYTHON = Path('/workspace/.venvs/gpaw-mpi/bin/python')
GPAW_EXE = Path('/workspace/.venvs/gpaw-mpi/bin/gpaw')
MPI_EXE = shutil.which('mpiexec') or '/workspace/.venvs/gpaw-mpi/bin/mpiexec'
CONFIG = ROOT / 'diagnostic_config.json'
LOCK = ROOT / 'EXECUTION_LOCK.json'
WORKER = ROOT / 'attempt_worker.py'
PROFILE = ROOT / 'RESOURCE_PROFILE.jsonl'


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def atomic_json(path, data):
    path = Path(path)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(data, indent=2, sort_keys=True) + '\n')
    os.replace(temp, path)


def read_int(name):
    value = Path('/sys/fs/cgroup', name).read_text().strip()
    return 2**63 - 1 if value == 'max' else int(value)


def memory_stats():
    out = {}
    for line in Path('/sys/fs/cgroup/memory.stat').read_text().splitlines():
        key, value = line.split()
        out[key] = int(value)
    return {key: out.get(key, 0) for key in
            ('anon', 'file', 'inactive_file', 'active_file', 'shmem', 'kernel', 'slab')}


def geom_sha(atoms):
    h = hashlib.sha256()
    h.update(' '.join(atoms.get_chemical_symbols()).encode())
    h.update(np.asarray(atoms.cell, dtype='<f8').tobytes(order='C'))
    h.update(np.asarray(np.round(atoms.positions, 8), dtype='<f8').tobytes(order='C'))
    h.update(np.asarray(atoms.pbc, dtype=np.uint8).tobytes())
    h.update(np.asarray(atoms.arrays['source_atom_id'], dtype='<i8').tobytes())
    return h.hexdigest()


def active_compute_processes():
    found = []
    for proc in Path('/proc').iterdir():
        if not proc.name.isdigit():
            continue
        try:
            comm = (proc / 'comm').read_text().strip().lower()
            state = next(x.split()[1] for x in (proc / 'status').read_text().splitlines()
                         if x.startswith('State:'))
            if state != 'Z' and comm in {'gpaw', 'mpiexec', 'mpirun', 'mpiexec.hydra', 'orted', 'prte'}:
                found.append({'pid': int(proc.name), 'comm': comm, 'state': state})
        except (OSError, StopIteration):
            pass
    return found


def preflight(expected_lock_sha):
    cfg = json.loads(CONFIG.read_text())
    lock = json.loads(LOCK.read_text())
    input_path = BASE / 'inputs' / f'{TASK}.extxyz'
    failed_path = BASE / 'results' / TASK
    checks = {}
    checks['execution_authorized_by_current_user'] = lock.get('execution_authorized') is True
    checks['manual_attempt_id'] = lock.get('attempt_id') == 'attempt-0002'
    checks['automatic_retry_disabled'] = lock.get('automatic_retry') is False
    checks['lock_sha256'] = sha256(LOCK) == expected_lock_sha
    checks['config_sha256'] = sha256(CONFIG) == lock.get('config_sha256')
    checks['worker_sha256'] = sha256(WORKER) == lock.get('worker_sha256')
    checks['supervisor_sha256'] = sha256(Path(__file__)) == lock.get('supervisor_sha256')
    checks['input_sha256'] = sha256(input_path) == EXPECTED_INPUT_SHA == lock.get('input_sha256')
    atoms = read(input_path, format='extxyz')
    comp = dict(Counter(atoms.get_chemical_symbols()))
    checks['geometry_identity'] = (
        len(atoms) == 100 and comp == {'Ti': 48, 'Si': 20, 'C': 32}
        and atoms.pbc.tolist() == [True, True, True]
        and np.array_equal(atoms.arrays.get('source_atom_id'), np.arange(1, 101))
        and geom_sha(atoms) == cfg.get('geometry_sha256')
        and atoms.info.get('structure_id') == TASK
        and atoms.info.get('surface_termination') == 'Si'
        and 'energy' not in atoms.info and 'free_energy' not in atoms.info
        and 'forces' not in atoms.arrays and 'stress' not in atoms.info)
    old_status = json.loads((failed_path / 'RUN_STATUS.json').read_text())
    checks['failed_attempt_preserved_no_checkpoint'] = (
        old_status.get('status') == 'FAILED'
        and old_status.get('stop_reason') == 'WALLTIME_LIMIT_18000_SECONDS'
        and old_status.get('checkpoint_written') is False
        and not (failed_path / 'RESULT.json').exists())
    old_sha = subprocess.run(['sha256sum', '-c', 'SHA256SUMS.txt'], cwd=failed_path,
                             text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    checks['failed_attempt_sha256'] = old_sha.returncode == 0
    package_sha = subprocess.run(['sha256sum', '-c', 'SHA256SUMS.txt'], cwd=BASE,
                                 text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    checks['input_package_sha256'] = package_sha.returncode == 0
    checks['GPAW_version'] = gpaw.__version__ == cfg['software']['GPAW'] == '26.7.0'
    checks['ASE_version'] = ase.__version__ == cfg['software']['ASE'] == '3.29.0'
    checks['gpaw_data_version'] = importlib.metadata.version('gpaw-data') == cfg['software']['gpaw-data'] == '1.2.1'
    setup_path = Path(gpaw_data.datapath()).resolve()
    checks['PAW_setups'] = setup_path.is_dir() and all(
        sha256(setup_path / name) == digest for name, digest in cfg['paw_sha256'].items())
    cpu = Path('/sys/fs/cgroup/cpu.max').read_text().split()
    checks['cpu_quota_at_least_4'] = len(cpu) == 2 and int(cpu[0]) / int(cpu[1]) >= 4
    memory_limit = read_int('memory.max')
    memory_current = read_int('memory.current')
    checks['memory_limit_32GiB'] = memory_limit >= 32 * 1024**3
    checks['memory_below_28GiB_guard'] = memory_current < GUARD
    free = shutil.disk_usage('/workspace').free
    checks['workspace_start_free_at_least_16GiB'] = free >= MIN_START_FREE
    checks['no_active_compute_process'] = not active_compute_processes()
    checks['attempt2_output_not_started'] = not (ROOT / 'RUN_MANIFEST.json').exists() and not (ROOT / 'RESULT.json').exists()
    cp = json.loads((ROOT / 'checkpoint_preflight' / 'checkpoint_preflight.json').read_text())
    checks['checkpoint_write_read_preflight'] = (
        cp.get('status') == 'PASS' and cp.get('wavefunctions_reloaded') is True
        and cp.get('checkpoint_mode') == 'all' and cp.get('checkpoint_precision') == 'single'
        and cp['write']['peak_memory_current_bytes'] < GUARD
        and cp['read']['peak_memory_current_bytes'] < GUARD
        and cp.get('test_checkpoint_bytes', 2**63) < 4 * 1024**3)
    checks['base_input_commit_ancestor'] = subprocess.run(
        ['git', 'merge-base', '--is-ancestor', lock['input_package_commit'], 'HEAD'], cwd=REPO).returncode == 0
    subprocess.run(['git', 'fetch', '--no-tags', 'origin', lock['C_review_commit']], cwd=REPO, check=True)
    review_path = 'research/mace-v12-transfer/coordination/reports/window-c/tsc_surface_input_review_20261017/C_INPUT_REVIEW.md'
    decision_path = 'research/mace-v12-transfer/coordination/reports/window-c/tsc_surface_input_review_20261017/C_INPUT_DECISION.json'
    review_bytes = subprocess.check_output(['git', 'show', f"{lock['C_review_commit']}:{review_path}"], cwd=REPO)
    decision_bytes = subprocess.check_output(['git', 'show', f"{lock['C_review_commit']}:{decision_path}"], cwd=REPO)
    checks['C_review_artifact_hashes'] = (
        hashlib.sha256(review_bytes).hexdigest() == lock['C_review_sha256']
        and hashlib.sha256(decision_bytes).hexdigest() == lock['C_decision_sha256'])
    subprocess.run(['git', 'fetch', '--no-tags', 'origin',
                    f'+refs/heads/{BRANCH}:refs/remotes/origin/{BRANCH}'], cwd=REPO, check=True)
    local_head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip()
    remote_head = subprocess.check_output(['git', 'rev-parse', f'refs/remotes/origin/{BRANCH}'], cwd=REPO, text=True).strip()
    checks['branch_fresh_and_tracked_tree_clean'] = (
        local_head == remote_head
        and not subprocess.check_output(['git', 'status', '--porcelain', '--untracked-files=no'], cwd=REPO, text=True).strip()
        and subprocess.run(['git', 'diff', '--cached', '--quiet'], cwd=REPO).returncode == 0)
    details = {
        'task_id': TASK, 'attempt_id': 'attempt-0002', 'branch': BRANCH,
        'local_head': local_head, 'remote_head': remote_head,
        'workspace_free_bytes': free, 'memory_current_bytes': memory_current,
        'memory_limit_bytes': memory_limit, 'memory_stat': memory_stats(),
        'cpu_max': ' '.join(cpu), 'active_compute_processes': active_compute_processes(),
        'GPAW': gpaw.__version__, 'ASE': ase.__version__,
        'gpaw_data': importlib.metadata.version('gpaw-data'),
        'checkpoint_preflight_status': cp.get('status'),
        'checkpoint_preflight_peak_bytes': max(cp['write']['peak_memory_current_bytes'], cp['read']['peak_memory_current_bytes']),
        'checkpoint_preflight_wavefunctions_reloaded': cp.get('wavefunctions_reloaded'),
        'C_review_artifact_hashes_verified': checks['C_review_artifact_hashes'],
        'old_attempt_sha_check_output': old_sha.stdout,
        'input_package_sha_check_output': package_sha.stdout,
    }
    return checks, details


def sample(elapsed):
    row = {'utc': utc(), 'elapsed_seconds': round(elapsed, 3),
           'memory_current_bytes': read_int('memory.current'),
           'memory_stat': memory_stats(),
           'workspace_free_bytes': shutil.disk_usage('/workspace').free}
    with PROFILE.open('a') as f:
        f.write(json.dumps(row, sort_keys=True) + '\n')
        f.flush()
    return row


def terminate_group(proc, wait_seconds=30):
    if proc.poll() is not None:
        return proc.wait()
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        return proc.wait(timeout=wait_seconds)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        return proc.wait()


def write_small_checksums():
    entries = []
    for path in sorted(ROOT.rglob('*')):
        if not path.is_file() or path == ROOT / 'SHA256SUMS.txt':
            continue
        rel = path.relative_to(ROOT)
        if ('checkpoints' in rel.parts or '__pycache__' in rel.parts
                or path.suffix in {'.gpw', '.pyc'} or path.name.endswith('.tmp')):
            continue
        entries.append(f'{sha256(path)}  ./{rel.as_posix()}')
    (ROOT / 'SHA256SUMS.txt').write_text('\n'.join(entries) + '\n')


def archive_to_branch():
    subprocess.run(['git', 'fetch', '--no-tags', 'origin',
                    f'+refs/heads/{BRANCH}:refs/remotes/origin/{BRANCH}'], cwd=REPO, check=True)
    local_head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip()
    remote_head = subprocess.check_output(['git', 'rev-parse', f'refs/remotes/origin/{BRANCH}'], cwd=REPO, text=True).strip()
    if local_head != remote_head:
        raise RuntimeError(f'remote moved during attempt: local={local_head} remote={remote_head}; preserving local outputs')
    write_small_checksums()
    rel_paths = []
    for path in sorted(ROOT.rglob('*')):
        if not path.is_file():
            continue
        rel = path.relative_to(ROOT)
        if ('checkpoints' in rel.parts or '__pycache__' in rel.parts
                or path.suffix in {'.gpw', '.pyc'} or path.name.endswith('.tmp')):
            continue
        rel_paths.append(path.relative_to(REPO).as_posix())
    if not rel_paths:
        raise RuntimeError('no small attempt artifacts to archive')
    subprocess.run(['git', 'add', '--', *rel_paths], cwd=REPO, check=True)
    msg = f'[E] Archive Si surface DFT reference attempt-0002 {TASK}'
    subprocess.run(['git', 'commit', '--only', '-m', msg, '--', *rel_paths], cwd=REPO, check=True)
    subprocess.run(['git', 'push', 'origin', f'HEAD:refs/heads/{BRANCH}'], cwd=REPO, check=True)
    subprocess.run(['git', 'fetch', '--no-tags', 'origin',
                    f'+refs/heads/{BRANCH}:refs/remotes/origin/{BRANCH}'], cwd=REPO, check=True)
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip()
    remote = subprocess.check_output(['git', 'rev-parse', f'refs/remotes/origin/{BRANCH}'], cwd=REPO, text=True).strip()
    if head != remote:
        raise RuntimeError(f'remote tip mismatch after push: local={head} remote={remote}')
    with __import__('tempfile').TemporaryDirectory(prefix='si-ref-archive-check-') as tmp:
        archive = subprocess.check_output(['git', 'archive', '--format=tar', 'HEAD',
                                           ROOT.relative_to(REPO).as_posix()], cwd=REPO)
        import tarfile
        import io
        tarfile.open(fileobj=io.BytesIO(archive), mode='r:').extractall(tmp)
        extracted = Path(tmp) / ROOT.relative_to(REPO)
        subprocess.run(['sha256sum', '-c', 'SHA256SUMS.txt'], cwd=extracted, check=True,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return {'commit': head, 'remote_tip': remote, 'archived_small_files': len(rel_paths),
            'fresh_fetch_sha256_verified': True,
            'checkpoint_binary_uploaded': False}


def execute(expected_lock_sha):
    cfg = json.loads(CONFIG.read_text())
    lock = json.loads(LOCK.read_text())
    input_path = BASE / 'inputs' / f'{TASK}.extxyz'
    attempt_input = ROOT / 'input.extxyz'
    if attempt_input.exists():
        if sha256(attempt_input) != EXPECTED_INPUT_SHA:
            raise RuntimeError('attempt-0002 input already exists with mismatching SHA')
    else:
        shutil.copy2(input_path, attempt_input)
    checkpoint_dir = ROOT / 'checkpoints'
    checkpoint_dir.mkdir(exist_ok=True)
    if not (ROOT / 'RUN_MANIFEST.json').exists():
        atoms = read(input_path, format='extxyz')
        manifest = {
            'structure_id': TASK, 'attempt_id': 'attempt-0002',
            'parent_attempt_id': 'attempt-0001', 'role': lock['role'],
            'input_sha256': EXPECTED_INPUT_SHA,
            'geometry_sha256': cfg['geometry_sha256'],
            'base_config_sha256': lock['base_config_sha256'],
            'diagnostic_config_sha256': lock['config_sha256'],
            'execution_lock_sha256': expected_lock_sha,
            'worker_sha256': lock['worker_sha256'],
            'supervisor_sha256': lock['supervisor_sha256'],
            'source_input_commit': lock['input_package_commit'],
            'C_review_commit': lock['C_review_commit'],
            'created_utc': utc(), 'walltime_limit_seconds': LIMIT,
            'max_scf_iterations': 160, 'checkpoint_enabled': True,
            'checkpoint_policy': lock['checkpoint_policy'], 'automatic_retry': False,
            'config_diff': lock['config_diff'],
            'environment': {
                'hostname': __import__('socket').gethostname(),
                'python_version': sys.version,
                'GPAW_version': gpaw.__version__, 'ASE_version': ase.__version__,
                'gpaw_data_version': importlib.metadata.version('gpaw-data'),
                'GPAW_SETUP_PATH': str(Path(gpaw_data.datapath()).resolve()),
                'GPAW_MPI_BACKEND': 'cgpaw', 'MPI_ranks': 4, 'OMP_NUM_THREADS': 1,
                'cpu_max': Path('/sys/fs/cgroup/cpu.max').read_text().strip(),
                'memory_limit_bytes': read_int('memory.max'),
                'memory_current_bytes_at_start': read_int('memory.current'),
                'workspace_free_bytes_at_start': shutil.disk_usage('/workspace').free,
                'paw_sha256': cfg['paw_sha256']}}
        atomic_json(ROOT / 'RUN_MANIFEST.json', manifest)
    atomic_json(ROOT / 'RUN_STATUS.json', {
        'structure_id': TASK, 'attempt_id': 'attempt-0002', 'status': 'RUNNING',
        'started_utc': utc(), 'walltime_limit_seconds': LIMIT,
        'automatic_retry': False, 'checkpoint_written': False,
        'result_written_only_after_scf_convergence': True,
        'labels_written': False})
    env = os.environ.copy()
    env.update({'GPAW_MPI_BACKEND': 'cgpaw', 'OMP_NUM_THREADS': '1',
                'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1',
                'GPAW_SETUP_PATH': str(Path(gpaw_data.datapath()).resolve()),
                'TMPDIR': str(checkpoint_dir), 'TMP': str(checkpoint_dir),
                'ATTEMPT_ID': 'attempt-0002',
                'TASK_DEADLINE_MONOTONIC': str(time.monotonic() + LIMIT)})
    cmd = [MPI_EXE, '--bind-to', 'core', '--map-by', 'core', '-n', '4',
           str(GPAW_EXE), 'python', '--', str(WORKER),
           '--config', str(CONFIG), '--config-sha256', lock['config_sha256'],
           '--execution-lock', str(LOCK), '--execution-lock-sha256', expected_lock_sha]
    started_mono = time.monotonic()
    started_utc = utc()
    stdout = (ROOT / 'stdout.log').open('wb')
    stderr = (ROOT / 'stderr.log').open('wb')
    proc = subprocess.Popen(cmd, cwd=REPO, env=env, stdin=subprocess.DEVNULL,
                            stdout=stdout, stderr=stderr, start_new_session=True)
    atomic_json(ROOT / 'supervisor_started.json', {
        'utc': started_utc, 'pid': proc.pid, 'command': cmd,
        'walltime_limit_seconds': LIMIT, 'memory_guard_bytes': GUARD,
        'workspace_min_free_bytes': MIN_RUNTIME_FREE})
    stop_reason = None
    stop_requested = False
    min_free = shutil.disk_usage('/workspace').free
    max_mem = read_int('memory.current')
    max_stats = memory_stats()
    sample(0)
    try:
        while proc.poll() is None:
            elapsed = time.monotonic() - started_mono
            row = sample(elapsed)
            min_free = min(min_free, row['workspace_free_bytes'])
            max_mem = max(max_mem, row['memory_current_bytes'])
            if row['memory_current_bytes'] >= GUARD:
                stop_reason = 'CGROUP_MEMORY_GUARD_28_GIB'
            elif row['workspace_free_bytes'] < MIN_RUNTIME_FREE:
                stop_reason = 'DISK_HEADROOM_BELOW_6_GIB'
            elif elapsed >= LIMIT:
                stop_reason = 'WALLTIME_LIMIT_18000_SECONDS'
            elif elapsed >= LIMIT - STOP_RESERVE and not stop_requested:
                (ROOT / 'STOP_REQUEST').write_text('safe walltime stop requested\n')
                stop_requested = True
                atomic_json(ROOT / 'SAFE_STOP_REQUEST.json', {
                    'utc': utc(), 'elapsed_seconds': round(elapsed, 3),
                    'remaining_walltime_seconds': round(LIMIT - elapsed, 3),
                    'reason': 'reserve time to finish one SCF step and preserve latest checkpoint'})
            if stop_reason:
                break
            time.sleep(1)
    except KeyboardInterrupt:
        stop_reason = 'USER_INTERRUPT_NO_AUTOMATIC_RETRY'
    if stop_reason:
        exit_code = terminate_group(proc, 35)
    else:
        exit_code = proc.wait()
    stdout.close()
    stderr.close()
    elapsed = time.monotonic() - started_mono
    partial = checkpoint_dir / 'latest.partial.gpw'
    partial_size = partial.stat().st_size if partial.exists() else 0
    if partial.exists():
        partial.unlink()
    checkpoint_status = None
    cp_meta = ROOT / 'CHECKPOINT_STATUS.json'
    if cp_meta.exists():
        checkpoint_status = json.loads(cp_meta.read_text())
    result_exists = (ROOT / 'RESULT.json').is_file()
    scf_done = False
    conv_path = ROOT / 'SCF_CONVERGENCE.json'
    if conv_path.exists():
        scf_done = json.loads(conv_path.read_text()).get('scf_converged') is True
    if stop_reason is None and (ROOT / 'SCF_STOPPED.json').exists():
        stop_reason = 'SAFE_WALLTIME_STOP_NO_CONVERGENCE'
    if stop_reason is None and exit_code != 0:
        stop_reason = 'SCF_FAILURE_OR_NONCONVERGENCE'
    success = exit_code == 0 and result_exists and scf_done and stop_reason is None
    row = sample(elapsed)
    max_mem = max(max_mem, row['memory_current_bytes'])
    final_status = {
        'structure_id': TASK, 'attempt_id': 'attempt-0002',
        'status': 'SUCCESS' if success else 'FAILED',
        'scf_converged': bool(scf_done), 'exit_code': exit_code,
        'stop_reason': stop_reason, 'elapsed_seconds': round(elapsed, 3),
        'walltime_limit_seconds': LIMIT, 'automatic_retry': False,
        'result_written_only_after_scf_convergence': True,
        'labels_written': bool(result_exists and scf_done),
        'checkpoint_written': checkpoint_status is not None,
        'latest_checkpoint': checkpoint_status.get('path') if checkpoint_status else None,
        'latest_checkpoint_step': checkpoint_status.get('step') if checkpoint_status else None,
        'latest_checkpoint_bytes': checkpoint_status.get('bytes') if checkpoint_status else None,
        'latest_checkpoint_sha256': checkpoint_status.get('sha256') if checkpoint_status else None,
        'latest_checkpoint_precision': checkpoint_status.get('precision') if checkpoint_status else None,
        'partial_checkpoint_removed_bytes': partial_size,
        'minimum_workspace_free_bytes_observed': min_free,
        'maximum_cgroup_memory_current_bytes_observed': max_mem,
        'final_memory_stat': memory_stats(), 'ended_utc': utc()}
    atomic_json(ROOT / 'RUN_STATUS.json', final_status)
    write_small_checksums()
    return final_status


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--lock-sha256', required=True)
    parser.add_argument('--preflight-only', action='store_true')
    args = parser.parse_args()
    checks, details = preflight(args.lock_sha256)
    result = {'status': 'PREFLIGHT_PASS' if all(checks.values()) else 'PREFLIGHT_HOLD',
              'checks': checks, 'details': details,
              'walltime_limit_seconds': LIMIT,
              'memory_guard_bytes': GUARD,
              'minimum_start_free_bytes': MIN_START_FREE,
              'minimum_runtime_free_bytes': MIN_RUNTIME_FREE,
              'checkpoint_policy': 'single-precision mode=all, includes wavefunctions, first step then every 10, atomic latest replacement'}
    atomic_json(ROOT / 'TARGET_PREFLIGHT.json', result)
    print(json.dumps(result, indent=2))
    if not all(checks.values()) or args.preflight_only:
        return 0 if all(checks.values()) else 2
    final_status = execute(args.lock_sha256)
    try:
        archive = archive_to_branch()
    except Exception as exc:
        atomic_json(ROOT / 'ARCHIVE_STATUS.json', {
            'status': 'ARCHIVE_FAILED', 'exception_type': type(exc).__name__,
            'message': str(exc), 'utc': utc()})
        write_small_checksums()
        raise
    print(json.dumps({'RUN_STATUS': final_status, 'ARCHIVE_STATUS': archive}, indent=2))
    return 0 if final_status['status'] == 'SUCCESS' else 3


if __name__ == '__main__':
    raise SystemExit(main())
