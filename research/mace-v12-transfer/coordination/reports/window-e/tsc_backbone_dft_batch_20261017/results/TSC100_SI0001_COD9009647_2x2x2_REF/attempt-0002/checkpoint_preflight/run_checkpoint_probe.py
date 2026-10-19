#!/usr/bin/env python3
"""Supervise serial GPAW checkpoint write/read probes under the run guards."""
import datetime as dt
import hashlib
import json
import os
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RUN_DIR = ROOT.parent
PROBE = ROOT / 'checkpoint_probe.py'
SOURCE = Path('/workspace/-/research/mace-v12-transfer/coordination/reports/window-e/tsc_backbone_dft_batch_20261017/results/TSC100_TI0001_COD9009647_2x2x2_REF/attempt-0003/checkpoints/checkpoint-step-00020.gpw')
TARGET = ROOT / 'single_precision_full_wfs_test.gpw'
GUARD = 28 * 1024**3
DISK_HEADROOM = 6 * 1024**3
PROFILE = ROOT / 'memory_profile.jsonl'


def cgroup_int(name):
    return int(Path('/sys/fs/cgroup', name).read_text().strip())


def sha256_file(path):
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


def stat_row():
    d = {}
    for line in Path('/sys/fs/cgroup/memory.stat').read_text().splitlines():
        k, v = line.split()
        d[k] = int(v)
    return {'utc': dt.datetime.now(dt.timezone.utc).isoformat(),
            'memory_current': cgroup_int('memory.current'),
            'memory_stat': {k: d.get(k, 0) for k in
                            ('anon', 'file', 'inactive_file', 'active_file', 'shmem', 'kernel', 'slab')},
            'workspace_free': shutil.disk_usage('/workspace').free}


def run_stage(stage):
    gpaw_exe = shutil.which('gpaw')
    mpi_exe = shutil.which('mpiexec')
    if not gpaw_exe or not mpi_exe:
        raise RuntimeError('gpaw or mpiexec unavailable in activated venv')
    env = os.environ.copy()
    env.update({'GPAW_MPI_BACKEND': 'cgpaw', 'OMP_NUM_THREADS': '1'})
    cmd = [mpi_exe, '-n', '4', gpaw_exe, 'python', '--', str(PROBE),
           stage, str(SOURCE), str(TARGET)]
    out = (ROOT / f'{stage}.stdout').open('wb')
    err = (ROOT / f'{stage}.stderr').open('wb')
    proc = subprocess.Popen(cmd, cwd=RUN_DIR, env=env, stdout=out, stderr=err,
                            start_new_session=True)
    started = time.monotonic()
    reason = None
    while proc.poll() is None:
        row = stat_row()
        with PROFILE.open('a') as f:
            f.write(json.dumps({'stage': stage, **row}) + '\n')
        if row['memory_current'] >= GUARD:
            reason = 'CGROUP_MEMORY_GUARD_28_GIB'
        elif row['workspace_free'] < DISK_HEADROOM:
            reason = 'WORKSPACE_HEADROOM_BELOW_6_GIB'
        elif time.monotonic() - started > 1800:
            reason = 'CHECKPOINT_PREFLIGHT_TIMEOUT_1800_SECONDS'
        if reason:
            os.killpg(proc.pid, signal.SIGTERM)
            try:
                proc.wait(timeout=15)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()
            break
        time.sleep(0.5)
    code = proc.wait()
    out.close()
    err.close()
    return {'stage': stage, 'exit_code': code, 'stop_reason': reason,
            'elapsed_seconds': round(time.monotonic() - started, 3),
            'peak_memory_current_bytes': max(
                [json.loads(x)['memory_current'] for x in PROFILE.read_text().splitlines()
                 if json.loads(x)['stage'] == stage] or [cgroup_int('memory.current')]),
            'peak_workspace_free_bytes': min(
                [json.loads(x)['workspace_free'] for x in PROFILE.read_text().splitlines()
                 if json.loads(x)['stage'] == stage] or [shutil.disk_usage('/workspace').free])}


result = {'started_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
          'source_checkpoint': str(SOURCE),
          'source_sha256_expected': '38d99727d5b5b3f5152e918685e600790b4a36d819400b552b02519ab92eb2cb',
          'checkpoint_mode': 'all', 'checkpoint_precision': 'single',
          'memory_guard_bytes': GUARD, 'disk_headroom_bytes': DISK_HEADROOM}
if not SOURCE.is_file():
    result.update(status='HOLD', reason='source full-wavefunction checkpoint missing')
else:
    write_result = run_stage('write')
    result['write'] = write_result
    if write_result['exit_code'] == 0 and write_result['stop_reason'] is None and TARGET.is_file():
        read_result = run_stage('read')
        result['read'] = read_result
        if read_result['exit_code'] == 0 and read_result['stop_reason'] is None:
            result.update(status='PASS', wavefunctions_reloaded=True,
                          test_checkpoint_sha256=sha256_file(TARGET),
                          test_checkpoint_bytes=TARGET.stat().st_size)
        else:
            result.update(status='HOLD', reason='checkpoint readback failed or hit a guard')
    else:
        result.update(status='HOLD', reason='checkpoint write failed or hit a guard')
result['ended_utc'] = dt.datetime.now(dt.timezone.utc).isoformat()
(ROOT / 'checkpoint_preflight.json').write_text(json.dumps(result, indent=2) + '\n')
if TARGET.exists():
    TARGET.unlink()
print(json.dumps(result, indent=2))
sys.exit(0 if result['status'] == 'PASS' else 2)
