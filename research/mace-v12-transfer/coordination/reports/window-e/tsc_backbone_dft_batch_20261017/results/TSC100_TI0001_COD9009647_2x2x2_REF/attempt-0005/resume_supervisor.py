#!/usr/bin/env python3
"""One-shot 5-hour, 28-GiB-guarded continuation supervisor."""
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import time

ROOT = Path(__file__).resolve().parent
REPO = Path(subprocess.check_output(['git', 'rev-parse', '--show-toplevel'], cwd=ROOT, text=True).strip())
CG = Path('/sys/fs/cgroup')
CHECKPOINT = ROOT.parent / 'attempt-0003/checkpoints/checkpoint-step-00020.gpw'
EXPECTED = '38d99727d5b5b3f5152e918685e600790b4a36d819400b552b02519ab92eb2cb'
LIMIT = 9100.0
MEMORY_STOP = 30064771072
DISK_STOP = 6 * 1024**3
CACHE_RECLAIM_TRIGGER = 24 * 1024**3
cmd = ['mpiexec', '--bind-to', 'core', '--map-by', 'core', '-n', '4',
       '/workspace/.venvs/gpaw-mpi/bin/gpaw', 'python', '--', str(ROOT / 'resume_from_step20_save_early.py')]
env = os.environ.copy()
env.update({
    'GPAW_MPI_BACKEND': 'cgpaw', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1',
    'GPAW_SETUP_PATH': '/workspace/.venvs/gpaw-mpi/lib/python3.12/site-packages/gpaw_data/setups',
    'LD_LIBRARY_PATH': '/workspace/.local/gpaw-mpi-deps/usr/lib/x86_64-linux-gnu/openmpi/lib:'
                       '/workspace/.local/gpaw-mpi-deps/usr/lib/x86_64-linux-gnu:'
                       '/workspace/.local/gpaw-libs/usr/lib/x86_64-linux-gnu:'
                       '/workspace/.local/openblas/usr/lib/x86_64-linux-gnu/openblas-pthread'
                       + (':' + os.environ['LD_LIBRARY_PATH'] if os.environ.get('LD_LIBRARY_PATH') else ''),
    'TMPDIR': str(ROOT), 'TMP': str(ROOT),
})

def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def memstat():
    out = {}
    for line in (CG / 'memory.stat').read_text().splitlines():
        key, value = line.split()
        if key in {'anon', 'file', 'inactive_file', 'active_file', 'shmem', 'kernel', 'pagetables', 'sock', 'slab', 'file_mapped', 'file_dirty', 'file_writeback'}:
            out[key] = int(value)
    return out

def proc_table():
    table = {}
    for p in Path('/proc').glob('[0-9]*'):
        try:
            lines = (p / 'status').read_text().splitlines()
            ppid = int(next(x.split()[1] for x in lines if x.startswith('PPid:')))
            comm = (p / 'comm').read_text().strip()
            table[int(p.name)] = (ppid, comm, lines)
        except (OSError, ValueError, StopIteration):
            pass
    return table

def descendants(rootpid, table):
    found = {rootpid}
    change = True
    while change:
        change = False
        for pid, (ppid, _, _) in table.items():
            if ppid in found and pid not in found:
                found.add(pid); change = True
    return sorted(found)

def process_memory(rootpid, table):
    rows = []
    for pid in descendants(rootpid, table):
        if pid not in table:
            continue
        ppid, comm, lines = table[pid]
        d = {'pid': pid, 'ppid': ppid, 'comm': comm}
        for line in lines:
            if line.startswith(('VmRSS:', 'RssAnon:', 'RssFile:', 'RssShmem:', 'VmSwap:')):
                key, value, *_ = line.split(); d[key] = int(value) * 1024
        rows.append(d)
    return rows

def atomic_json(path, data):
    tmp = Path(str(path) + '.tmp')
    tmp.write_text(json.dumps(data, indent=2, sort_keys=True) + '\n')
    os.replace(tmp, path)

active = []
for p in Path('/proc').glob('[0-9]*'):
    try:
        if (p / 'stat').read_text().split()[2] == 'Z':
            continue
        if (p / 'comm').read_text().strip() in {'gpaw', 'mpiexec', 'mpirun', 'prterun'}:
            active.append(int(p.name))
    except (OSError, IndexError):
        pass
assert not active, f'concurrent GPAW/MPI tasks: {active}'
assert CHECKPOINT.is_file() and CHECKPOINT.stat().st_size == 6475186475
assert shutil.disk_usage(ROOT).free > DISK_STOP
assert int((CG / 'memory.max').read_text()) == 34359738368
assert int((CG / 'memory.current').read_text()) < MEMORY_STOP
initial_hash = sha256(CHECKPOINT)
assert initial_hash == EXPECTED
fd = os.open(CHECKPOINT, os.O_RDONLY)
try:
    os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
finally:
    os.close(fd)
started = dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')
started_mono = time.monotonic()
max_current = 0
peak_sample = None
per_process_max = {}
reclaim_events = []
stop_reason = None
samples = 0
profile_path = ROOT / 'continuation_memory_profile.jsonl'
with profile_path.open('w') as profile, (ROOT / 'continuation_stdout.log').open('wb') as out, (ROOT / 'continuation_stderr.log').open('wb') as err:
    proc = subprocess.Popen(cmd, cwd=REPO, env=env, stdin=subprocess.DEVNULL,
                            stdout=out, stderr=err, start_new_session=True)
    while proc.poll() is None:
        current = int((CG / 'memory.current').read_text())
        stats = memstat()
        free = shutil.disk_usage(ROOT).free
        rows = process_memory(proc.pid, proc_table())
        samples += 1
        if current > max_current:
            max_current = current
            peak_sample = {'utc': dt.datetime.now(dt.timezone.utc).isoformat(timespec='milliseconds'),
                           'memory_current': current, 'memory_stat': stats,
                           'processes': rows, 'workspace_free': free}
            atomic_json(ROOT / 'continuation_memory_peak_sample.json', peak_sample)
        for row in rows:
            pid = str(row['pid'])
            agg = per_process_max.setdefault(pid, {'comm': row['comm'], 'peak_rss_bytes': 0,
                                                    'peak_anon_bytes': 0, 'peak_file_bytes': 0,
                                                    'peak_shmem_bytes': 0})
            agg['peak_rss_bytes'] = max(agg['peak_rss_bytes'], row.get('VmRSS:', 0))
            agg['peak_anon_bytes'] = max(agg['peak_anon_bytes'], row.get('RssAnon:', 0))
            agg['peak_file_bytes'] = max(agg['peak_file_bytes'], row.get('RssFile:', 0))
            agg['peak_shmem_bytes'] = max(agg['peak_shmem_bytes'], row.get('RssShmem:', 0))
        rec = {'utc': dt.datetime.now(dt.timezone.utc).isoformat(timespec='milliseconds'),
               'memory_current': current, 'memory_stat': stats, 'workspace_free': free, 'processes': rows}
        profile.write(json.dumps(rec, separators=(',', ':')) + '\n'); profile.flush()
        elapsed = time.monotonic() - started_mono
        if current >= CACHE_RECLAIM_TRIGGER and stats.get('inactive_file', 0) > 512 * 1024**2 and not reclaim_events:
            fd = os.open(CHECKPOINT, os.O_RDONLY)
            try:
                os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
            finally:
                os.close(fd)
            reclaim_events.append({'utc': dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds'),
                                   'memory_current_before': current, 'file_before': stats.get('file', 0),
                                   'inactive_file_before': stats.get('inactive_file', 0)})
            stats = memstat(); current = int((CG / 'memory.current').read_text())
            reclaim_events[-1].update({'memory_current_after': current, 'file_after': stats.get('file', 0),
                                       'inactive_file_after': stats.get('inactive_file', 0)})
        if current >= MEMORY_STOP:
            stop_reason = 'MEMORY_LIMIT_GUARD'
        elif free < DISK_STOP:
            stop_reason = 'DISK_HEADROOM_GUARD'
        elif elapsed >= LIMIT:
            stop_reason = 'NEW_5H_WALLTIME_LIMIT_WITH_300S_SHUTDOWN_RESERVE'
        if stop_reason:
            try:
                os.killpg(proc.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                proc.wait(timeout=20)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(proc.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                proc.wait()
            break
        time.sleep(1)
    rc = proc.wait()
ended = dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')
final_hash = sha256(CHECKPOINT)
assert final_hash == EXPECTED, 'source checkpoint changed'
result_path = ROOT / 'FINAL_RESULT.json'
converged = result_path.is_file() and json.loads(result_path.read_text()).get('scf_converged') is True
summary = {
    'status': 'CONVERGED' if converged else ('NOT_CONVERGED' if stop_reason else 'WORKER_FAILED'),
    'exit_code': rc, 'stop_reason': stop_reason or 'PROCESS_EXIT', 'started_utc': started,
    'ended_utc': ended, 'elapsed_seconds': round(time.monotonic() - started_mono, 3),
    'new_window_walltime_limit_seconds': 18000, 'supervisor_limit_seconds': LIMIT,
    'memory_guard_bytes': MEMORY_STOP, 'disk_guard_bytes': DISK_STOP,
    'memory_current_peak_bytes': max_current, 'peak_sample': peak_sample,
    'per_process_peak_rss': per_process_max, 'memory_samples': samples,
    'cache_reclaim_events': reclaim_events, 'source_checkpoint_sha256_before': initial_hash,
    'source_checkpoint_sha256_after': final_hash, 'source_checkpoint_preserved': True,
    'new_intermediate_gpw_written': False, 'logical_step_limit': 160,
    'starting_checkpoint_iteration': 20, 'scf_converged': converged,
    'result_files_present': sorted(p.name for p in ROOT.iterdir() if p.is_file() and p.name in
                                    {'CONVERGED_ENERGY_AND_EVIDENCE.json', 'FORCES_EVA.json', 'STRESS_EVA3.json', 'FINAL_RESULT.json'}),
}
atomic_json(ROOT / 'continuation_memory_summary.json', {
    k: summary[k] for k in ('memory_guard_bytes', 'memory_current_peak_bytes', 'peak_sample',
                              'per_process_peak_rss', 'memory_samples', 'cache_reclaim_events')})
atomic_json(ROOT / 'CONTINUATION_SUPERVISOR_SUMMARY.json', summary)
if not converged:
    atomic_json(ROOT / 'RUN_STATUS.json', {
        'status': summary['status'], 'scf_converged': False,
        'attempt_id': 'attempt-0005', 'source_checkpoint_sha256': EXPECTED,
        'source_checkpoint_preserved': True, 'elapsed_seconds': summary['elapsed_seconds'],
        'stop_reason': summary['stop_reason'], 'labels_written': False,
        'latest_progress': json.loads((ROOT / 'continuation_progress.json').read_text())
                           if (ROOT / 'continuation_progress.json').is_file() else None,
        'ended_utc': ended})
print(json.dumps(summary, indent=2), flush=True)
raise SystemExit(0 if converged else 1)
