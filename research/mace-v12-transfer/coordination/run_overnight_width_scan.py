"""Resume registered width jobs in series; preserve unfinished work on failure."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import sync_tasks as sync

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('task', choices=['window-a', 'window-b'])
parser.add_argument('--hours', type=float, default=9)
args = parser.parse_args()
assert 0 < args.hours <= 24
TASK = args.task
sync.ensure_owner(TASK, sync.load_task(TASK))
lock = open('/workspace/.v19-overnight-' + TASK + '.lock', 'a')
fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
deadline = time.time() + args.hours * 3600
BASE = sync.REPO / 'research/mace-v12-transfer'
ROOTS = {width: BASE / ('periodic_interface_v4/' + name) for width, name in [
    ('005', 'pbe_interface_v19_cif_smearing_pilot'),
    ('020', 'pbe_interface_v19_cif_smearing020_pilot')]}
PY = '/workspace/.venvs/gpaw-mpi/bin/python'
k = 5 if TASK == 'window-a' else 6
labels = {width: f'AgSi_COD9009647_pilot_k{k}x{k}_sigma0p' + ('05' if width == '005' else '20') for width in ROOTS}
state_path = Path('/workspace/.setup/v19-overnight-' + TASK + '.json')
state_path.parent.mkdir(exist_ok=True)
deps = '/workspace/.local/gpaw-mpi-deps/usr/lib/x86_64-linux-gnu'
os.environ.update(GPAW_MPI_BACKEND='cgpaw', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1',
                  XDG_CACHE_HOME='/workspace/.cache', MPLCONFIGDIR='/workspace/.cache/matplotlib', FC_CACHEDIR='/workspace/.cache/fontconfig')
os.environ['LD_LIBRARY_PATH'] = ':'.join([deps + '/openmpi/lib', deps, '/workspace/.local/gpaw-libs/usr/lib/x86_64-linux-gnu',
                                       '/workspace/.local/openblas/usr/lib/x86_64-linux-gnu/openblas-pthread', os.environ.get('LD_LIBRARY_PATH', '')])

def state(stage, **detail):
    record = dict(task=TASK, stage=stage, updated_utc=sync.now(), deadline_epoch_s=deadline, **detail)
    state_path.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record), flush=True)

def run(*command):
    return subprocess.check_output(command, cwd=sync.REPO, text=True)

def wait():
    if time.time() >= deadline:
        raise TimeoutError('Monitoring window ended; preserve running jobs and resume this driver later')
    time.sleep(min(30, max(0, deadline - time.time())))

def active():
    for p in Path('/proc').iterdir():
        if not p.name.isdigit():
            continue
        try:
            names = [Path(x.decode(errors='ignore')).name for x in (p / 'cmdline').read_bytes().split(b'\0')]
        except OSError:
            continue
        if any(n in ('run_queue.py', 'run_pw_reference.py', 'run_pw_reference_draft.py') for n in names):
            return True
    return False

def verify(root, label):
    data = json.loads(run(PY, str(root / 'verify_result.py'), label, str(root / 'calculations' / label)))
    assert data['status'] == 'PASS'

def synchronize():
    assert not sync.git('diff', '--name-only').stdout.strip(), 'Tracked changes require diagnosis'
    assert not sync.git('diff', '--cached', '--name-only').stdout.strip(), 'Index changes require diagnosis'
    sync.git('fetch', 'origin', 'main')
    sync.git('merge', '--ff-only', 'origin/main')
    sync.ensure_owner(TASK, sync.load_task(TASK))

try:
    for width, root in ROOTS.items():
        label = labels[width]
        state('WAITING_FOR_EXISTING_JOB_OR_VERIFIED_ARCHIVE', label=label)
        while active():
            wait()
        # Continue after our own archive verifies; another owner is not a
        # prerequisite for the next already-registered fixed-mesh width job.
        synchronize()
        if not (root / 'calculations' / label / 'summary.json').exists():
            local = Path('/workspace') / ('mace_v19_cif_smearing_' + TASK) / label
            assert not local.exists(), 'Unarchived same-job directory exists; recover instead of rerunning'
            assert time.time() < deadline
            state('RUNNING_EXACT_REGISTERED_JOB', label=label)
            subprocess.run(['python3', '-u', str(root / 'run_queue.py'), TASK, label], cwd=sync.REPO, check=True)
        verify(root, label)
        # Queue publishes archives, this additionally marks the task label completed.
        subprocess.run(['python3', str(sync.HERE / 'sync_tasks.py'), 'publish', TASK, '--label', label], cwd=sync.REPO, check=True)
        state('VERIFIED_PUBLISHED', label=label)

    root = ROOTS['020']
    pairs = [(f'AgSi_COD9009647_pilot_k{k}x{k}', labels['005']),
             (f'AgSi_COD9009647_pilot_k{k}x{k}', labels['020']),
             (labels['005'], labels['020'])]
    comparisons = []
    for first, second in pairs:
        comparisons.append(json.loads(run(PY, str(root / 'compare_mesh_results.py'), first, second)))
    subprocess.run(['python3', str(sync.HERE / 'sync_tasks.py'), 'publish', TASK], cwd=sync.REPO, check=True)
    state('OWN_FIXED_MESH_WIDTH_COMPARISONS_PUBLISHED')
    # Independently verify both owner's 0.20 archives before comparing meshes.
    other = f'AgSi_COD9009647_pilot_k{6 if k == 5 else 5}x{6 if k == 5 else 5}_sigma0p20'
    while True:
        synchronize()
        if (root / 'calculations' / other / 'summary.json').exists():
            break
        state('WAITING_FOR_OTHER_OWNER_SIGMA020')
        wait()
    verify(root, other)
    comparisons.append(json.loads(run(PY, str(root / 'compare_mesh_results.py'), 'AgSi_COD9009647_pilot_k5x5_sigma0p05', 'AgSi_COD9009647_pilot_k6x6_sigma0p05')))
    comparisons.append(json.loads(run(PY, str(root / 'compare_mesh_results.py'), 'AgSi_COD9009647_pilot_k5x5_sigma0p20', 'AgSi_COD9009647_pilot_k6x6_sigma0p20')))
    report_dir = sync.HERE / 'reports' / TASK / 'v19_overnight_width_scan'
    report_dir.mkdir(parents=True, exist_ok=True)
    report = {'status': 'COMPARISONS_VERIFIED', 'comparisons': comparisons,
              'scientific_conclusion': 'Assess width and mesh differences separately with unchanged budgets. Broader smearing is not automatically a more accurate target.',
              'limitations': ['One frozen geometry only.', 'No convergence claim from SCF success.', 'No training/validation-label integration.', 'Vacuum/dipole remain a separate controlled check.']}
    (report_dir / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    subprocess.run(['python3', str(sync.HERE / 'sync_tasks.py'), 'publish', TASK], cwd=sync.REPO, check=True)
    state('WIDTH_SCAN_RESULTS_VERIFIED_PUBLISHED')
except Exception as exc:
    state('STOPPED_PRESERVE_SAME_RUN', error=str(exc))
    raise
