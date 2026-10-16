import datetime, json, os, signal, subprocess, time
from pathlib import Path
from zoneinfo import ZoneInfo
BASE = Path('/workspace/h-dft-runs/p1-20261016')
WORK = Path('/workspace/-')
REL = Path('research/mace-v12-transfer/coordination/reports/chatgpt-nep-ace-dft-repair-20261016')
OUT = BASE / 'Ti_gap_32_k1x6x4'
assert not OUT.exists(), 'Refusing to overwrite existing output'
assert not (BASE / 'launcher.log').exists(), 'Refusing to overwrite prior launcher log'
command = [
    'mpiexec', '-n', '4', '/workspace/.venvs/gpaw-mpi/bin/gpaw', 'python', '--',
    str(WORK / REL / 'run_32atom_kmesh_MPIOPT.py'), '--kmesh', '1x6x4', '--execute',
    '--paw-manifest', str(BASE / 'paw_manifest.json'), '--output', str(OUT)
]
def timestamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
def rss_tree(root_pid):
    rows = subprocess.check_output(['ps', '-eo', 'pid=,ppid=,rss='], text=True).splitlines()
    procs = [tuple(map(int, line.split())) for line in rows]
    descendants = {root_pid}
    while True:
        expanded = descendants | {pid for pid, ppid, rss in procs if ppid in descendants}
        if expanded == descendants:
            break
        descendants = expanded
    return sum(rss * 1024 for pid, ppid, rss in procs if pid in descendants)
started_at = timestamp()
start = time.monotonic()
peak_rss = 0
stop_reason = None
last_report = 0
with (BASE / 'launcher.log').open('x') as log:
    process = subprocess.Popen(command, cwd=WORK, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
    while True:
        status = process.poll()
        elapsed = time.monotonic() - start
        peak_rss = max(peak_rss, rss_tree(process.pid))
        disk = os.statvfs(BASE)
        free = disk.f_bavail * disk.f_frsize
        record = {
            'started_at_utc': started_at,
            'pid': process.pid,
            'state': 'RUNNING' if status is None else 'EXITED',
            'elapsed_seconds': elapsed,
            'peak_sampled_process_tree_rss_bytes': peak_rss,
            'available_disk_bytes': free,
            'returncode': status,
            'stop_reason': stop_reason,
            'command': command,
            'sampling_interval_seconds': 2,
            'memory_metric_note': 'sampled sum RSS of MPI process tree; shared pages may be counted more than once'
        }
        if status is None and (elapsed > 12 * 3600 or free < 2_000_000_000):
            stop_reason = '12h wall limit' if elapsed > 12 * 3600 else 'less than 2 GB free disk'
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=20)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            continue
        if status is not None:
            record['ended_at_utc'] = timestamp()
        (BASE / 'EXECUTION_MONITOR.json').write_text(json.dumps(record, indent=2) + '\n')
        if status is not None:
            print(json.dumps(record), flush=True)
            break
        if elapsed - last_report >= 60:
            print(f'RUNNING elapsed={elapsed:.0f}s peak_RSS={peak_rss / 1024**3:.2f}GiB free_disk={free / 1e9:.2f}GB', flush=True)
            last_report = elapsed
        time.sleep(2)
raise SystemExit(status)
