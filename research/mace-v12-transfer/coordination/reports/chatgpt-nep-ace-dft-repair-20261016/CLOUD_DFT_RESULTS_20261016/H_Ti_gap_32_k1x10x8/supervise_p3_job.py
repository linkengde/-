#!/usr/bin/env python3
import datetime, json, os, signal, subprocess, time
from pathlib import Path
BASE=Path(__file__).resolve().parent
WORK=Path('/workspace/-')
RUNNER=BASE/'runner_pinned.py'
OUT=BASE/'calculation'
MANIFEST=BASE/'paw_manifest.json'
assert not OUT.exists(), 'Refusing to overwrite existing output'
assert not (BASE/'launcher.log').exists(), 'Refusing to overwrite prior launcher log'
command=['mpiexec','-n','4','/workspace/.venvs/gpaw-mpi/bin/gpaw','python','--',str(RUNNER),'--kmesh','1x10x8','--execute','--paw-manifest',str(MANIFEST),'--output',str(OUT)]
def stamp(): return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
def rss_tree(root):
 rows=subprocess.check_output(['ps','-eo','pid=,ppid=,rss='],text=True).splitlines()
 procs=[tuple(map(int,line.split())) for line in rows]
 descendants={root}
 while True:
  expanded=descendants|{pid for pid,ppid,rss in procs if ppid in descendants}
  if expanded==descendants: break
  descendants=expanded
 return sum(rss*1024 for pid,ppid,rss in procs if pid in descendants)
started=stamp(); start=time.monotonic(); peak=0; reason=None; last=0
with (BASE/'launcher.log').open('x') as log:
 p=subprocess.Popen(command,cwd=WORK,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
 while True:
  rc=p.poll(); elapsed=time.monotonic()-start; peak=max(peak,rss_tree(p.pid))
  sv=os.statvfs(BASE); free=sv.f_bavail*sv.f_frsize
  rec={'started_at_utc':started,'pid':p.pid,'state':'RUNNING' if rc is None else 'EXITED','elapsed_seconds':elapsed,'peak_sampled_process_tree_rss_bytes':peak,'available_disk_bytes':free,'returncode':rc,'stop_reason':reason,'command':command,'sampling_interval_seconds':2,'memory_metric_note':'sampled sum RSS of MPI process tree; shared pages may be counted more than once'}
  if rc is None and (elapsed>12*3600 or free<2_000_000_000):
   reason='12h wall limit' if elapsed>12*3600 else 'less than 2 GB free disk'
   os.killpg(p.pid,signal.SIGTERM)
   try:p.wait(timeout=20)
   except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait()
   continue
  if rc is not None: rec['ended_at_utc']=stamp()
  (BASE/'EXECUTION_MONITOR.json').write_text(json.dumps(rec,indent=2)+'\n')
  if rc is not None:
   print(json.dumps(rec),flush=True);break
  if elapsed-last>=60:
   print(f'RUNNING elapsed={elapsed:.0f}s peak_RSS={peak/1024**3:.2f}GiB free_disk={free/1e9:.2f}GB',flush=True);last=elapsed
  time.sleep(2)
raise SystemExit(rc)
