import datetime, json, os, signal, subprocess, time
from pathlib import Path
from zoneinfo import ZoneInfo
BASE=Path(__file__).resolve().parent
WORK=Path('/workspace/h-dft-repair')
REL=Path('research/mace-v12-transfer/coordination/reports/chatgpt-nep-ace-dft-repair-20261016')
OUT=BASE/'Ti_gap_32_k1x4x2'
assert not OUT.exists(), 'Refuse to overwrite existing result'
command=['mpiexec','-n','4','/workspace/.venvs/gpaw-mpi/bin/gpaw','python','--',str(WORK/REL/'run_32atom_kmesh_MPIOPT.py'),'--kmesh','1x4x2','--execute','--paw-manifest',str(BASE/'paw_manifest.json'),'--output',str(OUT)]
def timestamp():return datetime.datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(timespec='seconds')
def rss_tree(pid):
 rows=subprocess.check_output(['ps','-eo','pid=,ppid=,rss='],text=True).splitlines()
 procs=[tuple(map(int,line.split())) for line in rows]
 descendants={pid}
 while True:
  updated=descendants|{p for p,pp,r in procs if pp in descendants}
  if updated==descendants:break
  descendants=updated
 return sum(r*1024 for p,pp,r in procs if p in descendants)
started=timestamp();start=time.monotonic();peak=0;reason=None;last=0
with (BASE/'launcher.log').open('x') as log:
 p=subprocess.Popen(command,cwd=WORK,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
 while True:
  status=p.poll();elapsed=time.monotonic()-start
  peak=max(peak,rss_tree(p.pid))
  disk=os.statvfs(BASE);free=disk.f_bavail*disk.f_frsize
  data={'started_at':started,'pid':p.pid,'state':'RUNNING' if status is None else 'EXITED','elapsed_seconds':elapsed,'peak_sampled_process_tree_rss_bytes':peak,'available_disk_bytes':free,'returncode':status,'stop_reason':reason,'command':command,'sampling_interval_seconds':2,'memory_metric_note':'sampled sum RSS of MPI process tree; shared pages may be counted more than once'}
  if status is None and (elapsed>12*3600 or free<2_000_000_000):
   reason='12h wall limit' if elapsed>12*3600 else 'less than 2 GB free disk'
   os.killpg(p.pid,signal.SIGTERM)
   try:p.wait(timeout=20)
   except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait()
   continue
  if status is not None:data['ended_at']=timestamp()
  (BASE/'EXECUTION_MONITOR.json').write_text(json.dumps(data,indent=2)+'\n')
  if status is not None:
   print(json.dumps(data),flush=True);break
  if elapsed-last>=60:
   print(f'RUNNING elapsed={elapsed:.0f}s peak_RSS={peak/1024**3:.2f}GiB free_disk={free/1e9:.2f}GB',flush=True);last=elapsed
  time.sleep(2)
raise SystemExit(status)
