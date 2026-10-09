"""Exact registered draft slab control only; disabled until A integration."""
import fcntl,json,os,subprocess,sys,time,re
from pathlib import Path
from entry import ROOT,REPO,load,sha
assert len(sys.argv)==2,'usage: run_queue.py EXACT_LABEL'
label=sys.argv[1];m,r=load(label,launch=True)
lock=open('/workspace/.mace-v19-dft.lock','a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
for p in Path('/proc').iterdir():
 if not p.name.isdigit():continue
 try:names=[Path(x.decode(errors='ignore')).name for x in (p/'cmdline').read_bytes().split(b'\0')]
 except OSError:continue
 assert not any(n in ('run_pw_reference.py','run_pw_reference_draft.py') for n in names),'Active GPAW runner: no duplicate'
out=Path('/workspace/mace_v19_gap_scan')/label;archive=ROOT/'calculations'/label
assert not out.exists() and not archive.exists(),'Existing run/archive: preserve and review'
# Future integration must retain status publication before launch.
sys.path.insert(0,str(REPO/'research/mace-v12-transfer/coordination'));import sync_tasks as sync
sync.claim(r['owner_task'])
subprocess.run([sys.executable,str(sync.HERE/'sync_tasks.py'),'progress',r['owner_task'],'--job',label,'--state','running','--iteration','0','--note','Starting registered rigid Ag-layer gap point; four MPI ranks'],check=True)
env=os.environ.copy();deps='/workspace/.local/gpaw-mpi-deps/usr/lib/x86_64-linux-gnu'
env.update(GPAW_MPI_BACKEND='cgpaw',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',XDG_CACHE_HOME='/workspace/.cache',MPLCONFIGDIR='/workspace/.cache/matplotlib',FC_CACHEDIR='/workspace/.cache/fontconfig')
env['LD_LIBRARY_PATH']=':'.join([deps+'/openmpi/lib',deps,'/workspace/.local/gpaw-libs/usr/lib/x86_64-linux-gnu','/workspace/.local/openblas/usr/lib/x86_64-linux-gnu/openblas-pthread',env.get('LD_LIBRARY_PATH','')])
py='/workspace/.venvs/gpaw-mpi/bin/python';out.parent.mkdir(exist_ok=True)
info=subprocess.check_output([str(Path(py).with_name('gpaw')),'info'],env=env,text=True)
assert re.search(r'MPI enabled\s+yes',info) and re.search(r'scalapack\s+yes',info)
launcher=ROOT/'calculations'/f'{label}.launcher.log';launcher.parent.mkdir(exist_ok=True)
with launcher.open('w') as log:
 child=subprocess.Popen(['mpirun','--bind-to','core','--map-by','core','-n','4',py,str(ROOT/'run_pw_reference.py'),label,str(ROOT/r['input']),str(out)],env=env,stdout=log,stderr=subprocess.STDOUT)
 last=0
 while child.poll() is None:
  time.sleep(45)
  progress_file=out/'progress.json'
  if progress_file.exists():
   try:n=json.loads(progress_file.read_text()).get('iteration',0)
   except json.JSONDecodeError:continue
   if n>=last+10:
    subprocess.run([sys.executable,str(sync.HERE/'sync_tasks.py'),'progress',r['owner_task'],'--job',label,'--state','running','--iteration',str(n),'--note','Rigid interface gap SCF active; four MPI ranks'],check=True);last=n
 if child.returncode:
  subprocess.run([sys.executable,str(sync.HERE/'sync_tasks.py'),'progress',r['owner_task'],'--job',label,'--state','failed','--iteration',str(last),'--note',f'Calculator exit{child.returncode}; preserve same run/checkpoint'],check=True)
  raise SystemExit(child.returncode)
def checksum():
 (archive/'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in sorted(archive.iterdir()) if p.is_file() and p.name!='SHA256SUMS.txt' and p.suffix!='.gpw'))
checksum();v=subprocess.check_output([py,str(ROOT/'verify_result.py'),label,str(archive)],env=env,text=True)
assert json.loads(v)['status']=='PASS';(archive/'verification.json').write_text(v);checksum()
subprocess.run([py,str(ROOT/'verify_result.py'),label,str(archive)],env=env,check=True)
subprocess.run([sys.executable,str(sync.HERE/'sync_tasks.py'),'publish',r['owner_task'],'--label',label],check=True)
print('Assigned gap point verified and published',flush=True)
