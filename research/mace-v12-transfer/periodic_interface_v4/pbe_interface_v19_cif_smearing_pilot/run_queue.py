"""Registered-owner V19 queue; fail closed on unfinished runs or publication failure."""
import fcntl,hashlib,json,os,shutil,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[3]
COORD=REPO/'research/mace-v12-transfer/coordination'
sys.path.insert(0,str(COORD));import sync_tasks as sync
TASK=sys.argv[1] if len(sys.argv) in (2,3) else ''
if TASK not in ('window-a','window-b'):raise SystemExit('usage: run_queue.py window-a|window-b [EXACT_REGISTERED_LABEL]')
for line in (ROOT/'SHA256SUMS.txt').read_text().splitlines():
 digest,name=line.split('  ',1);p=ROOT/name
 assert p.resolve().is_relative_to(ROOT.resolve()) and hashlib.sha256(p.read_bytes()).hexdigest()==digest, 'Entry inventory mismatch: '+name
task=sync.load_task(TASK);sync.ensure_owner(TASK,task)
manifest=json.loads((ROOT/'input_manifest.json').read_text());records=[r for r in manifest['records'] if r['owner_task']==TASK]
if len(sys.argv)==3:records=[r for r in records if r['label']==sys.argv[2]]
assert len(records)==1 and manifest['launch_enabled'] is True, 'Choose one exact registered label; never restart the completed Gamma run'
for r in records:assert r['owner_instance']==sync.identity() and r['label'] in task['labels'] and hashlib.sha256((ROOT/r['input']).read_bytes()).hexdigest()==r['input_sha256']
sync.claim(TASK)
lock=open('/workspace/.mace-v19-dft.lock','a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
PYTHON=os.environ.get('GPAW_PYTHON','/workspace/.venvs/gpaw-mpi/bin/python')
env=os.environ.copy();env.update(GPAW_MPI_BACKEND='cgpaw',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',OMPI_ALLOW_RUN_AS_ROOT='1',OMPI_ALLOW_RUN_AS_ROOT_CONFIRM='1')
deps='/workspace/.local/gpaw-mpi-deps/usr/lib/x86_64-linux-gnu'
env['LD_LIBRARY_PATH']=':'.join([deps+'/openmpi/lib',deps,'/workspace/.local/gpaw-libs/usr/lib/x86_64-linux-gnu','/workspace/.local/openblas/usr/lib/x86_64-linux-gnu/openblas-pthread',env.get('LD_LIBRARY_PATH','')])
cores=len(os.sched_getaffinity(0));quota=Path('/sys/fs/cgroup/cpu.max').read_text().split()
if quota[0]!='max':cores=min(cores,int(quota[0])//int(quota[1]))
assert cores>=4
info=subprocess.check_output([str(Path(PYTHON).with_name('gpaw')),'info'],env=env,text=True)
import re
assert re.search(r'MPI enabled\s+yes',info) and re.search(r'scalapack\s+yes',info)
# Refuse coexistence with another actual GPAW runner on this machine.
for p in Path('/proc').iterdir():
 if p.name.isdigit():
  try:args=(p/'cmdline').read_bytes().split(b'\0')
  except OSError:continue
  if any(Path(x.decode(errors='ignore')).name in ('run_pw_reference.py','run_pw_reference_draft.py') for x in args):raise SystemExit('Existing GPAW runner detected; inspect instead of duplicate launch')
runs=Path('/workspace')/('mace_v19_cif_smearing_'+TASK);runs.mkdir(exist_ok=True)
def progress(label,state,iteration,note):
 subprocess.run([sys.executable,str(COORD/'sync_tasks.py'),'progress',TASK,'--job',label,'--state',state,'--iteration',str(iteration),'--note',note],check=True)
def verify(label,folder):return subprocess.check_output([PYTHON,str(ROOT/'verify_result.py'),label,str(folder)],env=env,text=True)
for r in records:
 label=r['label'];archive=ROOT/'calculations'/label;out=runs/label
 if archive.exists():
  if out.exists():raise SystemExit("Archive and local run coexist; inspect same completed job before queue reuse")
  verify(label,archive)
  continue
 else:
  if out.exists():raise SystemExit(f'Unarchived run exists: {out}; inspect/recover SAME job before proceeding')
  if shutil.disk_usage('/workspace').free<manifest["disk_budget"]["minimum_start_bytes"]:raise SystemExit('Insufficient conservative disk reserve; preserve checkpoints and resolve storage before next label')
  progress(label,'running',0,'V19 four-rank PW-PBE starting; fixed geometry, one BLAS thread/rank')
  launcher=ROOT/'calculations'/f'{label}.launcher.log';launcher.parent.mkdir(exist_ok=True)
  with launcher.open('w') as log:
   child=subprocess.Popen(['mpirun','--bind-to','core','--map-by','core','-np','4',PYTHON,str(ROOT/'run_pw_reference.py'),label,str(ROOT/r['input']),str(out)],env=env,stdout=log,stderr=subprocess.STDOUT)
   last=0
   while child.poll() is None:
    time.sleep(45)
    p=out/'progress.json'
    if p.exists():
     try:n=json.loads(p.read_text()).get('iteration',0)
     except json.JSONDecodeError:continue
     if n>=last+10:progress(label,'running',n,'SCF active;4 MPI ranks');last=n
   if child.returncode:progress(label,'failed',last,f'GPAW exit{child.returncode}; preserve same run/checkpoint');raise SystemExit(child.returncode)
  verification=verify(label,archive)
 (archive/'verification.json').write_text(verification)
 files=[p for p in archive.iterdir() if p.is_file() and p.suffix!='.gpw' and p.name!='SHA256SUMS.txt']
 (archive/'SHA256SUMS.txt').write_text(''.join(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n' for p in sorted(files)))
 with sync.LOCK.open('a') as gitlock:
  fcntl.flock(gitlock,fcntl.LOCK_EX)
  sync.git('add','--',str(archive.relative_to(REPO)),str((ROOT/'calculations'/f'{label}.launcher.log').relative_to(REPO)))
  if sync.git('diff','--cached','--name-only').stdout.strip():sync.git('commit','-m',f'Archive verified V19 {label}');sync.push_or_stop(TASK,True)
 summary=json.loads((archive/'summary.json').read_text());progress(label,'completed',summary['scf_iterations'],'Verified finite PW-PBE label and compact archive published; checkpoint retained locally')
print('Assigned CIF convergence pilot complete',flush=True)
