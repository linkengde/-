"""Exact registered draft slab control only; disabled until A integration."""
import fcntl,json,os,subprocess,sys
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
out=Path('/workspace/mace_v19_slab_controls')/label;archive=ROOT/'calculations'/label
assert not out.exists() and not archive.exists(),'Existing run/archive: preserve and review'
# Future integration must retain status publication before launch.
sys.path.insert(0,str(REPO/'research/mace-v12-transfer/coordination'));import sync_tasks as sync
sync.progress(r['owner_task'],label,'running',0,None,'Starting exact registered slab PBC/Poisson control; four MPI ranks')
env=os.environ.copy();deps='/workspace/.local/gpaw-mpi-deps/usr/lib/x86_64-linux-gnu'
env.update(GPAW_MPI_BACKEND='cgpaw',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',XDG_CACHE_HOME='/workspace/.cache',MPLCONFIGDIR='/workspace/.cache/matplotlib',FC_CACHEDIR='/workspace/.cache/fontconfig')
env['LD_LIBRARY_PATH']=':'.join([deps+'/openmpi/lib',deps,'/workspace/.local/gpaw-libs/usr/lib/x86_64-linux-gnu','/workspace/.local/openblas/usr/lib/x86_64-linux-gnu/openblas-pthread',env.get('LD_LIBRARY_PATH','')])
py='/workspace/.venvs/gpaw-mpi/bin/python';out.parent.mkdir(exist_ok=True)
subprocess.run(['mpirun','--bind-to','core','--map-by','core','-n','4',py,str(ROOT/'run_pw_reference.py'),label,str(ROOT/r['input']),str(out)],env=env,check=True)
def checksum():
 (archive/'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in sorted(archive.iterdir()) if p.is_file() and p.name!='SHA256SUMS.txt' and p.suffix!='.gpw'))
checksum();v=subprocess.check_output([py,str(ROOT/'verify_result.py'),label,str(archive)],env=env,text=True)
assert json.loads(v)['status']=='PASS';(archive/'verification.json').write_text(v);checksum()
subprocess.run([py,str(ROOT/'verify_result.py'),label,str(archive)],env=env,check=True)
print('Verified draft archive; A integration must wire compact publication and task completion. No implicit broad queue.',flush=True)
