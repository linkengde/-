#!/usr/bin/env python3
"""Non-launching resource, runtime and input preflight for the v22 DFT panel."""
from __future__ import annotations
import fcntl,hashlib,importlib.metadata as metadata,json,os,re,shutil,subprocess
from pathlib import Path

ROOT=Path(__file__).resolve().parent
CGROUP=Path('/sys/fs/cgroup')
PYTHON='/workspace/.venvs/gpaw-mpi/bin/python'
LIMIT=34359738368
START_LIMIT=3*1024**3
MIN_DISK=5*1024**3
LOCK=Path('/workspace/.mace-v19-dft.lock')
DEP='/workspace/.local/gpaw-mpi-deps/usr/lib/x86_64-linux-gnu'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
 manifest=json.loads((ROOT/'input_manifest.json').read_text())
 assert len(manifest['preflight_jobs'])==6
 records=manifest['records']
 assert len(records)==6
 assert all(r['geometry_screens']['dev3_geometry_only']['exact_hits']==0 and r['geometry_screens']['dev3_geometry_only']['near_hits']==0 for r in records)
 assert all(r['geometry_screens']['train43_geometry_only']['exact_hits']==0 for r in records)
 for job in manifest['preflight_jobs']:
  p=ROOT/job['input']
  assert p.is_file() and sha(p)==job['input_sha256'],job['label']
  assert job['launch_enabled'] is True and job['kpts']==[6,6,1]
 assert int((CGROUP/'memory.max').read_text())==LIMIT
 memory=int((CGROUP/'memory.current').read_text())
 assert memory<START_LIMIT,f'starting memory {memory} exceeds {START_LIMIT}'
 free=shutil.disk_usage(ROOT).free
 assert free>=MIN_DISK,f'free disk {free} below {MIN_DISK}'
 active=[]
 for proc in Path('/proc').iterdir():
  if proc.name.isdigit():
   try: cmd=(proc/'cmdline').read_bytes().decode(errors='ignore')
   except OSError: continue
   if 'run_preflight_label.py' in cmd and int(proc.name)!=os.getpid(): active.append(proc.name)
 assert not active,f'active DFT runner(s): {active}'
 env=os.environ.copy(); env.update({'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','NUMEXPR_NUM_THREADS':'1'})
 env['LD_LIBRARY_PATH']=':'.join([DEP+'/openmpi/lib',DEP,'/workspace/.local/gpaw-libs/usr/lib/x86_64-linux-gnu','/workspace/.local/openblas/usr/lib/x86_64-linux-gnu/openblas-pthread',env.get('LD_LIBRARY_PATH','')])
 info=subprocess.check_output([str(Path(PYTHON).with_name('gpaw')),'info'],env=env,text=True)
 assert re.search(r'MPI enabled\s+yes',info) and re.search(r'scalapack\s+yes',info),info
 method=manifest['method_candidate']
 probe="import hashlib,importlib.metadata as m,gpaw_data,pathlib; p=pathlib.Path(gpaw_data.datapath()); print(m.version('gpaw'),m.version('ase'),m.version('gpaw-data')); [print(n,hashlib.sha256((p/n).read_bytes()).hexdigest()) for n in "+repr(list(method['PAW_setup_sha256']))+"]"
 rows=subprocess.check_output([PYTHON,'-c',probe],env=env,text=True).splitlines()
 assert rows[0]=='26.7.0 3.29.0 1.2.1',rows
 paws=dict(x.split() for x in rows[1:]); assert paws==method['PAW_setup_sha256'],paws
 with LOCK.open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  fcntl.flock(lock,fcntl.LOCK_UN)
 out={'status':'V22_GEOMETRY_AND_ENVIRONMENT_PREFLIGHT_PASS','candidate_count':6,'train43_sha256':manifest['reference_datasets']['train43_sha256'],'dev3_sha256':manifest['reference_datasets']['dev3_sha256'],'all_input_hashes_pass':True,'all_dev3_geometry_screens_clear':True,'cgroup_memory_limit_bytes':LIMIT,'memory_current_bytes':memory,'minimum_free_disk_bytes':MIN_DISK,'free_disk_bytes':free,'mpi_enabled':True,'scalapack_enabled':True,'GPAW_version':rows[0].split()[0],'ASE_version':rows[0].split()[1],'PAW_setup_hashes_pass':True,'shared_DFT_lock_available':True,'maximum_concurrent_jobs':1,'automatic_retries':0}
 (ROOT/'preflight_environment.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(out,indent=2))

if __name__=='__main__': main()
