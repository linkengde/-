#!/usr/bin/env python3
"""Run registered v25 Ti3SiC2 bulk labels sequentially with memory and lock guards."""
from __future__ import annotations
import fcntl, hashlib, json, os, re, shutil, signal, subprocess, sys, time
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parent
REPO=Path('/workspace/jeitschko-main-run')
PYTHON='/workspace/.venvs/gpaw-mpi/bin/python'
CGROUP=Path('/sys/fs/cgroup')
LIMIT=32*1024**3
STOP=30*1024**3
LOCK=Path('/workspace/.mace-v19-dft.lock')
MANIFEST=json.loads((ROOT/'input_manifest.json').read_text())
STATUS_PATH=ROOT/'sequence_status.json'
CALC=ROOT/'calculations'
CHECKPOINT=CALC/'latest_checkpoint'
PAUSE=ROOT/'pause.request'

def utc(): return datetime.now(timezone.utc).isoformat()
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def mem(): return int((CGROUP/'memory.current').read_text())
def active_runners():
    result=[]
    for p in Path('/proc').iterdir():
        if not p.name.isdigit(): continue
        try:
            cmd=(p/'cmdline').read_bytes().decode(errors='ignore')
            state=(p/'stat').read_text().rsplit(')',1)[1].strip().split()[0]
        except (OSError,IndexError): continue
        if 'run_preflight_label.py' in cmd and state!='Z': result.append(int(p.name))
    return result
def write_status(d):
    tmp=STATUS_PATH.with_suffix('.json.tmp')
    tmp.write_text(json.dumps(d,indent=2)+'\n')
    tmp.replace(STATUS_PATH)
def write_sums(folder):
    rows=[]
    for p in sorted(folder.iterdir()):
        if p.is_file() and p.name!='SHA256SUMS.txt': rows.append(f'{digest(p)}  {p.name}\n')
    (folder/'SHA256SUMS.txt').write_text(''.join(rows))
def verify(label):
    p=subprocess.run([PYTHON,str(ROOT/'verify_label.py'),label],cwd=REPO,text=True,capture_output=True)
    if p.returncode: raise RuntimeError(p.stderr or p.stdout)
    return json.loads(p.stdout)
def move_failed(label,reason,code,started,peak):
    scratch=CALC/f'.running_{label}'
    if scratch.exists():
        (scratch/'state.next.gpw').unlink(missing_ok=True)
        failed_root=CALC/'failed_attempts'; failed_root.mkdir(exist_ok=True)
        failed=failed_root/f'{label}__attempt_{int(started)}'
        if failed.exists(): failed=failed_root/f'{label}__attempt_{int(started)}_{int(time.time())}'
        scratch.rename(failed)
        (failed/'runtime.json').write_text(json.dumps({
            'label':label,'return_code':code,'reason':reason,'elapsed_s':time.time()-started,
            'memory_high_water_observed_bytes':peak,'cgroup_limit_bytes':LIMIT,'memory_stop_threshold_bytes':STOP,
        },indent=2)+'\n')
def preflight(env):
    if int((CGROUP/'memory.max').read_text())!=LIMIT: raise RuntimeError('unexpected cgroup memory limit')
    if mem()>=4*1024**3: raise RuntimeError(f'prelaunch memory exceeds 4 GiB ({mem()})')
    if shutil.disk_usage(ROOT).free<5*1024**3: raise RuntimeError('less than 5 GiB free disk')
    if active_runners(): raise RuntimeError(f'another DFT runner is active: {active_runners()}')
    info=subprocess.check_output(['/workspace/.venvs/gpaw-mpi/bin/gpaw','info'],env=env,text=True)
    if not re.search(r'MPI enabled\s+yes',info) or not re.search(r'scalapack\s+yes',info): raise RuntimeError('GPAW MPI/ScaLAPACK preflight failed')
    probe="import hashlib,importlib.metadata as m,gpaw_data,pathlib; p=pathlib.Path(gpaw_data.datapath()); print(m.version('gpaw'),m.version('ase'),m.version('gpaw-data')); [print(n,hashlib.sha256((p/n).read_bytes()).hexdigest()) for n in " + repr(list(MANIFEST['method']['PAW_setup_sha256'])) + "]"
    rows=subprocess.check_output([PYTHON,'-c',probe],env=env,text=True).splitlines()
    if rows[0] != '26.7.0 3.29.0 1.2.1': raise RuntimeError(f'unexpected GPAW environment: {rows[0]}')
    setups=dict(x.split() for x in rows[1:])
    if setups!=MANIFEST['method']['PAW_setup_sha256']: raise RuntimeError('PAW setup hashes changed')

def run_job(job,status):
    label=job['label']; target=CALC/label; scratch=CALC/f'.running_{label}'
    if target.exists():
        result=verify(label)
        receipt=target/'verification.json'
        if not receipt.exists():
            receipt.write_text(json.dumps(result,indent=2)+'\n')
            write_sums(target)
            result=verify(label)
        status.setdefault('completed_labels',[]).append({'label':label,'archive_verification':'PASS','archive_sha256':result['archive_sha256']})
        status['current_label']=None; status['last_reported_scf_iteration']=None; status['last_job_progress']=None
        write_status(status); return
    if scratch.exists():
        failed_root=CALC/'failed_attempts'; failed_root.mkdir(exist_ok=True)
        old=failed_root/f'{label}__stale_{int(time.time())}'; scratch.rename(old)
    status['current_label']=label; status['status']='running'; write_status(status)
    CALC.mkdir(parents=True,exist_ok=True)
    env=os.environ.copy()
    deps='/workspace/.local/gpaw-mpi-deps/usr/lib/x86_64-linux-gnu'
    env.update({'GPAW_MPI_BACKEND':'cgpaw','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','NUMEXPR_NUM_THREADS':'1'})
    env['LD_LIBRARY_PATH']=':'.join([deps+'/openmpi/lib',deps,'/workspace/.local/gpaw-libs/usr/lib/x86_64-linux-gnu','/workspace/.local/openblas/usr/lib/x86_64-linux-gnu/openblas-pthread',env.get('LD_LIBRARY_PATH','')])
    env.pop('V25_RESTART_CHECKPOINT',None); env.pop('V25_RESTART_ITERATION_OFFSET',None)
    cp=CHECKPOINT/'latest_state.gpw'; meta_path=CHECKPOINT/'latest_checkpoint.json'
    if cp.is_file() and meta_path.is_file():
        try:
            meta=json.loads(meta_path.read_text())
            if meta.get('label')==label and meta.get('input_sha256')==job['input_sha256'] and digest(cp)==meta.get('checkpoint_sha256'):
                env['V25_RESTART_CHECKPOINT']=str(cp)
                env['V25_RESTART_ITERATION_OFFSET']=str(meta.get('last_reported_scf_iteration') or 0)
        except (OSError,ValueError,json.JSONDecodeError): pass
    preflight(env)
    log_path=CALC/f'{label}.supervisor.log'
    started=time.time(); peak=mem(); stopped_memory=False
    with LOCK.open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if active_runners(): raise RuntimeError('another DFT runner acquired the slot after preflight')
        with log_path.open('w') as log:
            child=subprocess.Popen(['mpirun','--bind-to','core','--map-by','core','-n','4',PYTHON,str(ROOT/'run_preflight_label.py'),label,str(scratch)],cwd=REPO,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            while child.poll() is None:
                current=mem(); peak=max(peak,current)
                if current>=STOP:
                    stopped_memory=True; os.killpg(child.pid,signal.SIGTERM)
                    try: child.wait(timeout=30)
                    except subprocess.TimeoutExpired: os.killpg(child.pid,signal.SIGKILL)
                    break
                try:
                    progress=json.loads((scratch/'progress.json').read_text())
                    status['last_reported_scf_iteration']=progress.get('iteration',0)
                    status['last_job_progress']=progress.get('status')
                    write_status(status)
                    if progress.get('status')=='paused_checkpointed':
                        os.killpg(child.pid,signal.SIGTERM)
                except (OSError,json.JSONDecodeError): pass
                time.sleep(10)
            code=child.wait()
    try: progress=json.loads((scratch/'progress.json').read_text())
    except (OSError,json.JSONDecodeError): progress={}
    if progress.get('status')=='paused_checkpointed' or PAUSE.exists():
        move_failed(label,'user_requested_pause',code,started,peak)
        status.update({'status':'paused_by_user','current_label':label,'pause':{
            'label':label,'last_reported_scf_iteration':progress.get('iteration'),
            'checkpoint_available':bool(cp.is_file() and meta_path.is_file() and json.loads(meta_path.read_text()).get('label')==label),
            'requested_utc':utc()}})
        write_status(status)
        raise SystemExit(0)
    if stopped_memory or code!=0:
        move_failed(label,'memory_stop_guard' if stopped_memory else 'DFT_runner_nonzero_exit',code,started,peak)
        status.update({'status':'failed_stop_no_retry','failed_label':label,'current_label':label,'launcher_exit_code':code,'finished_utc':utc()})
        write_status(status)
        raise RuntimeError(f'DFT failed or hit memory guard for {label}; preserved attempt and checkpoint')
    if not (scratch/'summary.json').is_file(): raise RuntimeError(f'missing completed summary for {label}')
    (scratch/'state.gpw').unlink(missing_ok=True)
    (scratch/'state.next.gpw').unlink(missing_ok=True)
    runtime={'label':label,'return_code':0,'elapsed_s':time.time()-started,'memory_high_water_observed_bytes':peak,'cgroup_limit_bytes':LIMIT,'memory_stop_threshold_bytes':STOP}
    (scratch/'runtime.json').write_text(json.dumps(runtime,indent=2)+'\n')
    scratch.rename(target); write_sums(target)
    first=verify(label)
    (target/'verification.json').write_text(json.dumps(first,indent=2)+'\n')
    write_sums(target); result=verify(label)
    status.setdefault('completed_labels',[]).append({'label':label,'archive_verification':'PASS','archive_sha256':result['archive_sha256']})
    status['current_label']=None; status['last_reported_scf_iteration']=None; status['last_job_progress']=None
    write_status(status)

def main():
    args=set(sys.argv[1:])
    status=json.loads(STATUS_PATH.read_text()) if STATUS_PATH.exists() else {'status':'inputs_prepared','completed_labels':[]}
    if PAUSE.exists() and '--resume' not in args: raise SystemExit('pause.request exists; use resume_sequence.py to continue')
    if '--resume' in args:
        PAUSE.unlink(missing_ok=True)
        status['status']='running'; status.pop('pause',None)
    if status.get('status')=='paused_by_user' and '--resume' not in args: raise SystemExit('sequence is paused; run resume_sequence.py explicitly')
    status.setdefault('completed_labels',[])
    status['started_or_resumed_utc']=utc(); status['status']='running'; write_status(status)
    CALC.mkdir(parents=True,exist_ok=True)
    with LOCK.open('a') as lock:
        # The job-level lock is taken by run_job; this open records that the lock path is available.
        pass
    for job in MANIFEST['preflight_jobs']:
        completed={r['label'] for r in status['completed_labels'] if r.get('archive_verification')=='PASS'}
        if job['label'] in completed: continue
        run_job(job,status)
        if status.get('status')=='paused_by_user': return
    status.update({'status':'all_skeleton_archives_verified','current_label':None,'finished_utc':utc()})
    write_status(status)

if __name__=='__main__':
    try: main()
    except SystemExit: raise
    except Exception as exc:
        try:
            d=json.loads(STATUS_PATH.read_text()); d.update({'status':'failed_stop_no_retry','supervisor_error':str(exc),'finished_utc':utc()}); write_status(d)
        except Exception: pass
        print(f'STOPPED WITHOUT RETRY: {exc}',file=sys.stderr,flush=True)
        raise
