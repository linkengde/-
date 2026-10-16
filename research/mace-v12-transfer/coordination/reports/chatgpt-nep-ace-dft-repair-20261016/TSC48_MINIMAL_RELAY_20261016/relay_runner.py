"""Experimental two-GPW adapter to frozen V4. Offline qualified only, not real-restart validated."""
import hashlib
import os
from pathlib import Path
import subprocess
import sys
import time
import urllib.request
import frozen_v4 as base

V4_SHA='23250241f51004380850233cabbbed443be38087c47b286977bdead680ed5c95'
PREVIOUS_ADAPTER_SHA='fc32f931dc9ac336ee1ab11057293ae7df013b8be781a2de4d4d30a8e4a15dbb'
ROTATION_BRANCH='h-tsc48-local-only-rotation-20261017'
ROTATION_DIR='research/mace-v12-transfer/coordination/reports/chatgpt-nep-ace-dft-repair-20261016/TSC48_LOCAL_ONLY_ROTATION_20261017'
base.need(base.sha256(base.__file__)==V4_SHA, 'FROZEN_V4_MODIFIED')

def bounded_sha256(path):
    """Hash in bounded memory and advise Linux to discard clean read pages."""
    p=Path(path)
    h=hashlib.sha256()
    fd=os.open(p,os.O_RDONLY)
    offset=0
    try:
        while True:
            block=os.read(fd,4*1024*1024)
            if not block: break
            h.update(block)
            if hasattr(os,'posix_fadvise'):
                os.posix_fadvise(fd,offset,len(block),os.POSIX_FADV_DONTNEED)
            offset+=len(block)
    finally:
        os.close(fd)
    return h.hexdigest()

base.sha256=bounded_sha256

def read_current(folder, lock):
    pointer=folder/'current.json'
    if not pointer.exists(): return []
    m=base.read_json(pointer)
    base.need(all(m.get(k)==v for k,v in lock.items()),'RELAY_LOCK_MISMATCH')
    name='checkpoint.%08d.gpw'%m['sequence']
    base.need(m.get('file')==name and m.get('checkpoint_finalized') is True
              and m.get('mode')=='all','RELAY_MANIFEST_INVALID')
    p=folder/name
    base.need(p.is_file() and p.stat().st_size==m['bytes']
              and base.sha256(p)==m['sha256'],'RELAY_CHECKPOINT_CORRUPT')
    return [(m,p)]

def ensure_two_file_budget(folder, lock):
    good=read_current(folder,lock)
    allowed={good[0][1].name} if good else set()
    extras=[p.name for p in folder.glob('*.gpw') if p.name not in allowed]
    base.need(not extras,'PREVIOUS_OR_ORPHAN_GPW_REQUIRES_VERIFIED_RETIREMENT: '+repr(extras))
    return good

def local_only_enabled():
    return os.environ.get('TSC48_LOCAL_ONLY_ROTATION')=='YES'

def git_run(repo,*args):
    p=subprocess.run(['git','-C',str(repo),*args],capture_output=True,text=True)
    base.need(p.returncode==0,'GIT_COMMAND_FAILED '+repr(args)+': '+p.stderr.strip())
    return p.stdout.strip()

def archive_record(record):
    """Commit a small progress record and verify its immutable remote bytes."""
    repo=Path(os.environ.get('TSC48_ARCHIVE_REPO','')).resolve()
    branch=os.environ.get('TSC48_ARCHIVE_BRANCH','')
    base.need(local_only_enabled() and repo.is_dir()
              and branch==ROTATION_BRANCH,'ARCHIVE_REPO_OR_BRANCH_INVALID')
    base.need(git_run(repo,'branch','--show-current')==branch,'ARCHIVE_BRANCH_MISMATCH')
    remote=git_run(repo,'remote','get-url','origin')
    base.need('github.com/linkengde/-' in remote,'ARCHIVE_REPOSITORY_MISMATCH')
    base.need(not git_run(repo,'status','--porcelain'),'ARCHIVE_WORKTREE_DIRTY')
    folder=repo/ROTATION_DIR/'records'
    folder.mkdir(parents=True,exist_ok=True)
    path=folder/(record['event_id']+'.json')
    base.need(not path.exists(),'ARCHIVE_RECORD_ALREADY_EXISTS')
    base.write_atomic(path,record)
    rel=path.relative_to(repo).as_posix()
    git_run(repo,'add','--',rel)
    git_run(repo,'commit','-m','Archive TSC48 local checkpoint '+record['event_id'])
    commit=git_run(repo,'rev-parse','HEAD')
    git_run(repo,'push','origin','HEAD:refs/heads/'+branch)
    ref=git_run(repo,'ls-remote','origin','refs/heads/'+branch).split()
    base.need(len(ref)>=2 and ref[0]==commit,'ARCHIVE_REMOTE_REF_MISMATCH')
    url='https://raw.githubusercontent.com/linkengde/-/'+commit+'/'+rel
    with urllib.request.urlopen(url,timeout=60) as response:
        remote_bytes=response.read()
    local_bytes=path.read_bytes()
    base.need(hashlib.sha256(remote_bytes).digest()==hashlib.sha256(local_bytes).digest(),
              'ARCHIVE_REMOTE_READBACK_MISMATCH')
    return {'commit':commit,'path':rel,'url':url,
            'sha256':hashlib.sha256(local_bytes).hexdigest(),'bytes':len(local_bytes)}

def open_process_holds(path):
    target=path.stat()
    for pid_path in Path('/proc').iterdir():
        if not pid_path.name.isdigit():
            continue
        try:
            entries=list((pid_path/'fd').iterdir())
        except OSError:
            continue
        for fd in entries:
            try:
                st=fd.stat()
            except OSError:
                continue
            if (st.st_dev,st.st_ino)==(target.st_dev,target.st_ino):
                return pid_path.name
    return None

def retire_previous_local(folder,lock,old_meta,new_meta):
    """LOCAL_ONLY_ROTATION: archive evidence, then retire a superseded GPW."""
    if not local_only_enabled() or old_meta is None:
        return
    current=read_current(folder,lock)
    base.need(bool(current) and current[0][0]['sha256']==new_meta['sha256'],
              'NEW_CHECKPOINT_NOT_CURRENT')
    base.need(old_meta['file']=='checkpoint.%08d.gpw'%(new_meta['sequence']-1),
              'NONADJACENT_LOCAL_PREDECESSOR')
    old_path=folder/old_meta['file']
    base.need(not old_path.is_symlink() and old_path.is_file()
              and old_path.stat().st_size==old_meta['bytes']
              and base.sha256(old_path)==old_meta['sha256'],
              'OLD_CHECKPOINT_INVALID_BEFORE_RETIREMENT')
    active_path=folder/'ACTIVE_RUN.json'
    base.need(active_path.is_file(),'ACTIVE_RUN_LOCK_MISSING')
    active=base.read_json(active_path)
    base.need(active.get('latest_checkpoint_sha256',
                        active.get('source_checkpoint_sha256'))==old_meta['sha256'],
              'ACTIVE_RUN_STILL_REFERENCES_DIFFERENT_CHECKPOINT')
    active['latest_checkpoint_sha256']=new_meta['sha256']
    active['latest_checkpoint_sequence']=new_meta['sequence']
    active['updated_utc']=base.datetime.datetime.now(
        base.datetime.timezone.utc).isoformat()
    base.write_atomic(active_path,active)
    before_free=__import__('shutil').disk_usage(folder).free
    event_time=base.datetime.datetime.now(base.datetime.timezone.utc).isoformat()
    pre={'event_id':'retire-before-%08d'%old_meta['sequence'],
         'event':'LOCAL_ONLY_RETIRE_BEFORE',
         'accepted_risk':'H1_LOCAL_ONLY_WINDOW_LOSS_MAY_PREVENT_RECOVERY',
         'retiring_checkpoint':old_meta,'replacement_checkpoint':new_meta,
         'disk_free_bytes_before':before_free,'input_sha256':lock['input_sha256'],
         'h_config_sha256':lock['h_config_sha256'],
         'runner_sha256':lock['runner_sha256'],
         'local_rotation_adapter_sha256':lock['local_rotation_adapter_sha256'],
         'utc':event_time}
    pre_receipt=archive_record(pre)
    base.need(open_process_holds(old_path) is None,'OLD_CHECKPOINT_OPEN_BY_PROCESS')
    base.need(not (folder/'checkpoint.writing.gpw').exists(),
              'CHECKPOINT_WRITE_STILL_ACTIVE')
    current=read_current(folder,lock)
    base.need(current[0][1].name==new_meta['file']
              and current[0][0]['sha256']==new_meta['sha256'],
              'CURRENT_POINTER_CHANGED_BEFORE_RETIREMENT')
    active=base.read_json(active_path)
    base.need(active.get('latest_checkpoint_sha256')==new_meta['sha256'],
              'ACTIVE_RUN_POINTER_NOT_ADVANCED')
    allocated=old_path.stat().st_blocks*512
    old_path.unlink()
    base.sync_directory(folder)
    after_free=__import__('shutil').disk_usage(folder).free
    retired_utc=base.datetime.datetime.now(base.datetime.timezone.utc).isoformat()
    current_meta=current[0][0]
    current_meta.update({'previous_file':None,'previous_sha256':None,
                         'retired_predecessor_sha256':old_meta['sha256'],
                         'local_only_rotation':True,
                         'retired_predecessor_utc':retired_utc,
                         'retirement_archive_commit':pre_receipt['commit']})
    base.write_atomic(folder/'current.json',current_meta)
    base.write_atomic(folder/('checkpoint.%08d.json'%new_meta['sequence']),current_meta)
    base.append_ledger(folder,{'event':'LOCAL_ONLY_CHECKPOINT_RETIRED',
        'retired_sha256':old_meta['sha256'],'retired_bytes':old_meta['bytes'],
        'replacement_sha256':new_meta['sha256'],'allocated_bytes_estimate':allocated,
        'disk_free_bytes_before':before_free,'disk_free_bytes_after':after_free,
        'disk_free_delta_bytes':after_free-before_free,
        'pre_delete_archive_commit':pre_receipt['commit'],'utc':retired_utc})
    post={'event_id':'retire-after-%08d'%old_meta['sequence'],
          'event':'LOCAL_ONLY_RETIRE_AFTER',
          'accepted_risk':'H1_LOCAL_ONLY_WINDOW_LOSS_MAY_PREVENT_RECOVERY',
          'retired_checkpoint_sha256':old_meta['sha256'],
          'retired_checkpoint_bytes':old_meta['bytes'],
          'allocated_bytes_estimate':allocated,
          'replacement_checkpoint_sha256':new_meta['sha256'],
          'replacement_checkpoint_sequence':new_meta['sequence'],
          'scf_steps':new_meta['cumulative_logged_scf_steps'],
          'disk_free_bytes_before':before_free,'disk_free_bytes_after':after_free,
          'disk_free_delta_bytes':after_free-before_free,
          'pre_delete_archive_commit':pre_receipt['commit'],
          'input_sha256':lock['input_sha256'],'h_config_sha256':lock['h_config_sha256'],
          'runner_sha256':lock['runner_sha256'],
          'local_rotation_adapter_sha256':lock['local_rotation_adapter_sha256'],
          'utc':retired_utc}
    archive_record(post)

def retire_previous(folder,lock,receipt_path):
    """Explicit maintenance after externally verified remote roundtrip, while no runner active.
    Receipt is an auditable operator attestation, not independent proof of remote availability.
    Never removes the current GPW; caller must not fabricate an attestation.
    """
    base.need(not (folder/'ACTIVE_RUN.json').exists(),'RUNNER_ACTIVE')
    m,p=read_current(folder,lock)[0]; r=base.read_json(receipt_path)
    base.need(r.get('independent_reader_verified') is True and r.get('sha256')==m['sha256']
              and r.get('bytes')==m['bytes'] and bool(r.get('reader_identity'))
              and bool(r.get('remote_immutable_id')),'REMOTE_RECEIPT_NOT_VERIFIED')
    old=m.get('previous_file')
    if old:
        base.need(old=='checkpoint.%08d.gpw'%(m['sequence']-1),'REVIEW_NONADJACENT_PREDECESSOR')
        q=folder/old
        if q.exists():
            base.need(base.sha256(q)==m['previous_sha256'],'OLD_CHECKPOINT_MISMATCH')
            q.unlink(); base.sync_directory(folder)

def mpi_save(calc,world,folder,lock,sequence,iteration,source_sha,segment_number,
             segment_started,prior_elapsed,cumulative_steps,log_name):
    # All ranks inspect the same shared filesystem, then collectively write exactly once.
    base.root_action(world,folder,lambda:ensure_two_file_budget(folder,lock))
    stage=folder/'checkpoint.writing.gpw'
    base.collective_stage_clear(world,stage)
    start=time.monotonic()
    try: calc.write(str(stage),mode='all')
    except Exception:
        if hasattr(world,'abort'): world.abort(9)
        raise
    finished=time.monotonic(); world.barrier()
    def finish():
        # publish expects a writer: already collectively written file is moved aside only
        # logically, with no additional GPW copy and no extra disk allocation.
        old=read_current(folder,lock)
        final=folder/('checkpoint.%08d.gpw'%sequence)
        base.need(not final.exists() and stage.stat().st_size>0,'INVALID_STAGE')
        with stage.open('rb') as f: os.fsync(f.fileno())
        m={**lock,'sequence':sequence,'slot':'CURRENT','file':final.name,'mode':'all',
           'sha256':base.sha256(stage),'bytes':stage.stat().st_size,'checkpoint_finalized':True,
           'checkpoint_write_elapsed_seconds':finished-start,'iteration_in_segment':iteration,
           'source_checkpoint_sha256':source_sha,'segment_number':segment_number,
           'cumulative_logged_scf_steps':cumulative_steps,'scf_log_file':log_name,
           'previous_file':old[0][1].name if old else None,
           'previous_sha256':old[0][0]['sha256'] if old else None}
        os.replace(stage,final); base.sync_directory(folder)
        base.write_atomic(folder/('checkpoint.%08d.json'%sequence),m)
        base.append_ledger(folder,{'event':'RELAY_CHECKPOINT_PUBLISHING',**m})
        m['cumulative_elapsed_seconds']=prior_elapsed+time.monotonic()-segment_started
        base.write_atomic(folder/'current.json',m)
        if local_only_enabled():
            retire_previous_local(folder,lock,old[0][0] if old else None,m)
    base.root_action(world,folder,finish)
    import numpy as np
    elapsed=np.zeros(1)
    if world.rank==0: elapsed[0]=prior_elapsed+time.monotonic()-segment_started
    world.broadcast(elapsed,0)
    return float(elapsed[0])

_original_validate=base.validate_input

def validate(source,config):
    return {**_original_validate(source,config),
            'relay_adapter_sha256':PREVIOUS_ADAPTER_SHA,
            'local_rotation_adapter_sha256':base.sha256(__file__)}

def compatible_factory(factory, *args, **kwargs):
    # GPAW 26.7 parameter objects fall through Eigensolver.from_param to legacy.
    # Preserve every PPCG default explicitly while expressing it as a named dict.
    from gpaw.dft import PPCG
    solver=kwargs.get('eigensolver')
    if isinstance(solver,PPCG):
        kwargs['eigensolver']={'name':'ppcg',**solver.todict()}
    calc=factory(*args,**kwargs)
    base.need(hasattr(calc,'hooks'),'SCF_HOOK_INTERFACE_UNAVAILABLE')
    return calc

_original_decision=base.shared_decision

def periodic_decision(folder, world, *args):
    # Pause AFTER writing the second file, before any third-file requirement.
    # Continuing regular saves requires a separately verified remote retirement.
    import numpy as np
    save,stop,blocked=_original_decision(world,*args)
    pause=np.zeros(1,dtype=np.int32)
    if world.rank==0 and save and not blocked:
        pause[0]=int(len(list(folder.glob('checkpoint.*.gpw')))==1)
    world.broadcast(pause,0)
    return save,stop or bool(pause[0]),blocked

def main():
    base.validate_input=validate
    base.verified_slots=read_current
    base.save_checkpoint=mpi_save
    import gpaw
    factory=gpaw.GPAW
    gpaw.GPAW=lambda *a,**kw:compatible_factory(factory,*a,**kw)
    # Keep SCF iteration ceiling unchanged; require an explicit runtime checkpoint
    # cadence selected from measured write time and actual external wallclock.
    if '--segment-iters' not in sys.argv:sys.argv.extend(['--segment-iters',str(2**30)])
    if len(sys.argv)>1 and sys.argv[1] in ('start','resume'):
        base.need('--checkpoint-every' in sys.argv,'EXPLICIT_CHECKPOINT_CADENCE_REQUIRED')
        folder=Path(sys.argv[sys.argv.index('--workdir')+1])
        if local_only_enabled():
            base.need(os.environ.get('TSC48_ARCHIVE_REPO') and
                      os.environ.get('TSC48_ARCHIVE_BRANCH')==ROTATION_BRANCH,
                      'LOCAL_ROTATION_ARCHIVE_NOT_CONFIGURED')
            base.shared_decision=_original_decision
        else:
            base.shared_decision=lambda world,*args:periodic_decision(folder,world,*args)
    base.main()

if __name__=='__main__': main()
