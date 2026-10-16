"""Experimental two-GPW adapter to frozen V4. Offline qualified only, not real-restart validated."""
import hashlib
import os
from pathlib import Path
import sys
import time
import frozen_v4 as base

V4_SHA='23250241f51004380850233cabbbed443be38087c47b286977bdead680ed5c95'
base.need(base.sha256(base.__file__)==V4_SHA, 'FROZEN_V4_MODIFIED')

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

def publish(folder,lock,sequence,writer,metadata):
    """Old GPW is untouched until a fully hashed candidate and atomic pointer exist.
    Failure before pointer publication leaves old pointer usable; orphan files block new writes.
    """
    old=ensure_two_file_budget(folder,lock)
    if old: base.need(sequence>old[0][0]['sequence'],'SEQUENCE_NOT_INCREASING')
    final=folder/('checkpoint.%08d.gpw'%sequence)
    stage=folder/'checkpoint.writing.gpw'
    base.need(not final.exists() and not stage.exists(),'STALE_CHECKPOINT')
    t=time.monotonic(); writer(stage); write_elapsed=time.monotonic()-t
    base.need(stage.is_file() and stage.stat().st_size>0,'EMPTY_CHECKPOINT')
    with stage.open('rb') as f: os.fsync(f.fileno())
    info={**metadata,**lock,'sequence':sequence,'slot':'CURRENT','file':final.name,
          'mode':'all','bytes':stage.stat().st_size,'sha256':base.sha256(stage),
          'checkpoint_finalized':True,'checkpoint_write_elapsed_seconds':write_elapsed,
          'previous_file':old[0][1].name if old else None,
          'previous_sha256':old[0][0]['sha256'] if old else None}
    os.replace(stage,final); base.sync_directory(folder)
    base.write_atomic(folder/('checkpoint.%08d.json'%sequence),info)
    if 'prior_elapsed' in info:
        info['cumulative_elapsed_seconds']=info.pop('prior_elapsed')+time.monotonic()-info.pop('segment_started')
    # Pointer publication is the transaction commit. Never replace the sole old GPW.
    base.write_atomic(folder/'current.json',info)
    return info

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
    base.root_action(world,folder,finish)
    import numpy as np
    elapsed=np.zeros(1)
    if world.rank==0: elapsed[0]=prior_elapsed+time.monotonic()-segment_started
    world.broadcast(elapsed,0)
    return float(elapsed[0])

_original_validate=base.validate_input

def validate(source,config):
    return {**_original_validate(source,config),'relay_adapter_sha256':base.sha256(__file__)}

def main():
    base.validate_input=validate
    base.verified_slots=read_current
    base.save_checkpoint=mpi_save
    # Continuous SCF by default. User wallclock and resource evidence remain mandatory.
    for option in ('--segment-iters','--checkpoint-every'):
        if option not in sys.argv: sys.argv.extend([option,str(2**30)])
    base.main()

if __name__=='__main__': main()
