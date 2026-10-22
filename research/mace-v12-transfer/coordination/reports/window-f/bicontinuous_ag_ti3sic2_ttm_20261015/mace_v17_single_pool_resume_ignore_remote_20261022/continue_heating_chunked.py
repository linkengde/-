#!/usr/bin/env python3
"""Memory-bounded MACE V17 heating continuation from the last saved checkpoint.

The original trajectory/checkpoints are untouched. Resume starts from the saved
0.4 ps heating state and restores its NumPy Langevin RNG state. Each worker
process handles 0.05 ps, then exits so memory is reclaimed before the next chunk.
"""
from __future__ import annotations
import json, os, subprocess, sys, time
from pathlib import Path

ROOT=Path(__file__).resolve().parent
SCRIPT=ROOT/'run_compact_mace.py'
CHECK=ROOT/'HEATING_RESUME_CHECKPOINT.extxyz'
STATE=ROOT/'HEATING_RESUME_STATE.json'
OUTPUT=ROOT/'HEAT_COOL_OVITO_RESUMED.extxyz'
CHUNK=100
HEAT_MAX=8000
COOL_MAX=4000
WALL_LIMIT=None  # no wall-clock cap; the physical 4 ps heat cap remains


def json_safe(x):
    if hasattr(x,'item'): return x.item()
    if hasattr(x,'tolist'): return x.tolist()
    raise TypeError(type(x).__name__)


def prepare():
    if STATE.exists() and CHECK.exists(): return
    from ase.io import iread, write
    start=ROOT/'HEATING_CHECKPOINT_00800.extxyz'
    state0=json.loads((ROOT/'HEATING_CONTINUATION_STATE.json').read_text())
    CHECK.write_bytes(start.read_bytes())
    cfg={**state0,'stage':'heating','completed_cooling_steps':0,'melt_established':False,
         'stop_reason':'running','recent_rows':[],'resumed_wall_seconds':0.0}
    if OUTPUT.exists(): raise RuntimeError(f'{OUTPUT.name} exists without resume state; refusing to overwrite')
    for a in iread(ROOT/'HEAT_COOL_OVITO.extxyz',index=':'):
        step=int(a.info.get('step',-1))
        if step>100+int(state0['completed_heating_steps']): break
        write(OUTPUT,a,format='extxyz',append=OUTPUT.exists())
        if a.info.get('stage')=='heating': cfg['recent_rows'].append(dict(a.info))
    cfg['recent_rows']=cfg['recent_rows'][-3:]
    STATE.write_text(json.dumps(cfg,indent=2,default=json_safe)+'\n')
    print(f"RESUME_START heat_ps={cfg['heating_time_ps']}ps global_step={100+cfg['completed_heating_steps']} frames={len(cfg['recent_rows'])}",flush=True)


def worker(chunk_steps:int):
    import numpy as np
    from ase import units
    from ase.io import read, write
    import run_compact_mace as m
    cfg=json.loads(STATE.read_text())
    a=read(CHECK)
    rng=np.random.default_rng()
    rng.bit_generator.state=cfg['random_generator_state']
    m.attach(a)
    ref=read(ROOT/'BASELINE_CHECKPOINT.extxyz')
    refpos=ref.positions.copy()
    ag,tsc,hot,sink=m.masks(a)
    stage=cfg['stage']
    if stage=='heating':
        temperature=np.full((len(a),1),m.COLD_K);friction=np.zeros((len(a),1))
        temperature[hot]=m.HOT_K
        friction[hot]=m.FRICTION_FS/units.fs
        friction[sink]=m.FRICTION_FS/units.fs
        done=int(cfg['completed_heating_steps'])
        global_start=100+done
    elif stage=='cooling':
        temperature=np.full((len(a),1),m.COLD_K);friction=np.zeros((len(a),1))
        friction[sink]=m.FRICTION_FS/units.fs
        done=int(cfg['completed_cooling_steps'])
        global_start=100+int(cfg['completed_heating_steps'])+done
    else:
        raise RuntimeError(f'unknown stage {stage}')
    dyn=m.Langevin(a,timestep=m.DT_FS*units.fs,temperature_K=temperature,
                   friction=friction,fixcm=False,rng=rng)
    recent=list(cfg.get('recent_rows',[]))
    chunk_path=ROOT/f"RESUME_CHUNK_{stage}_{global_start:06d}.extxyz"
    if chunk_path.exists(): chunk_path.unlink()
    reason='running'
    chunk_done=0
    for n in range(1,chunk_steps+1):
        dyn.run(1);chunk_done=n
        if n%10: continue
        global_step=global_start+n
        row,f,fn,q,coord,disp,cand=m.metrics(a,refpos,global_step,stage)
        m.write_frame(chunk_path,a,row,f,fn,q,coord,disp,cand)
        recent.append(row);recent=recent[-3:]
        print(json.dumps({'stage':stage,'step':global_step,
              'heat_time_ps':(int(cfg['completed_heating_steps'])+(n if stage=='heating' else 0))*m.DT_FS/1000,
              'cool_time_ps':(int(cfg['completed_cooling_steps'])+(n if stage=='cooling' else 0))*m.DT_FS/1000,
              'T_hot_Ag_K':row['T_hot_Ag_K'],'Fmax_eV_A':row['Fmax_eV_A'],
              'hot_liquid_candidate_Ag':row['hot_liquid_candidate_Ag'],
              'hot_max_liquid_cluster_Ag':row['hot_max_liquid_cluster_Ag'],
              'outside_hot_max_liquid_cluster_Ag':row['outside_hot_max_cluster_Ag'],
              'q6_hot_median':row['q6_hot_median'],'hot_Ag_displacement_rms_A':row['hot_Ag_displacement_rms_A']},separators=(',',':'),default=json_safe),flush=True)
        if m.severe(row):
            reason='severe_force_or_short_contact';break
        if stage=='heating':
            if m.persistent_melt(recent,int(hot.sum())):
                cfg['melt_established']=True
                cfg['melt_start_heat_ps']=float((recent[-3]['step']-100)*m.DT_FS/1000)
                cfg['stage']='cooling';cfg['melt_reason']='persistent_localized_Ag_liquid_cluster';cfg['stop_reason']='running'
                reason='melt_established'
                break
            # Outside-hot candidate clusters remain metrics only per user directive;
            # they do not stop the run or count as proof of the target melt pool.
        else:
            if m.persistent_solidification(recent):
                cfg['resolidification_established']=True
                reason='resolidification_established';break
    if stage=='heating':
        cfg['completed_heating_steps']=int(cfg['completed_heating_steps'])+chunk_done
        cfg['heating_time_ps']=cfg['completed_heating_steps']*m.DT_FS/1000
    else:
        cfg['completed_cooling_steps']=int(cfg['completed_cooling_steps'])+chunk_done
        cfg['cooling_time_ps']=cfg['completed_cooling_steps']*m.DT_FS/1000
    if reason in {'severe_force_or_short_contact','remote_liquid_cluster_competes_with_single_pool'}:
        cfg['stop_reason']=reason
    if reason=='resolidification_established':
        cfg['resolidification_established']=True;cfg['stop_reason']=reason
    cfg['recent_rows']=recent
    cfg['random_generator_state']=rng.bit_generator.state
    nxt=ROOT/'HEATING_RESUME_CHECKPOINT_NEXT.extxyz'
    write(nxt,a,format='extxyz')
    tmp=ROOT/'HEATING_RESUME_STATE_NEXT.json'
    tmp.write_text(json.dumps(cfg,indent=2,default=json_safe)+'\n')
    print(f'CHUNK_DONE stage={stage} steps={chunk_done} reason={reason}',flush=True)
    return 0


def supervisor():
    prepare()
    started=time.monotonic()
    log=(ROOT/'HEATING_RESUME.log').open('a',buffering=1)
    env=os.environ.copy();env.update({'OMP_NUM_THREADS':'4','MKL_NUM_THREADS':'4','OPENBLAS_NUM_THREADS':'1'})
    while True:
        cfg=json.loads(STATE.read_text())
        if cfg.get('stop_reason') not in {'running','safety_limit'}:
            break
        stage=cfg['stage']
        if stage=='heating' and int(cfg['completed_heating_steps'])>=HEAT_MAX:
            cfg['stop_reason']='4ps_heat_limit';STATE.write_text(json.dumps(cfg,indent=2,default=json_safe)+'\n');break
        if stage=='cooling' and int(cfg['completed_cooling_steps'])>=COOL_MAX:
            cfg['stop_reason']='2ps_cooling_limit';STATE.write_text(json.dumps(cfg,indent=2,default=json_safe)+'\n');break
        n=CHUNK
        if stage=='heating': n=min(n,HEAT_MAX-int(cfg['completed_heating_steps']))
        else: n=min(n,COOL_MAX-int(cfg['completed_cooling_steps']))
        cmd=[sys.executable,str(Path(__file__).resolve()),'worker',str(n)]
        proc=subprocess.Popen(cmd,cwd=ROOT,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1)
        assert proc.stdout
        for line in proc.stdout:
            print(line.rstrip(),flush=True);log.write(line)
        code=proc.wait()
        if code!=0:
            cfg=json.loads(STATE.read_text());cfg['stop_reason']=f'worker_exit_{code}';cfg['last_error']='chunk worker exited; last committed checkpoint preserved'
            STATE.write_text(json.dumps(cfg,indent=2,default=json_safe)+'\n')
            print(f'WORKER_FAILED exit={code} checkpoint_retained={CHECK.name}',flush=True);break
        chunk_path=ROOT/f"RESUME_CHUNK_{stage}_{100+int(cfg['completed_heating_steps'])+int(cfg['completed_cooling_steps']):06d}.extxyz"
        # Worker names the chunk using the global step at chunk start.
        candidates=sorted(ROOT.glob(f'RESUME_CHUNK_{stage}_*.extxyz'),key=lambda p:p.stat().st_mtime)
        if candidates:
            chunk_path=candidates[-1]
            from ase.io import iread,write
            for a in iread(chunk_path,index=':'):
                write(OUTPUT,a,format='extxyz',append=True)
            chunk_path.unlink()
        os.replace(ROOT/'HEATING_RESUME_CHECKPOINT_NEXT.extxyz',CHECK)
        os.replace(ROOT/'HEATING_RESUME_STATE_NEXT.json',STATE)
        cfg=json.loads(STATE.read_text())
        cfg['resumed_wall_seconds']=time.monotonic()-started
        STATE.write_text(json.dumps(cfg,indent=2,default=json_safe)+'\n')
        print(f"CHECKPOINT stage={cfg['stage']} heat_ps={cfg.get('heating_time_ps')} cool_ps={cfg.get('cooling_time_ps',0)} reason={cfg.get('stop_reason')}",flush=True)
        if cfg.get('stop_reason') not in {'running','safety_limit'}: break
    log.close()

if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='worker': sys.exit(worker(int(sys.argv[2])))
    supervisor()
