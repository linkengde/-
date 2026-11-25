#!/usr/bin/env python3
"""Compare 1.109 keV localized Ag pulse outcomes across two random seeds."""
from pathlib import Path
import csv, json
import numpy as np

ROOT=Path(__file__).resolve().parent
SEEDS=(20261131,20261132)
THERMO_FIELDS=['step','time_ps','TAg_K','Tspot_K','PE_eV','KE_eV','Etot_eV','Fmax_eV_A',
               'Q6mean','highQ6_fraction','MSD_A2','zmax_A','N_Ag_above220','N_Ag_above250','N_atoms']

def read_single_dump(path):
    with path.open() as f:
        while f.readline().strip()!='ITEM: NUMBER OF ATOMS': pass
        n=int(f.readline()); f.readline()
        bounds=np.asarray([[float(v) for v in f.readline().split()[:2]] for _ in range(3)])
        names=f.readline().split()[2:]; idx={k:i for i,k in enumerate(names)}
        atoms=np.loadtxt(f,max_rows=n)
    return idx,atoms,bounds

def spatial(path):
    idx,a,_=read_single_dump(path); ag=a[a[:,idx['type']]==4]
    x,y,z=(ag[:,idx[k]] for k in ('x','y','z'))
    r=np.hypot((x-100+100)%200-100,(y-100+100)%200-100)
    out={'N_Ag':int(len(ag)),'zmax_A':float(z.max()),'N_z_gt_220':int((z>220).sum()),'N_z_gt_250':int((z>250).sum())}
    for radius in (20,25,40):
        m=r<radius
        out[f'r_lt_{radius}_N_180_220']=int((m&(z>=180)&(z<220)).sum())
        out[f'r_lt_{radius}_N_170_220']=int((m&(z>=170)&(z<220)).sum())
        out[f'r_lt_{radius}_N_220_250']=int((m&(z>=220)&(z<250)).sum())
        out[f'r_lt_{radius}_N_gt_250']=int((m&(z>=250)).sum())
    for label,rlo,rhi in [('center',0,20),('annulus',20,40)]:
        heights=[]
        for x0 in range(0,200,10):
            for y0 in range(0,200,10):
                rr=np.hypot((x0+5-100+100)%200-100,(y0+5-100+100)%200-100)
                if not rlo<=rr<rhi: continue
                m=(x>=x0)&(x<x0+10)&(y>=y0)&(y<y0+10)&(z>=150)&(z<230)
                if m.sum()>=3: heights.append(float(np.quantile(z[m],.95)))
        out[label+'_surface_bin_count']=len(heights)
        if heights:
            out[label+'_surface_q50_A']=float(np.quantile(heights,.5))
            out[label+'_surface_q10_A']=float(np.quantile(heights,.1))
    return out

def read_thermo(path):
    rows=[]
    for line in path.read_text().splitlines():
        p=line.split()
        if len(p)!=len(THERMO_FIELDS) or not p[0].isdigit(): continue
        try:
            r=dict(zip(THERMO_FIELDS,map(float,p)))
            r['elapsed_since_pulse_ps']=(r['step']-6800)*.00025
            rows.append(r)
        except ValueError: pass
    if not rows: raise RuntimeError(f'No thermo rows parsed: {path}')
    et=[r['Etot_eV'] for r in rows]
    return rows,{'records':len(rows),'elapsed_ps':rows[-1]['elapsed_since_pulse_ps'],
      'total_energy_range_eV':max(et)-min(et),'final_Ag_temp_K':rows[-1]['TAg_K'],
      'final_hotspot_temp_K':rows[-1]['Tspot_K'],'final_Fmax_eV_per_A':rows[-1]['Fmax_eV_A'],
      'final_zmax_A':rows[-1]['zmax_A'],'final_Ag_above220':int(rows[-1]['N_Ag_above220']),
      'final_Ag_above250':int(rows[-1]['N_Ag_above250']),
      'all_atoms_constant':all(int(r['N_atoms'])==458248 for r in rows)}

def trajectory_metrics(path):
    rows=[]; ever220=set(); ever250=set()
    with path.open() as f:
        while True:
            line=f.readline()
            if not line: break
            if line.strip()!='ITEM: TIMESTEP': continue
            step=int(f.readline()); assert f.readline().strip()=='ITEM: NUMBER OF ATOMS'
            n=int(f.readline()); f.readline()
            for _ in range(3): f.readline()
            names=f.readline().split()[2:]; idx={k:i for i,k in enumerate(names)}
            a=np.loadtxt(f,max_rows=n)
            x,y,z=(a[:,idx[k]] for k in ('x','y','z')); ids=a[:,idx['id']].astype(np.int64)
            r=np.hypot((x-100+100)%200-100,(y-100+100)%200-100)
            local=r<25; m220=local&(z>220); m250=local&(z>250)
            ever220.update(ids[m220].tolist()); ever250.update(ids[m250].tolist())
            mcenter=(r<20)&(z>=180)&(z<220)
            rows.append({'seed':path.stem.split('_seed')[-1].split('.')[0],'step':step,
              'elapsed_ps':(step-6800)*.00025,'N_Ag':n,'N_z_gt_220':int((z>220).sum()),
              'N_z_gt_250':int((z>250).sum()),'N_source_r25_z_gt_220':int(m220.sum()),
              'N_source_r25_z_gt_250':int(m250.sum()),'max_z_A':float(z.max()),
              'max_source_r25_z_A':float(z[local].max()),'N_center_r20_z180_220':int(mcenter.sum())})
    if not rows: raise RuntimeError(f'No frames parsed: {path}')
    return rows,{'frames':len(rows),'peak_source_z_A':max(r['max_source_r25_z_A'] for r in rows),
      'peak_source_r25_N_z_gt_220':max(r['N_source_r25_z_gt_220'] for r in rows),
      'peak_source_r25_N_z_gt_250':max(r['N_source_r25_z_gt_250'] for r in rows),
      'unique_source_ids_ever_z_gt_220':len(ever220),'unique_source_ids_ever_z_gt_250':len(ever250),
      'center_r20_N_180_220_start_end':[rows[0]['N_center_r20_z180_220'],rows[-1]['N_center_r20_z180_220']]}

all_traj=[]; results={}
for seed in SEEDS:
    tag=f'1109eV_seed{seed}'
    start=spatial(ROOT/'prepared'/f'local_gaussian_pulse_1109eV_seed{seed}_start.dump')
    final=spatial(ROOT/f'trajectory_surface_local_1109eV_seed{seed}_2ps'/f'final_all_atoms_1109eV_seed{seed}.dump')
    trows,thermo=read_thermo(ROOT/f'log.local_gaussian_pulse_1109eV_seed{seed}_2ps.lammps')
    traj_path=ROOT/f'trajectory_surface_local_1109eV_seed{seed}_2ps'/f'ag_hotspot_1109eV_seed{seed}.dump'
    tr,track=trajectory_metrics(traj_path); all_traj.extend(tr)
    results[tag]={'start':start,'final':final,'thermo':thermo,'trajectory':track}

OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)
(OUT/'1109eV_two_seed_analysis.json').write_text(json.dumps(results,indent=2)+'\n')
with (OUT/'thermo_timeseries_1109eV_two_seeds.csv').open('w',newline='') as f:
    fields=['seed']+['time_since_pulse_ps']+THERMO_FIELDS
    w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n'); w.writeheader()
    for seed in SEEDS:
        rs,_=read_thermo(ROOT/f'log.local_gaussian_pulse_1109eV_seed{seed}_2ps.lammps')
        for r in rs:
            row={k:v for k,v in r.items() if k!='elapsed_since_pulse_ps'}
            w.writerow({'seed':seed,'time_since_pulse_ps':r['elapsed_since_pulse_ps'],**row})
with (OUT/'trajectory_metrics_1109eV_two_seeds.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(all_traj[0]),lineterminator='\n'); w.writeheader(); w.writerows(all_traj)
print(json.dumps(results,indent=2))
