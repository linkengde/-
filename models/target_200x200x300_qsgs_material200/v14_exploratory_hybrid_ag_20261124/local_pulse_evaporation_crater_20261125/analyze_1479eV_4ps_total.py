#!/usr/bin/env python3
"""Summarize the 2 ps continuation runs and compare crater metrics."""
from pathlib import Path
import csv, json
import numpy as np

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)
SEEDS=(20261133,20261134)
THERMO_FIELDS=['step','time_ps','TAg_K','Tspot_K','PE_eV','KE_eV','Etot_eV','Fmax_eV_A',
               'Q6mean','highQ6_fraction','MSD_A2','zmax_A','N_Ag_above220','N_Ag_above250','N_atoms']

def frames(path):
    with path.open() as f:
        while True:
            line=f.readline()
            if not line: break
            if line.strip()!='ITEM: TIMESTEP': continue
            step=int(f.readline()); assert f.readline().strip()=='ITEM: NUMBER OF ATOMS'
            n=int(f.readline()); f.readline()
            bounds=np.asarray([[float(v) for v in f.readline().split()[:2]] for _ in range(3)])
            names=f.readline().split()[2:]; idx={k:i for i,k in enumerate(names)}
            a=np.loadtxt(f,max_rows=n)
            yield step,idx,a,bounds

def surface_bins(x,y,z,rlo,rhi):
    heights=[]
    for x0 in range(0,200,10):
        for y0 in range(0,200,10):
            rr=np.hypot((x0+5-100+100)%200-100,(y0+5-100+100)%200-100)
            if not rlo<=rr<rhi: continue
            m=(x>=x0)&(x<x0+10)&(y>=y0)&(y<y0+10)&(z>=150)&(z<230)
            if m.sum()>=3: heights.append(float(np.quantile(z[m],.95)))
    return heights

def spatial(path):
    _,idx,a,_=next(frames(path)); ag=a[a[:,idx['type']]==4]
    x,y,z=(ag[:,idx[k]] for k in ('x','y','z'))
    r=np.hypot((x-100+100)%200-100,(y-100+100)%200-100)
    h0=surface_bins(x,y,z,0,20); h1=surface_bins(x,y,z,20,40)
    c=((r<20)&(z>=180)&(z<220)); a2040=((r>=20)&(r<40)&(z>=180)&(z<220))
    return {'N_Ag':int(len(ag)),'zmax_A':float(z.max()),'N_Ag_above220':int((z>220).sum()),
      'N_Ag_above250':int((z>250).sum()),'N_center_r20_z180_220':int(c.sum()),
      'N_annulus_r20_40_z180_220':int(a2040.sum()),
      'center_density_ratio_to_annulus':float((c.sum()/400)/(a2040.sum()/1200)) if a2040.sum() else None,
      'center_surface_bins':len(h0),'annulus_surface_bins':len(h1),
      'center_surface_q50_A':float(np.quantile(h0,.5)) if h0 else None,
      'annulus_surface_q50_A':float(np.quantile(h1,.5)) if h1 else None,
      'center_minus_annulus_q50_A':float(np.quantile(h0,.5)-np.quantile(h1,.5)) if h0 and h1 else None}

def read_thermo(path):
    out=[]
    for line in path.read_text().splitlines():
        p=line.split()
        if len(p)!=len(THERMO_FIELDS) or not p[0].isdigit(): continue
        try:
            r=dict(zip(THERMO_FIELDS,map(float,p)))
            r['elapsed_since_pulse_ps']=(r['step']-6800)*0.00025
            out.append(r)
        except ValueError: pass
    if not out: raise RuntimeError(f'No thermo rows in {path}')
    return out

all_rows=[]; report={}
previous=json.loads((OUT/'1479eV_two_seed_analysis.json').read_text())
for seed in SEEDS:
    name=f'1479eV_seed{seed}'
    ext_dump=ROOT/f'trajectory_surface_local_1479eV_seed{seed}_extension2ps/ag_extension_1479eV_seed{seed}.dump'
    ext_rows=[]; ever220=set(); ever250=set()
    for step,idx,a,_ in frames(ext_dump):
        x,y,z=(a[:,idx[k]] for k in ('x','y','z')); ids=a[:,idx['id']].astype(np.int64)
        r=np.hypot((x-100+100)%200-100,(y-100+100)%200-100)
        src=r<25; c=(r<20)&(z>=180)&(z<220); ann=(r>=20)&(r<40)&(z>=180)&(z<220)
        m220=src&(z>220); m250=src&(z>250)
        ever220.update(ids[m220].tolist()); ever250.update(ids[m250].tolist())
        hc=surface_bins(x,y,z,0,20); ha=surface_bins(x,y,z,20,40)
        cq=float(np.quantile(hc,.5)) if hc else float('nan'); aq=float(np.quantile(ha,.5)) if ha else float('nan')
        row={'seed':seed,'step':step,'elapsed_since_pulse_ps':(step-6800)*.00025,
          'N_Ag':int(len(a)),'N_Ag_above220':int((z>220).sum()),'N_Ag_above250':int((z>250).sum()),
          'N_source_r25_above220':int(m220.sum()),'N_source_r25_above250':int(m250.sum()),
          'unique_source_ids_ever_above220':len(ever220),'unique_source_ids_ever_above250':len(ever250),
          'N_center_r20_z180_220':int(c.sum()),'N_annulus_r20_40_z180_220':int(ann.sum()),
          'center_density_ratio_to_annulus':float((c.sum()/400)/(ann.sum()/1200)) if ann.sum() else float('nan'),
          'center_surface_bins':len(hc),'annulus_surface_bins':len(ha),
          'center_surface_q50_A':cq,'annulus_surface_q50_A':aq,'center_minus_annulus_q50_A':cq-aq,
          'max_z_A':float(z.max()),'max_source_r25_z_A':float(z[src].max())}
        ext_rows.append(row); all_rows.append(row)
    thermo=read_thermo(ROOT/f'log.continue_1479eV_seed{seed}_2ps.lammps')
    prior=read_thermo(ROOT/f'log.local_gaussian_pulse_1479eV_seed{seed}_2ps.lammps')
    allthermo=prior+thermo
    et=[r['Etot_eV'] for r in allthermo]
    ext_et=[r['Etot_eV'] for r in thermo]
    last_half=[r for r in ext_rows if r['elapsed_since_pulse_ps']>=3.5]
    report[name]={
      'two_ps':previous[name]['final'],
      'four_ps':spatial(ROOT/f'trajectory_surface_local_1479eV_seed{seed}_extension2ps/final_all_atoms_1479eV_seed{seed}_4ps.dump'),
      'continuation_thermo':{'records':len(thermo),'first_step':thermo[0]['step'],'last_step':thermo[-1]['step'],
        'total_energy_range_eV_0_to_4ps':max(et)-min(et),'final_Ag_temp_K':thermo[-1]['TAg_K'],
        'continuation_only_total_energy_range_eV':max(ext_et)-min(ext_et),
        'restart_readback_energy_offset_eV':thermo[0]['Etot_eV']-prior[-1]['Etot_eV'],
        'final_hotspot_temp_K':thermo[-1]['Tspot_K'],'final_Fmax_eV_A':thermo[-1]['Fmax_eV_A'],
        'all_atoms_constant':all(int(r['N_atoms'])==458248 for r in allthermo)},
      'extension_trajectory':{'frames':len(ext_rows),'unique_local_ids_ever_above220':len(ever220),
        'unique_local_ids_ever_above250':len(ever250),
        'peak_local_r25_count_above220':max(r['N_source_r25_above220'] for r in ext_rows),
        'peak_local_r25_count_above250':max(r['N_source_r25_above250'] for r in ext_rows),
        'max_local_r25_z_A':max(r['max_source_r25_z_A'] for r in ext_rows),
        'last_0p5ps_center_minus_annulus_q50_median_A':float(np.median([r['center_minus_annulus_q50_A'] for r in last_half])),
        'last_0p5ps_center_minus_annulus_q50_min_A':float(min(r['center_minus_annulus_q50_A'] for r in last_half)),
        'last_0p5ps_center_minus_annulus_q50_max_A':float(max(r['center_minus_annulus_q50_A'] for r in last_half)),
        'last_0p5ps_center_density_ratio_median':float(np.median([r['center_density_ratio_to_annulus'] for r in last_half]))}}

(OUT/'1479eV_4ps_total_analysis.json').write_text(json.dumps(report,indent=2)+'\n')
with (OUT/'trajectory_metrics_1479eV_extension_2ps.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(all_rows[0]),lineterminator='\n'); w.writeheader(); w.writerows(all_rows)
with (OUT/'thermo_timeseries_1479eV_extension_2ps.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['seed','elapsed_since_pulse_ps']+THERMO_FIELDS,lineterminator='\n'); w.writeheader()
    for seed in SEEDS:
        for row in read_thermo(ROOT/f'log.continue_1479eV_seed{seed}_2ps.lammps'):
            w.writerow({'seed':seed,'elapsed_since_pulse_ps':row['elapsed_since_pulse_ps'],**row})
print(json.dumps(report,indent=2))
