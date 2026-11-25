#!/usr/bin/env python3
"""Analyze the independent seed-2 replicate of the localized Ag pulse."""
from pathlib import Path
import csv, json
import numpy as np

ROOT = Path(__file__).resolve().parent
START = ROOT / 'prepared/local_gaussian_pulse_surface_local_33750_seed2_start.dump'
FINAL = ROOT / 'trajectory_surface_local_33750_seed2_2ps/final_all_atoms_surface_local_33750_seed2.dump'
LOG = ROOT / 'log.local_gaussian_pulse_surface_local_33750_seed2_2ps.lammps'

def read_dump(path):
    with path.open() as f:
        while f.readline().strip() != 'ITEM: NUMBER OF ATOMS':
            pass
        n = int(f.readline())
        f.readline()
        bounds = np.array([[float(v) for v in f.readline().split()[:2]] for _ in range(3)])
        names = f.readline().split()[2:]
        atoms = np.loadtxt(f, max_rows=n)
    return {name:i for i,name in enumerate(names)}, atoms, bounds

def spatial_metrics(path):
    idx, atoms, _ = read_dump(path)
    ag = atoms[atoms[:,idx['type']] == 4]
    x,y,z = (ag[:,idx[k]] for k in ('x','y','z'))
    r = np.hypot((x-100+100)%200-100,(y-100+100)%200-100)
    out={'N_Ag':int(len(ag)),'zmax_A':float(z.max()),'N_z_gt_220':int((z>220).sum()),
         'N_z_gt_250':int((z>250).sum())}
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

start,final=spatial_metrics(START),spatial_metrics(FINAL)
idx,atoms,_=read_dump(FINAL)
ag=atoms[atoms[:,idx['type']]==4]
high=[]
for atom in ag[ag[:,idx['z']]>250]:
    dx=(atom[idx['x']]-100+100)%200-100; dy=(atom[idx['y']]-100+100)%200-100
    high.append({'id':int(atom[idx['id']]),'x_A':float(atom[idx['x']]),'y_A':float(atom[idx['y']]),
                 'z_A':float(atom[idx['z']]),'r_from_source_A':float(np.hypot(dx,dy))})
fields=['step','lammps_time_ps','TAg_K','Tspot_K','PE_eV','KE_eV','Etot_eV','Fmax_eV_A',
        'Q6mean','highQ6_fraction','MSD_A2','zmax_A','N_Ag_above220','N_Ag_above250','N_atoms']
rows=[]
for line in LOG.read_text().splitlines():
    p=line.split()
    if len(p)==len(fields) and p[0].isdigit():
        try:
            row=dict(zip(fields,map(float,p)))
            row['time_since_pulse_ps']=(row['step']-6800)*.00025
            rows.append(row)
        except ValueError: pass
thermo={}
if rows:
    et=[r['Etot_eV'] for r in rows]
    thermo={'records':len(rows),'start_step':int(rows[0]['step']),'end_step':int(rows[-1]['step']),
            'elapsed_ps':rows[-1]['time_since_pulse_ps'],'total_energy_range_eV':max(et)-min(et),
            'final_Ag_temp_K':rows[-1]['TAg_K'],'final_hotspot_temp_K':rows[-1]['Tspot_K'],
            'final_Fmax_Ag_eV_per_A':rows[-1]['Fmax_eV_A'],'final_zmax_A':rows[-1]['zmax_A'],
            'final_Ag_above220':int(rows[-1]['N_Ag_above220']),
            'final_Ag_above250':int(rows[-1]['N_Ag_above250']),
            'all_atoms_constant':all(int(r['N_atoms'])==458248 for r in rows)}
    csv_path=ROOT/'results/thermo_timeseries_surface_local_33750_seed2_2ps.csv'
    with csv_path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['time_since_pulse_ps']+fields,lineterminator='\n'); w.writeheader(); w.writerows(rows)
    (ROOT/'results/thermo_audit_surface_local_33750_seed2_2ps.json').write_text(json.dumps(thermo,indent=2)+'\n')
result={'start':start,'final':final,'high_z_atoms':high,'thermo':thermo}
(ROOT/'results/surface_local_33750_seed2_2ps_analysis.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
