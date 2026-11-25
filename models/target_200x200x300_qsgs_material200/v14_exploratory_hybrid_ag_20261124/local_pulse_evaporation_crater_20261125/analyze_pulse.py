from pathlib import Path
import json
import numpy as np

ROOT=Path(__file__).resolve().parent

def read_dump(path):
    with open(path) as f:
        while True:
            line=f.readline()
            if not line: raise RuntimeError(f'No frame in {path}')
            if line.startswith('ITEM: NUMBER OF ATOMS'):
                n=int(f.readline()); break
        bounds=[]
        assert f.readline().startswith('ITEM: BOX BOUNDS')
        for _ in range(3): bounds.append([float(x) for x in f.readline().split()[:2]])
        header=f.readline().split()[2:]
        arr=np.loadtxt(f,max_rows=n)
    return header,arr,np.asarray(bounds)

def summarize(path):
    header,arr,bounds=read_dump(path)
    ix={k:i for i,k in enumerate(header)}
    ag=arr[arr[:,ix['type']]==4]
    x,y,z=[ag[:,ix[k]] for k in ('x','y','z')]
    r=np.hypot((x-100+100)%200-100,(y-100+100)%200-100)
    out={'file':str(path.name),'N_Ag':int(len(ag)),'zmax_A':float(z.max()),'zmin_A':float(z.min()),'above_220_Ag':int((z>220).sum()),'above_230_Ag':int((z>230).sum()),'above_250_Ag':int((z>250).sum())}
    for rad in [10,15,20,25,35,50]:
        m=r<rad
        out[f'r_lt_{rad}_z_180_200']=int((m&(z>=180)&(z<200)).sum())
        out[f'r_lt_{rad}_z_200_220']=int((m&(z>=200)&(z<220)).sum())
        out[f'r_lt_{rad}_z_220_250']=int((m&(z>=220)&(z<250)).sum())
        out[f'r_lt_{rad}_z_gt250']=int((m&(z>=250)).sum())
        out[f'r_lt_{rad}_z_180_220']=int((m&(z>=180)&(z<220)).sum())
    # Average surface heights: 95th and 99th percentile by 10 A xy bins, provided bin has Ag.
    for rad,label in [(20,'center_r20'),(50,'field_r50')]:
        heights=[]
        for x0 in range(0,200,10):
            for y0 in range(0,200,10):
                xc=x0+5; yc=y0+5
                rr=np.hypot((xc-100+100)%200-100,(yc-100+100)%200-100)
                if rr>=rad: continue
                m=(x>=x0)&(x<x0+10)&(y>=y0)&(y<y0+10)
                zz=z[m]
                zz=zz[zz>150]
                if len(zz)>=3: heights.append(float(np.quantile(zz,0.95)))
        out[label+'_n_bins']=len(heights)
        out[label+'_surface_q50_A']=float(np.quantile(heights,.5)) if heights else None
        out[label+'_surface_q95_A']=float(np.quantile(heights,.95)) if heights else None
    # local count in the pulse source cylinder and same-height annular comparison
    for rlo,rhi,label in [(0,20,'center'),(20,40,'annulus')]:
        m=(r>=rlo)&(r<rhi)
        out[label+'_N_z180_220']=int((m&(z>=180)&(z<220)).sum())
        out[label+'_N_z170_220']=int((m&(z>=170)&(z<220)).sum())
    return out

paths=[ROOT.parent/'simple_ag_fixed_skeleton_20261125/reheat_melt_assessment_20261125/trajectories/ag_reheat_q6_ptm_msd_final.dump',ROOT/'trajectory_1ps/final_all_atoms_1ps.dump',ROOT/'trajectory_9000K/final_all_atoms_9000K.dump',ROOT/'trajectory_22500K/final_all_atoms_22500K.dump']
metrics={p.name:summarize(p) for p in paths}
(ROOT/'results/summary_density_metrics.json').write_text(json.dumps(metrics,indent=2,ensure_ascii=False)+'\n')
# Record LAMMPS thermo rows with elapsed time reconstructed from the step counter.
import re, csv
names=['step','lammps_time_ps','TAg_K','Tspot_K','PE_eV','KE_eV','Etot_eV','Fmax_eV_A','Q6mean','highQ6_fraction','MSD_A2','zmax_A','Ag_above220','Ag_above250','N_atoms']
rows=[]
for line in (ROOT/'log.local_gaussian_pulse_1ps_nve.lammps').read_text().splitlines():
    parts=line.split()
    if len(parts)==len(names) and parts[0].isdigit():
        try:
            row=dict(zip(names,[float(v) for v in parts]))
            row['time_since_pulse_ps']=(row['step']-6800)*0.00025
            rows.append(row)
        except ValueError:
            pass
if rows:
    with (ROOT/'results/thermo_timeseries_1ps.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['time_since_pulse_ps']+names)
        w.writeheader(); w.writerows(rows)
    et=[r['Etot_eV'] for r in rows]
    audit={'thermo_records':len(rows),'start_step':int(rows[0]['step']),'end_step':int(rows[-1]['step']),'elapsed_ps':rows[-1]['time_since_pulse_ps'],'total_energy_range_eV':max(et)-min(et),'start_total_energy_eV':et[0],'end_total_energy_eV':et[-1],'final_Ag_temp_K':rows[-1]['TAg_K'],'final_hotspot_temp_K':rows[-1]['Tspot_K'],'final_Fmax_Ag_eV_per_A':rows[-1]['Fmax_eV_A'],'final_zmax_A':rows[-1]['zmax_A'],'final_Ag_above220':int(rows[-1]['Ag_above220']),'final_Ag_above250':int(rows[-1]['Ag_above250']),'all_atoms_constant':all(int(r['N_atoms'])==458248 for r in rows)}
    (ROOT/'results/thermo_audit_1ps.json').write_text(json.dumps(audit,indent=2)+'\n')
    print('THERMO AUDIT')
    print(json.dumps(audit,indent=2))
# A separately calibrated time origin is used in the stronger pulse input.
rows2=[]
for line in (ROOT/'log.local_gaussian_pulse_9000K_nve.lammps').read_text().splitlines():
    parts=line.split()
    if len(parts)==len(names) and parts[0].isdigit():
        try:
            row=dict(zip(names,[float(v) for v in parts])); row['time_since_pulse_ps']=(row['step']-6800)*0.00025; rows2.append(row)
        except ValueError: pass
if rows2:
    with (ROOT/'results/thermo_timeseries_9000K.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['time_since_pulse_ps']+names); w.writeheader(); w.writerows(rows2)
    et=[r['Etot_eV'] for r in rows2]
    audit2={'thermo_records':len(rows2),'start_step':int(rows2[0]['step']),'end_step':int(rows2[-1]['step']),'elapsed_ps':rows2[-1]['time_since_pulse_ps'],'total_energy_range_eV':max(et)-min(et),'start_total_energy_eV':et[0],'end_total_energy_eV':et[-1],'final_Ag_temp_K':rows2[-1]['TAg_K'],'final_hotspot_temp_K':rows2[-1]['Tspot_K'],'final_Fmax_Ag_eV_per_A':rows2[-1]['Fmax_eV_A'],'final_zmax_A':rows2[-1]['zmax_A'],'final_Ag_above220':int(rows2[-1]['Ag_above220']),'final_Ag_above250':int(rows2[-1]['Ag_above250']),'all_atoms_constant':all(int(r['N_atoms'])==458248 for r in rows2)}
    (ROOT/'results/thermo_audit_9000K.json').write_text(json.dumps(audit2,indent=2)+'\n')
    print('9000K THERMO AUDIT'); print(json.dumps(audit2,indent=2))

rows3=[]
for line in (ROOT/'log.local_gaussian_pulse_22500K_nve.lammps').read_text().splitlines():
    parts=line.split()
    if len(parts)==len(names) and parts[0].isdigit():
        try:
            row=dict(zip(names,[float(v) for v in parts])); row['time_since_pulse_ps']=(row['step']-6800)*0.00025; rows3.append(row)
        except ValueError: pass
if rows3:
    with (ROOT/'results/thermo_timeseries_22500K.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['time_since_pulse_ps']+names); w.writeheader(); w.writerows(rows3)
    et=[r['Etot_eV'] for r in rows3]
    audit3={'thermo_records':len(rows3),'start_step':int(rows3[0]['step']),'end_step':int(rows3[-1]['step']),'elapsed_ps':rows3[-1]['time_since_pulse_ps'],'total_energy_range_eV':max(et)-min(et),'start_total_energy_eV':et[0],'end_total_energy_eV':et[-1],'final_Ag_temp_K':rows3[-1]['TAg_K'],'final_hotspot_temp_K':rows3[-1]['Tspot_K'],'final_Fmax_Ag_eV_per_A':rows3[-1]['Fmax_eV_A'],'final_zmax_A':rows3[-1]['zmax_A'],'final_Ag_above220':int(rows3[-1]['Ag_above220']),'final_Ag_above250':int(rows3[-1]['Ag_above250']),'all_atoms_constant':all(int(r['N_atoms'])==458248 for r in rows3)}
    (ROOT/'results/thermo_audit_22500K.json').write_text(json.dumps(audit3,indent=2)+'\n')
    print('22500K THERMO AUDIT'); print(json.dumps(audit3,indent=2))
print('DENSITY METRICS'); print(json.dumps(metrics,indent=2,ensure_ascii=False))
