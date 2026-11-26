import json, math, csv
from pathlib import Path
import numpy as np
ROOT=Path(__file__).parent

def read_dump(path,wanted_type=4):
    frames=[]
    with open(path) as f:
        while True:
            line=f.readline()
            if not line: break
            if not line.startswith('ITEM: TIMESTEP'): continue
            step=int(f.readline().strip()); f.readline(); n=int(f.readline().strip())
            f.readline(); bounds=[tuple(map(float,f.readline().split()[:2])) for _ in range(3)]
            atom_header=f.readline().strip().split()[2:]; cols={c:i for i,c in enumerate(atom_header)}
            data=np.empty((n,len(atom_header)),dtype=float)
            for i in range(n): data[i]=np.fromstring(f.readline(),sep=' ')
            if 'type' in cols: data=data[data[:,cols['type']]==wanted_type]
            frames.append((step,cols,data))
    return frames

def stats(frame,threshold=230.,zmin=150.,binw=4.):
    step,c,d=frame; x,y,z=(d[:,c[k]] for k in ('x','y','z'))
    xr=x-100.; yr=y-100.; r=np.sqrt(xr*xr+yr*yr); active=(z>=zmin)&(z<=threshold)
    xb=np.floor((xr+50.)/binw).astype(int); yb=np.floor((yr+50.)/binw).astype(int); bins={}
    for i in np.where(active)[0]:
        key=(xb[i],yb[i]); bins[key]=max(bins.get(key,-1e9),z[i])
    center=[]; ring=[]
    for (ix,iy),h in bins.items():
        rr=math.hypot((ix+.5)*binw-50,(iy+.5)*binw-50)
        if rr<10: center.append(h)
        elif 20<=rr<40: ring.append(h)
    def q(a,p): return float(np.percentile(a,p)) if len(a) else None
    cm,rm=q(center,50),q(ring,50)
    return {'step':step,'time_ps':step*.00025,'max_z':float(np.max(z)),
      'n_z_gt_220':int(np.count_nonzero(z>220)),'n_z_gt_250':int(np.count_nonzero(z>250)),
      'n_z_gt_300':int(np.count_nonzero(z>300)),
      'n_local_r25_z_gt_220':int(np.count_nonzero((r<25)&(z>220))),
      'n_local_r25_z_gt_250':int(np.count_nonzero((r<25)&(z>250))),
      'n_local_r25_z_gt_300':int(np.count_nonzero((r<25)&(z>300))),
      'n_local_r20_z_150_220':int(np.count_nonzero((r<20)&(z>=150)&(z<=220))),
      'center_bins':len(center),'ring_bins':len(ring),'center_top_median_A':cm,
      'center_top_q90_A':q(center,90),'ring_top_median_A':rm,'ring_top_q90_A':q(ring,90),
      'center_minus_ring_A':cm-rm if cm is not None and rm is not None else None}

traj0=read_dump(ROOT/'seed41_local_pulse_1479eV_ag.dump')
traj1=read_dump(ROOT/'seed41_local_pulse_continue_ag.dump')
traj_by_step={fr[0]:fr for fr in traj0+traj1}
traj=[traj_by_step[k] for k in sorted(traj_by_step)]
start=read_dump(ROOT/'seed41_local_pulse_1479eV_start.dump')
final=read_dump(ROOT/'seed41_local_pulse_final_5p8ps_all.dump')
rows=[stats(fr) for fr in traj]
for row in rows: row['time_ps']=1.8+(row['step']-4400)*.00025
with open(ROOT/'seed41_local_pulse_frame_metrics.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
start_stats=stats(start[0]); start_stats['time_ps']=1.8
final_stats=stats(final[0]); final_stats['time_ps']=5.8

def unique_above(cut):
    ids=set()
    for _,c,d in traj:
        z=d[:,c['z']]
        for atomid in d[z>cut,c['id']].astype(int): ids.add(int(atomid))
    return len(ids)
thermo=[]
for logfile in [ROOT/'log.seed41_local_pulse_1479eV_2ps.lammps', ROOT/'log.seed41_local_pulse_continue_3p8_to5p8ps.lammps']:
 for line in logfile.read_text().splitlines():
  sp=line.split()
  if len(sp)>=19 and sp[0].isdigit():
   try:
    thermo.append({'step':int(sp[0]),'time_ps':float(sp[1]),'TAg_K':float(sp[2]),'Tspot_K':float(sp[3]),'PE_eV':float(sp[4]),'KE_eV':float(sp[5]),'Etot_eV':float(sp[6]),'Fmax_eV_A':float(sp[7]),'Q6mean':float(sp[8]),'Q6gt03':float(sp[9]),'MSD_A2':float(sp[10]),'zmax_A':float(sp[11]),'atoms':int(float(sp[-1]))})
   except Exception: pass
# Remove duplicate step at the restart boundary.
thermo=list({r['step']:r for r in thermo}.values()); thermo.sort(key=lambda r:r['step'])
summary={'run':'seed41 local velocity impulse, NVE follow-up; exploratory only','time_interval_ps':[1.8,5.8],
 'pulse_energy_eV':1479.2492958829243,'timestep_fs':.25,'total_atom_count_initial':469040,
 'Ag_count':int(len(traj[0][2])),'atom_loss_observed':False,'start_surface':start_stats,'final_surface':final_stats,
 'relative_surface_change':{'center_median_A':final_stats['center_top_median_A']-start_stats['center_top_median_A'],
 'ring_median_A':final_stats['ring_top_median_A']-start_stats['ring_top_median_A'],
 'center_minus_ring_change_A':final_stats['center_minus_ring_A']-start_stats['center_minus_ring_A'] if final_stats['center_minus_ring_A'] is not None and start_stats['center_minus_ring_A'] is not None else None,
 'local_ag_count_r20_z150_220_change':final_stats['n_local_r20_z_150_220']-start_stats['n_local_r20_z_150_220']},
 'sampled_unique_ejected_ids':{str(cut):unique_above(cut) for cut in (220,250,300)},
 'final_above_threshold':{str(cut):final_stats[f'n_z_gt_{cut}'] for cut in (220,250,300)},
 'final_local_above_threshold':{str(cut):final_stats[f'n_local_r25_z_gt_{cut}'] for cut in (220,250,300)},
 'thermo_first':thermo[0] if thermo else None,'thermo_last':thermo[-1] if thermo else None,
 'etotal_drift_eV':thermo[-1]['Etot_eV']-thermo[0]['Etot_eV'] if thermo else None,
 'thermo_records':len(thermo),'frames_sampled':len(traj),
 'caveat':'Atoms above z thresholds remain inside the 500 A cell; these are ejection candidates, not confirmed permanent evaporation. The backbone is frozen and Ag-backbone uses repulsive LJ/WCA; this is exploratory, not a physically validated interface/TTM model.'}
(ROOT/'seed41_local_pulse_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
print(json.dumps(summary,indent=2,ensure_ascii=False))
