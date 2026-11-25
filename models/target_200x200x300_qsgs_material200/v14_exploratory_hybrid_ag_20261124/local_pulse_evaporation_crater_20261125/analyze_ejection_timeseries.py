#!/usr/bin/env python3
"""Track surface-local Ag ejection during both narrow-pulse trajectories."""
from pathlib import Path
import csv
import numpy as np
ROOT=Path(__file__).resolve().parent
CASES={
 '474eV_sigma12_depth6':ROOT/'trajectory_surface_local_2ps/ag_hotspot_nve_surface_local.dump',
 '740eV_sigma12_depth6':ROOT/'trajectory_surface_local_33750_2ps/ag_hotspot_nve_surface_local_33750.dump',
}
fields=['case','step','elapsed_ps','N_Ag','N_z_gt_220','N_z_gt_250','N_source_r25_z_gt_220','N_source_r25_z_gt_250','max_z_A','max_source_r25_z_A','N_center_r20_z180_220']
rows=[]
seen={}
for label,path in CASES.items():
    with path.open() as f:
        while True:
            line=f.readline()
            if not line: break
            if line.strip()!='ITEM: TIMESTEP': continue
            step=int(f.readline())
            assert f.readline().strip()=='ITEM: NUMBER OF ATOMS'
            n=int(f.readline())
            f.readline()  # BOX BOUNDS
            bounds=[f.readline() for _ in range(3)]
            names=f.readline().split()[2:]
            idx={name:i for i,name in enumerate(names)}
            a=np.loadtxt(f,max_rows=n)
            x,y,z=(a[:,idx[k]] for k in ('x','y','z'))
            ids=a[:,idx['id']].astype(np.int64)
            r=np.hypot((x-100+100)%200-100,(y-100+100)%200-100)
            local=r<25
            high_local=local&(z>250)
            seen.setdefault(label,set()).update(ids[high_local].tolist())
            rows.append({'case':label,'step':step,'elapsed_ps':(step-6800)*0.00025,
              'N_Ag':n,'N_z_gt_220':int((z>220).sum()),'N_z_gt_250':int((z>250).sum()),
              'N_source_r25_z_gt_220':int((local&(z>220)).sum()),
              'N_source_r25_z_gt_250':int(high_local.sum()),'max_z_A':float(z.max()),
              'max_source_r25_z_A':float(z[local].max()),
              'N_center_r20_z180_220':int(((r<20)&(z>=180)&(z<220)).sum())})
        # two runs' frames are appended in order, no cross-file mixing
out=ROOT/'results/trajectory_metrics_surface_local_pulses.csv'
with out.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
for label,ids in seen.items():
    print(label,'frames=',sum(r['case']==label for r in rows),'unique source IDs ever above 250=',sorted(ids))
    rr=[r for r in rows if r['case']==label]
    print('max source z=',max(r['max_source_r25_z_A'] for r in rr), 'peak local counts >250=',max(r['N_source_r25_z_gt_250'] for r in rr))
