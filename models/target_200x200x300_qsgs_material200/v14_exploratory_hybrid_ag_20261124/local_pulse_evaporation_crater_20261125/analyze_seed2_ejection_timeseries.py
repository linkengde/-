#!/usr/bin/env python3
"""Track local/global Ag escape and center density through seed-2 trajectory."""
from pathlib import Path
import csv
import numpy as np

ROOT=Path(__file__).resolve().parent
SRC=ROOT/'trajectory_surface_local_33750_seed2_2ps/ag_hotspot_nve_surface_local_33750_seed2.dump'
OUT=ROOT/'results/trajectory_metrics_surface_local_33750_seed2_2ps.csv'
rows=[]; ever_source_over250=set(); ever_source_over220=set()
with SRC.open() as f:
    while True:
        line=f.readline()
        if not line: break
        if line.strip()!='ITEM: TIMESTEP': continue
        step=int(f.readline()); assert f.readline().strip()=='ITEM: NUMBER OF ATOMS'
        n=int(f.readline()); f.readline()
        for _ in range(3): f.readline()
        names=f.readline().split()[2:]; idx={name:i for i,name in enumerate(names)}
        a=np.loadtxt(f,max_rows=n)
        x,y,z=(a[:,idx[k]] for k in ('x','y','z'))
        ids=a[:,idx['id']].astype(np.int64)
        r=np.hypot((x-100+100)%200-100,(y-100+100)%200-100)
        local=r<25; source220=local&(z>220); source250=local&(z>250)
        ever_source_over220.update(ids[source220].tolist()); ever_source_over250.update(ids[source250].tolist())
        center=(r<20)&(z>=180)&(z<220)
        rows.append({'step':step,'elapsed_ps':(step-6800)*.00025,'N_Ag':n,
          'N_z_gt_220':int((z>220).sum()),'N_z_gt_250':int((z>250).sum()),
          'N_source_r25_z_gt_220':int(source220.sum()),'N_source_r25_z_gt_250':int(source250.sum()),
          'max_z_A':float(z.max()),'max_source_r25_z_A':float(z[local].max()),
          'N_center_r20_z180_220':int(center.sum()),'max_center_r20_z_A':float(z[r<20].max())})
fields=list(rows[0])
with OUT.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n'); w.writeheader(); w.writerows(rows)
print(f'frames={len(rows)} time={rows[0]["elapsed_ps"]:.3f}-{rows[-1]["elapsed_ps"]:.3f} ps')
print(f'peak_source_z={max(r["max_source_r25_z_A"] for r in rows):.3f} A; peak source r<25 counts >220/>250={max(r["N_source_r25_z_gt_220"] for r in rows)}/{max(r["N_source_r25_z_gt_250"] for r in rows)}')
print(f'unique source IDs ever >220/>250={len(ever_source_over220)}/{len(ever_source_over250)}: {sorted(ever_source_over250)}')
print(f'center r<20 count z180-220 start/end={rows[0]["N_center_r20_z180_220"]}/{rows[-1]["N_center_r20_z180_220"]}')
