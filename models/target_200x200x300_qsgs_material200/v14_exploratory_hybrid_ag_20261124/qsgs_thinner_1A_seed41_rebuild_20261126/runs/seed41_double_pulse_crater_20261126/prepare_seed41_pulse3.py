#!/usr/bin/env python3
"""Prepare a third same-dose localized Ag velocity pulse at 10.8 ps."""
from pathlib import Path
import json
import numpy as np

ROOT=Path(__file__).resolve().parent
SRC=ROOT/'seed41_afterpulse2_final_10p8ps_nomsd_all.dump'
STEP=40400
TIME_PS=10.8
CENTER=(100.0,100.0,180.09)
SIGMA_XY=12.0
DEPTH=6.0
TARGET_E=1479.2492958829243
SEED=20261143
MASS_AG=107.8682
KINETIC_FACTOR=1.03642696526e-4
KB=8.617333262145e-5
BOX=200.0

with SRC.open() as f:
    assert f.readline().strip()=='ITEM: TIMESTEP'
    source_step=int(f.readline().strip())
    f.readline(); n=int(f.readline().strip())
    f.readline(); bounds=[f.readline() for _ in range(3)]
    header=f.readline().strip().split()[2:]
    cols={k:i for i,k in enumerate(header)}
    required=('id','type','x','y','z','vx','vy','vz')
    assert all(k in cols for k in required)
    rows=np.empty((n,len(header)),dtype=np.float64)
    for i in range(n): rows[i]=np.fromstring(f.readline(),sep=' ')
assert source_step==STEP and n==469038
ag=rows[:,cols['type']].astype(int)==4
n_ag=int(ag.sum()); assert n_ag==110610
x,y,z=(rows[:,cols[k]] for k in ('x','y','z'))
dx=(x-CENTER[0]+BOX/2)%BOX-BOX/2
dy=(y-CENTER[1]+BOX/2)%BOX-BOX/2
wxy=np.exp(-0.5*(dx*dx+dy*dy)/SIGMA_XY**2)
wz=np.exp(-np.maximum(CENTER[2]-z,0.0)/DEPTH)
weights=np.clip(wxy*wz,0.0,1.0)[ag]
vcols=[cols[k] for k in ('vx','vy','vz')]
oldv=rows[ag][:,vcols].copy()
rng=np.random.default_rng(SEED)
nominal_dt=TARGET_E/(1.5*KB*weights.sum())
sigma_v=np.sqrt(KB*nominal_dt/(MASS_AG*KINETIC_FACTOR))
kick=rng.normal(size=(n_ag,3))*sigma_v*np.sqrt(weights[:,None])
kick-=kick.mean(axis=0,keepdims=True)
a=.5*MASS_AG*KINETIC_FACTOR*np.sum(kick*kick)
b=MASS_AG*KINETIC_FACTOR*np.sum(oldv*kick)
scale=(-b+np.sqrt(b*b+4*a*TARGET_E))/(2*a)
kick*=scale
newv=oldv+kick
delta_ke=.5*MASS_AG*KINETIC_FACTOR*np.sum(newv*newv-oldv*oldv)
assert abs(delta_ke-TARGET_E)<1e-7
rows[np.ix_(np.flatnonzero(ag),vcols)]=newv
out=ROOT/'seed41_pulse3_1479eV_start.dump'
with out.open('w') as f:
    f.write(f'ITEM: TIMESTEP\n{STEP}\nITEM: NUMBER OF ATOMS\n{n}\nITEM: BOX BOUNDS pp pp ff\n')
    f.writelines(bounds)
    f.write('ITEM: ATOMS id type x y z vx vy vz\n')
    for row in rows:
        f.write(f'{int(row[cols["id"]])} {int(row[cols["type"]])} '+' '.join(f'{row[j]:.15g}' for j in (cols['x'],cols['y'],cols['z'],*vcols))+'\n')
info={'source_dump':SRC.name,'source_time_ps':TIME_PS,'source_step':source_step,
      'total_atoms':n,'Ag_atoms':n_ag,'seed':SEED,'center_A':CENTER,
      'sigma_xy_A':SIGMA_XY,'attenuation_depth_A':DEPTH,
      'target_added_kinetic_energy_eV':TARGET_E,'actual_added_kinetic_energy_eV':float(delta_ke),
      'spatial_weight_sum':float(weights.sum()),'effective_weighted_atoms':float(weights.sum()**2/(weights@weights)),
      'nominal_peak_delta_temperature_K':float(nominal_dt),'kick_scale':float(scale),
      'pulse_kind':'third localized zero-net-momentum random velocity impulse; same dose/profile as pulse 2; no TTM'}
(ROOT/'seed41_pulse3_parameters.json').write_text(json.dumps(info,indent=2,ensure_ascii=False)+'\n')
print(json.dumps(info,indent=2,ensure_ascii=False))
