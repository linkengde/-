#!/usr/bin/env python3
"""Apply a second localized Ag velocity pulse to the 5.8 ps seed41 state."""
from pathlib import Path
import json
import numpy as np
ROOT=Path(__file__).resolve().parent
SRC=ROOT.parent/'seed41_local_pulse_crater_20261126/seed41_local_pulse_final_5p8ps_all.dump'
CENTER_X=CENTER_Y=100.0
SURFACE_Z=180.09
SIGMA_XY=12.0
DEPTH=6.0
TARGET_E=1479.2492958829243
SEED=20261142
N_EXPECTED=469040
N_AG_EXPECTED=110612
MASS_AG=107.8682
KB=8.617333262145e-5
KINETIC_FACTOR=1.03642696526e-4
BOX_XY=200.0
STEP=20400
TIME_PS=5.8
lines=SRC.read_text().splitlines()
assert lines[0]=='ITEM: TIMESTEP' and int(lines[1])==STEP
n=int(lines[3]); assert n==N_EXPECTED
header=next(i for i,line in enumerate(lines) if line.startswith('ITEM: ATOMS '))
cols=lines[header].split()[2:]; idx={k:i for i,k in enumerate(cols)}
required=('id','type','x','y','z','vx','vy','vz')
assert all(k in idx for k in required)
rows=np.asarray([[float(v) for v in line.split()] for line in lines[header+1:header+1+n]],dtype=float)
assert len(rows)==n
is_ag=rows[:,idx['type']].astype(int)==4
assert int(is_ag.sum())==N_AG_EXPECTED
x,y,z=(rows[:,idx[k]] for k in ('x','y','z'))
dx=(x-CENTER_X+BOX_XY/2)%BOX_XY-BOX_XY/2
dy=(y-CENTER_Y+BOX_XY/2)%BOX_XY-BOX_XY/2
wxy=np.exp(-0.5*(dx*dx+dy*dy)/SIGMA_XY**2)
wz=np.exp(-np.maximum(SURFACE_Z-z,0.0)/DEPTH)
weights=np.clip(wxy*wz,0.0,1.0)[is_ag]
velcols=[idx[k] for k in ('vx','vy','vz')]
oldv=rows[is_ag][:,velcols].copy()
rng=np.random.default_rng(SEED)
nominal_dt=TARGET_E/(1.5*KB*weights.sum())
sigma_v=np.sqrt(KB*nominal_dt/(MASS_AG*KINETIC_FACTOR))
kick=rng.normal(size=(N_AG_EXPECTED,3))*sigma_v*np.sqrt(weights[:,None])
kick-=kick.mean(axis=0,keepdims=True)
a=0.5*MASS_AG*KINETIC_FACTOR*np.sum(kick*kick)
b=MASS_AG*KINETIC_FACTOR*np.sum(oldv*kick)
scale=(-b+np.sqrt(b*b+4*a*TARGET_E))/(2*a)
kick*=scale
newv=oldv+kick
old_ke=0.5*MASS_AG*KINETIC_FACTOR*np.sum(oldv*oldv)
new_ke=0.5*MASS_AG*KINETIC_FACTOR*np.sum(newv*newv)
assert abs((new_ke-old_ke)-TARGET_E)<1e-7
rows[np.ix_(np.flatnonzero(is_ag),velcols)]=newv
out=ROOT/'seed41_pulse2_1479eV_start.dump'
with out.open('w') as f:
    f.write(f'ITEM: TIMESTEP\n{STEP}\nITEM: NUMBER OF ATOMS\n{n}\nITEM: BOX BOUNDS pp pp ff\n0 200\n0 200\n0 500\nITEM: ATOMS id type x y z vx vy vz\n')
    for row in rows:
        f.write(f'{int(row[idx["id"]])} {int(row[idx["type"]])} '+' '.join(f'{row[j]:.15g}' for j in (idx['x'],idx['y'],idx['z'],*velcols))+'\n')
info={'source_dump':str(SRC.relative_to(ROOT.parent)),'source_time_ps':TIME_PS,'timestep':STEP,
      'total_atoms':N_EXPECTED,'Ag_atoms':N_AG_EXPECTED,'seed':SEED,'center_A':[CENTER_X,CENTER_Y,SURFACE_Z],
      'sigma_xy_A':SIGMA_XY,'attenuation_depth_A':DEPTH,
      'target_added_kinetic_energy_eV':TARGET_E,'actual_added_kinetic_energy_eV':float(new_ke-old_ke),
      'spatial_weight_sum':float(weights.sum()),'effective_weighted_atoms':float(weights.sum()**2/(weights@weights)),
      'nominal_peak_delta_temperature_K':float(nominal_dt),'kick_scale':float(scale),
      'pulse_kind':'second localized zero-net-momentum random velocity impulse; no TTM'}
(ROOT/'seed41_pulse2_parameters.json').write_text(json.dumps(info,indent=2)+'\n')
print(json.dumps(info,indent=2))
