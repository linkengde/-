#!/usr/bin/env python3
"""Prepare equal-dose 1.109 keV surface-local Ag pulses for two seeds."""
from pathlib import Path
import json
import numpy as np

BASE=Path(__file__).resolve().parent
SRC=BASE.parent/'simple_ag_fixed_skeleton_20261125/reheat_melt_assessment_20261125/trajectories/ag_reheat_q6_ptm_msd_final.dump'
OUT=BASE/'prepared'
SEEDS=(20261131,20261132)
TARGET_E=739.6246479414622*1.5
CENTER_X=CENTER_Y=100.0
SURFACE_Z=203.0
SIGMA_XY=12.0
DEPTH=6.0
PEAK_DELTA_T=33750.0*1.5
MASS_AG=107.8682
KB=8.617333262145e-5
KINETIC_FACTOR=1.03642696526e-4
BOX_XY=200.0

lines=SRC.read_text().splitlines()
assert lines[0]=='ITEM: TIMESTEP' and int(lines[1])==1700
n=int(lines[3]); assert n==458248 and lines[8].startswith('ITEM: ATOMS ')
cols=lines[8].split()[2:]; idx={k:i for i,k in enumerate(cols)}
required=('id','type','x','y','z','vx','vy','vz')
assert all(k in idx for k in required)
base_rows=np.asarray([[float(v) for v in s.split()] for s in lines[9:9+n]],dtype=float)
assert len(base_rows)==n
is_ag=base_rows[:,idx['type']].astype(int)==4
assert int(is_ag.sum())==106396
x,y,z=(base_rows[:,idx[k]] for k in ('x','y','z'))
dx=(x-CENTER_X+BOX_XY/2)%BOX_XY-BOX_XY/2
dy=(y-CENTER_Y+BOX_XY/2)%BOX_XY-BOX_XY/2
wxy=np.exp(-.5*(dx*dx+dy*dy)/SIGMA_XY**2)
wz=np.exp(-np.maximum(SURFACE_Z-z,0.0)/DEPTH)
weights=np.clip(wxy*wz,0.0,1.0)[is_ag]
velcols=[idx[k] for k in ('vx','vy','vz')]
OUT.mkdir(parents=True,exist_ok=True)

for seed in SEEDS:
    rows=base_rows.copy()
    rng=np.random.default_rng(seed)
    sigma_v=np.sqrt(KB*PEAK_DELTA_T/(MASS_AG*KINETIC_FACTOR))
    kick=rng.normal(size=(int(is_ag.sum()),3))*sigma_v*np.sqrt(weights[:,None])
    kick-=kick.mean(axis=0,keepdims=True)
    oldv=rows[is_ag][:,velcols].copy()
    a=.5*MASS_AG*KINETIC_FACTOR*np.sum(kick*kick)
    b=MASS_AG*KINETIC_FACTOR*np.sum(oldv*kick)
    scale=(-b+np.sqrt(b*b+4*a*TARGET_E))/(2*a)
    kick*=scale
    newv=oldv+kick
    rows[np.ix_(np.flatnonzero(is_ag),velcols)]=newv
    old_ke=.5*MASS_AG*KINETIC_FACTOR*np.sum(oldv*oldv)
    new_ke=.5*MASS_AG*KINETIC_FACTOR*np.sum(newv*newv)
    assert abs((new_ke-old_ke)-TARGET_E)<1e-8
    out=OUT/f'local_gaussian_pulse_1109eV_seed{seed}_start.dump'
    with out.open('w') as f:
        f.write(f'ITEM: TIMESTEP\n1700\nITEM: NUMBER OF ATOMS\n{n}\nITEM: BOX BOUNDS pp pp ff\n0 200\n0 200\n0 500\nITEM: ATOMS id type x y z vx vy vz\n')
        for row in rows:
            f.write(f'{int(row[idx["id"]])} {int(row[idx["type"]])} '+' '.join(f'{row[j]:.15g}' for j in (idx['x'],idx['y'],idx['z'],*velcols))+'\n')
    info={'source_snapshot_timestep':1700,'source_snapshot_time_ps':1.7,'total_atoms':n,'Ag_atoms':int(is_ag.sum()),
          'seed':seed,'center_A':[CENTER_X,CENTER_Y,SURFACE_Z],'sigma_xy_A':SIGMA_XY,'attenuation_depth_A':DEPTH,
          'nominal_peak_delta_temperature_K':PEAK_DELTA_T,'target_added_kinetic_energy_eV':TARGET_E,
          'actual_added_kinetic_energy_eV':float(new_ke-old_ke),'kick_scale':float(scale),
          'profile':'exp(-r_xy^2/(2 sigma_xy^2))*exp(-max(surface_z-z,0)/depth)',
          'pulse_kind':'single instantaneous zero-net-momentum random velocity increment; no TTM',
          'prepared_box_z_hi_A':500.0}
    (OUT/f'pulse_1109eV_seed{seed}_parameters.json').write_text(json.dumps(info,indent=2)+'\n')
    print(json.dumps(info))
