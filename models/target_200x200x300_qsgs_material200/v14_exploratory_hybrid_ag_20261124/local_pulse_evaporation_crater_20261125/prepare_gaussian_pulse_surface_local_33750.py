#!/usr/bin/env python3
"""Add a one-time, zero-net-momentum Gaussian thermal kick to Ag velocities."""
from pathlib import Path
import json
import numpy as np

BASE = Path(__file__).resolve().parent
src = BASE.parent / 'simple_ag_fixed_skeleton_20261125/reheat_melt_assessment_20261125' / 'trajectories' / 'ag_reheat_q6_ptm_msd_final.dump'
out = BASE / 'prepared' / 'local_gaussian_pulse_surface_local_33750_start.dump'
center_x = center_y = 100.0
surface_z = 203.0
sigma_xy = 12.0
attenuation_depth = 6.0
peak_delta_T = 33750.0
seed = 20261129
mass_ag = 107.8682  # g/mol
kB = 8.617333262145e-5  # eV/K
kinetic_factor = 1.03642696526e-4  # eV per [g/mol * (A/ps)^2], multiplying 0.5*m*v^2
box_xy = 200.0

lines = src.read_text().splitlines()
assert lines[0] == 'ITEM: TIMESTEP'
step = int(lines[1])
assert lines[2] == 'ITEM: NUMBER OF ATOMS'
n = int(lines[3])
assert lines[4].startswith('ITEM: BOX BOUNDS')
header_i = 8
assert lines[header_i].startswith('ITEM: ATOMS ')
cols = lines[header_i].split()[2:]
idx = {name:i for i,name in enumerate(cols)}
for name in ['id','type','x','y','z','vx','vy','vz']:
    assert name in idx, name
rows = np.array([[float(v) for v in line.split()] for line in lines[header_i+1:header_i+1+n]], dtype=float)
assert len(rows) == n == 458248
is_ag = rows[:,idx['type']].astype(int) == 4
assert int(is_ag.sum()) == 106396
x,y,z = (rows[:,idx[k]] for k in ('x','y','z'))
dx = (x-center_x+box_xy/2) % box_xy - box_xy/2
dy = (y-center_y+box_xy/2) % box_xy - box_xy/2
r2 = dx*dx + dy*dy
wxy = np.exp(-0.5*r2/sigma_xy**2)
wz = np.exp(-np.maximum(surface_z-z,0.0)/attenuation_depth)
w = np.clip(wxy*wz,0.0,1.0)
weights = w[is_ag]
rng = np.random.default_rng(seed)
sigma_v = np.sqrt(kB*peak_delta_T/(mass_ag*kinetic_factor))
kicks = rng.normal(size=(int(is_ag.sum()),3))*sigma_v*np.sqrt(weights[:,None])
kicks -= kicks.mean(axis=0,keepdims=True)
velcols=[idx['vx'],idx['vy'],idx['vz']]
oldv = rows[is_ag][:,velcols].copy()
newv = oldv + kicks
rows[np.ix_(np.flatnonzero(is_ag),velcols)] = newv
old_ke = 0.5*mass_ag*kinetic_factor*np.sum(oldv*oldv)
new_ke = 0.5*mass_ag*kinetic_factor*np.sum(newv*newv)
r2_ag=r2[is_ag]; z_ag=z[is_ag]
spot_mask=(r2_ag < sigma_xy**2)&(z_ag > surface_z-attenuation_depth)
spot_old=0.5*mass_ag*kinetic_factor*np.sum(oldv[spot_mask]**2)
spot_new=0.5*mass_ag*kinetic_factor*np.sum(newv[spot_mask]**2)
info = {
  'source_snapshot_timestep': step,
  'source_snapshot_time_ps': 1.7,
  'Ag_atoms': int(is_ag.sum()),
  'source_profile': 'exp(-r_xy^2/(2 sigma_xy^2)) * exp(-max(z_surface-z,0)/attenuation_depth)',
  'center_A': [center_x,center_y,surface_z],
  'sigma_xy_A': sigma_xy,
  'attenuation_depth_A': attenuation_depth,
  'peak_delta_temperature_K': peak_delta_T,
  'seed': seed,
  'velocity_noise_sigma_A_per_ps_at_peak': float(sigma_v),
  'sum_gaussian_weights_Ag': float(weights.sum()),
  'pulse_added_kinetic_energy_eV': float(new_ke-old_ke),
  'source_zone_Ag_count': int(spot_mask.sum()),
  'source_zone_kinetic_energy_before_eV': float(spot_old),
  'source_zone_kinetic_energy_after_eV': float(spot_new),
  'initial_full_system_box_z_hi_A': 300.0,
  'prepared_box_z_hi_A': 500.0,
}
out.parent.mkdir(parents=True,exist_ok=True)
with out.open('w') as f:
    f.write('ITEM: TIMESTEP\n')
    f.write(f'{step}\nITEM: NUMBER OF ATOMS\n{n}\nITEM: BOX BOUNDS pp pp ff\n')
    f.write('0 200\n0 200\n0 500\n')
    f.write('ITEM: ATOMS id type x y z vx vy vz\n')
    keep=[idx[k] for k in ('id','type','x','y','z','vx','vy','vz')]
    for row in rows:
        f.write(f'{int(row[idx["id"]])} {int(row[idx["type"]])} ' + ' '.join(f'{row[j]:.12g}' for j in keep[2:]) + '\n')
(BASE/'prepared'/'pulse_surface_local_33750_parameters.json').write_text(json.dumps(info,indent=2)+'\n')
print(json.dumps(info,indent=2))
