#!/usr/bin/env python3
"""FIRE screen for the full-height 458k geometry using POD v12g."""
from __future__ import annotations

import os
import subprocess
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / 'Ag_Ti3SiC2_bicontinuous_200x200x300.data'
OUT = HERE / 'fire_v12g_458k'
OUT.mkdir(exist_ok=True)
LMP = Path('/workspace/software/lammps-pod-snap-clean/build/lmp')
PODROOT = Path('/workspace/potential-fit/round7_pod_snap/models/v12_candidate_forceorder_fixed/v12g_si04_candidate/POD')
inp = OUT / 'in.fire_100'
inp.write_text(f'''clear
units metal
atom_style atomic
boundary p p m
processors 4 1 1
atom_modify sort 0 0
read_data {DATA}
# Source types Ti/Si/C/Ag -> model types Ag/C/Si/Ti.
group src_Ti type 1
group src_Si type 2
group src_C type 3
group src_Ag type 4
set group src_Ti type 4
set group src_Si type 3
set group src_C type 2
set group src_Ag type 1
mass 1 107.8682
mass 2 12.011
mass 3 28.085
mass 4 47.867
pair_style pod
pair_coeff * * {PODROOT/'AgTiSiC_POD_v12g.pod'} {PODROOT/'pod_v12g_coefficients.pod'} Ag C Si Ti
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
group skeleton type 2 3 4
region fixed_layer block INF INF INF INF 0.0 20.0 units box
group fixed_z region fixed_layer
group fixed intersect fixed_z skeleton
fix freeze fixed setforce 0.0 0.0 0.0
variable fmag atom sqrt(fx*fx+fy*fy+fz*fz)
compute fmax all reduce max v_fmag
thermo 10
thermo_style custom step atoms pe c_fmax
thermo_modify lost error flush yes format float %.16g
dump xyz all custom 10 {OUT/'xyz_by_10steps.dump'} id x y z
dump_modify xyz sort id
min_style fire
min_modify dmax 0.02
minimize 0.0 0.5 100 5000
undump xyz
write_dump all custom {OUT/'forces_minimized.dump'} id type x y z fx fy fz modify sort id
write_data {OUT/'minimized.data'}
''')
env = os.environ.copy()
env['LD_LIBRARY_PATH'] = '/workspace/software/openmpi-dev-5.0.7/usr/lib/x86_64-linux-gnu/openmpi/lib:' + str(LMP.parent) + ':' + env.get('LD_LIBRARY_PATH', '')
env['OMP_NUM_THREADS'] = '1'
cmd = ['/usr/bin/mpirun', '-np', '4', str(LMP), '-in', str(inp), '-log', 'log.lammps']
def lower_priority() -> None:
    os.nice(10)
start = time.time()
result = subprocess.run(cmd, cwd=OUT, env=env, stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT, text=True, preexec_fn=lower_priority)
(OUT/'stdout.log').write_text(result.stdout)
(OUT/'run_metadata.txt').write_text(
    f'POD v12g; full-height model N=458248; FIRE max=100 iterations, ftol=0.5 eV/A, '
    f'dmax=0.02 A; fixed bottom Ti3SiC2 skeleton z=0-20 A; MPI ranks=4; '
    f'wall={time.time()-start:.3f}s; rc={result.returncode}\n')
print(result.stdout[-5000:])
raise SystemExit(result.returncode)
