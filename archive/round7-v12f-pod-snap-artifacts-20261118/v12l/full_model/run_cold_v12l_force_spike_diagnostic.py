#!/usr/bin/env python3
"""Replay POD v12l cold trajectory from its 10 fs checkpoint and halt at Fmax 500 eV/A."""
import os, subprocess, time
from pathlib import Path
HERE=Path(__file__).resolve().parent
LMP=Path("/workspace/software/lammps-pod-snap-clean/build/lmp")
POD=Path("/workspace/potential-fit/round7_pod_snap/models/v12_candidate_forceorder_fixed/v12l_ticpair_candidate/POD")
DATA=HERE/"cold_v12l_initial_458k_1fs/end_0p01ps.data"
OUT=HERE/"cold_v12l_force_spike_diagnostic"
OUT.mkdir()
inp=OUT/"in.force_spike_diag"
dump=OUT/"first_force_over_500.dump"
inp.write_text(f'''clear
units metal
atom_style atomic
boundary p p m
processors 4 1 1
atom_modify sort 0 0
read_data {DATA}
mass 1 107.8682
mass 2 12.011
mass 3 28.085
mass 4 47.867
pair_style pod
pair_coeff * * {POD/'AgTiSiC_POD_v12l.pod'} {POD/'pod_v12l_coefficients.pod'} Ag C Si Ti
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
group skeleton type 2 3 4
region fixed_layer block INF INF INF INF 0.0 20.0 units box
region bath_layer block INF INF INF INF 20.0 40.0 units box
group fixed_z region fixed_layer
group fixed intersect skeleton fixed_z
group bath region bath_layer
group mobile subtract all fixed
velocity fixed set 0.0 0.0 0.0
fix freeze fixed setforce 0.0 0.0 0.0
compute Tmobile mobile temp
variable fmag atom sqrt(fx*fx+fy*fy+fz*fz)
compute fmax all reduce max v_fmag
variable fmax_now equal c_fmax
thermo 1
thermo_style custom step time atoms c_Tmobile pe ke etotal c_fmax
thermo_modify temp Tmobile lost error flush yes format float %.16g
timestep 0.001
fix integrate mobile nve
fix heat_bath bath langevin 300.0 300.0 0.1 639829 tally yes
fix stop all halt 1 v_fmax_now > 500.0 error soft message yes
run 40
write_dump all custom {dump} id type x y z vx vy vz fx fy fz modify sort id
write_data {OUT/'first_force_over_500.data'}
''')
env=os.environ.copy(); env["LD_LIBRARY_PATH"]="/workspace/software/openmpi-dev-5.0.7/usr/lib/x86_64-linux-gnu/openmpi/lib:"+str(LMP.parent)+":"+env.get("LD_LIBRARY_PATH",""); env["OMP_NUM_THREADS"]="1"
start=time.time()
res=subprocess.run(["/usr/bin/mpirun","-np","4",str(LMP),"-in",str(inp),"-log","log.lammps"],cwd=OUT,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,preexec_fn=lambda:os.nice(10))
(OUT/"stdout.log").write_text(res.stdout)
(OUT/"run_metadata.txt").write_text(f"POD v12l; exact 10-fs checkpoint replay; Langevin seed=639829; halts on first Fmax >500 eV/A; wall_s={time.time()-start:.3f}; rc={res.returncode}\n")
print(res.stdout[-5000:])
raise SystemExit(res.returncode)
