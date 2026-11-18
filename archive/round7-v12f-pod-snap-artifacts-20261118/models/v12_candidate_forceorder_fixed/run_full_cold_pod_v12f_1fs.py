#!/workspace/venvs/dft-labels/bin/python
"""0.01 ps no-laser cold integration screen for the full 347,464-atom model."""
import os,subprocess,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
DATA=Path('/workspace/potential-fit/round7_pod_snap/models/target_200x200x300_qsgs/md_runs/attempt_008_POD_TTM_mpi_big/start_wrapped.data')
LMP=Path('/workspace/software/lammps-pod-snap-clean/build/lmp')
OUT=HERE/'full_cold_347k_1fs'/'POD'; OUT.mkdir(parents=True,exist_ok=True)
inp=OUT/'in.nve_cold_0p01ps_1fs'
inp.write_text(f'''clear
units metal
atom_style atomic
boundary p p m
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
pair_coeff * * {HERE/'POD/AgTiSiC_POD_v12f.pod'} {HERE/'POD/pod_v12f_coefficients.pod'} Ag C Si Ti
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
region fixed_layer block INF INF INF INF 0.0 20.0 units box
region bath_layer block INF INF INF INF 20.0 40.0 units box
group fixed region fixed_layer
group bath region bath_layer
group mobile subtract all fixed
velocity mobile create 300.0 882733 mom yes rot no dist gaussian
velocity fixed set 0.0 0.0 0.0
fix freeze fixed setforce 0.0 0.0 0.0
compute Tmobile mobile temp
variable fmag atom sqrt(fx*fx+fy*fy+fz*fz)
compute fmax all reduce max v_fmag
timestep 0.001
thermo 1
thermo_style custom step time atoms c_Tmobile pe ke etotal c_fmax
thermo_modify temp Tmobile lost error flush yes format float %.16g
fix integrate mobile nve
fix heat_bath bath langevin 300.0 300.0 0.1 639827 tally yes
run 0
dump snap all custom 10 {OUT/'snapshots.dump'} id type x y z fx fy fz
dump_modify snap sort id
run 10
unfix heat_bath
unfix integrate
undump snap
write_data {OUT/'end_0p01ps.data'}
''')
env=os.environ.copy(); env['LD_LIBRARY_PATH']='/workspace/software/openmpi-dev-5.0.7/usr/lib/x86_64-linux-gnu/openmpi/lib:'+str(LMP.parent)+':'+env.get('LD_LIBRARY_PATH',''); env['OMP_NUM_THREADS']='1'
cmd=['/usr/bin/mpirun','-np','4',str(LMP),'-in',str(inp),'-log','log.lammps']
t=time.time(); p=subprocess.run(cmd,cwd=OUT,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
(OUT/'stdout.log').write_text(p.stdout); (OUT/'run_metadata.txt').write_text(f'POD v12f; full model N=347464; dt=0.001 ps (1 fs); duration=0.01 ps (10 steps); no laser; z=0-20 fixed, z=20-40 Langevin 300K, mobile remainder; ranks=4; wall={time.time()-t:.3f}s; rc={p.returncode}\n')
print(p.stdout[-6000:])
if p.returncode: raise SystemExit(p.returncode)
