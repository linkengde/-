#!/workspace/venvs/dft-labels/bin/python
"""Ten-step 1 fs no-heat NVE integration screen on the old Round7 crop."""
import os, subprocess, time
from pathlib import Path
HERE=Path(__file__).resolve().parent
SOURCE=Path('/workspace/potential-fit/round7_pod_snap/models/final_v9/ideal_sink_heating_ag4k/attempt_004/start.data')
LMP=Path('/workspace/software/lammps-pod-snap-clean/build/lmp')
OUT=HERE/'small_cold_v12h_14007_1fs'/'SNAP'
OUT.mkdir(parents=True,exist_ok=True)
inp=OUT/'in.nve_0p01ps_1fs'
inp.write_text(f'''clear
units metal
atom_style atomic
boundary f f f
atom_modify sort 0 0
read_data {SOURCE}
pair_style snap
pair_coeff * * {HERE/'SNAP/AgTiSiC_SNAP_v12h.snapcoeff'} {HERE/'SNAP/AgTiSiC_SNAP_v12h.snapparam'} Ag C Si Ti
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
region xlo block 0 7 INF INF INF INF units box
region xhi block 43 50 INF INF INF INF units box
region ylo block INF INF 0 7 INF INF units box
region yhi block INF INF 43 50 INF INF units box
group gxlo region xlo
group gxhi region xhi
group gylo region ylo
group gyhi region yhi
group side union gxlo gxhi gylo gyhi
velocity side set 0.0 0.0 0.0
fix rim side setforce 0.0 0.0 0.0
group mobile subtract all side
compute Tmobile mobile temp
variable fmag atom sqrt(fx*fx+fy*fy+fz*fz)
compute fmax all reduce max v_fmag
timestep 0.001
thermo 1
thermo_style custom step time atoms c_Tmobile pe ke etotal c_fmax
thermo_modify temp Tmobile lost error flush yes format float %.16g
dump traj all custom 1 {OUT/'trajectory.dump'} id type x y z vx vy vz fx fy fz
dump_modify traj sort id format float %.12g
run 0
velocity mobile create 300.0 738291 mom yes rot no dist gaussian
fix integ mobile nve
run 10
unfix integ
undump traj
write_data {OUT/'end_0p01ps.data'}
''')
env=os.environ.copy(); env['LD_LIBRARY_PATH']='/workspace/software/openmpi-dev-5.0.7/usr/lib/x86_64-linux-gnu/openmpi/lib:'+str(LMP.parent)+':'+env.get('LD_LIBRARY_PATH',''); env['OMP_NUM_THREADS']='1'
t=time.time(); p=subprocess.run([str(LMP),'-in',str(inp),'-log','log.lammps'],cwd=OUT,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
(OUT/'stdout.log').write_text(p.stdout)
(OUT/'run_metadata.txt').write_text(f'Potential: SNAP v12h\nSource: {SOURCE}\nAtoms: 14007\nMethod: run 0 then NVE, no heat source or thermostat; mobile atoms initialized at 300 K\nTimestep: 0.001 ps (1 fs)\nNVE duration: 0.01 ps (10 steps)\nWall time: {time.time()-t:.3f} s\nReturn code: {p.returncode}\n')
print(p.stdout[-5000:])
if p.returncode: raise SystemExit(p.returncode)
