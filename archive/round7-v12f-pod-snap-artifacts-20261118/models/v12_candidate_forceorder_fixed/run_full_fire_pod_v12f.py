#!/workspace/venvs/dft-labels/bin/python
"""Isolated 100-iteration FIRE pre-relaxation of the full 347k model."""
import os,subprocess,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
DATA=Path('/workspace/potential-fit/round7_pod_snap/models/target_200x200x300_qsgs/md_runs/attempt_008_POD_TTM_mpi_big/start_wrapped.data')
LMP=Path('/workspace/software/lammps-pod-snap-clean/build/lmp')
OUT=HERE/'full_fire_347k'/'POD'; OUT.mkdir(parents=True,exist_ok=True)
inp=OUT/'in.fire_100'
inp.write_text(f'''clear
units metal
atom_style atomic
boundary p p m
processors 5 1 1
atom_modify sort 0 0
read_data {DATA}
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
group fixed region fixed_layer
fix freeze fixed setforce 0.0 0.0 0.0
variable fmag atom sqrt(fx*fx+fy*fy+fz*fz)
compute fmax all reduce max v_fmag
thermo 10
thermo_style custom step atoms pe c_fmax
thermo_modify lost error flush yes format float %.16g
min_style fire
min_modify dmax 0.02
minimize 0.0 0.5 100 5000
write_dump all custom {OUT/'forces_minimized.dump'} id type x y z fx fy fz modify sort id
write_data {OUT/'minimized.data'}
''')
env=os.environ.copy(); env['LD_LIBRARY_PATH']='/workspace/software/openmpi-dev-5.0.7/usr/lib/x86_64-linux-gnu/openmpi/lib:'+str(LMP.parent)+':'+env.get('LD_LIBRARY_PATH',''); env['OMP_NUM_THREADS']='1'
cmd=['/usr/bin/mpirun','-np','5',str(LMP),'-in',str(inp),'-log','log.lammps']
t=time.time(); p=subprocess.run(cmd,cwd=OUT,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
(OUT/'stdout.log').write_text(p.stdout); (OUT/'run_metadata.txt').write_text(f'POD v12f; full model N=347464; FIRE 100 max iterations, ftol=0.5 eV/A, dmax=0.02 A; fixed z=0-20 A; MPI ranks=5; wall={time.time()-t:.3f}s; rc={p.returncode}\n')
print(p.stdout[-5000:])
if p.returncode: raise SystemExit(p.returncode)
