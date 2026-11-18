#!/workspace/venvs/dft-labels/bin/python
"""Continue the isolated FIRE relaxation from the prior 100-step checkpoint."""
import os,subprocess,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
SRC=HERE/'small_fire_14007/POD/minimized.data'
LMP=Path('/workspace/software/lammps-pod-snap-clean/build/lmp')
OUT=HERE/'small_fire_14007_iter2'/'POD'; OUT.mkdir(parents=True,exist_ok=True)
inp=OUT/'in.fire_100'
inp.write_text(f'''clear
units metal
atom_style atomic
boundary f f f
atom_modify sort 0 0
read_data {SRC}
pair_style pod
pair_coeff * * {HERE/'POD/AgTiSiC_POD_v12f.pod'} {HERE/'POD/pod_v12f_coefficients.pod'} Ag C Si Ti
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
fix rim side setforce 0.0 0.0 0.0
variable fmag atom sqrt(fx*fx+fy*fy+fz*fz)
compute fmax all reduce max v_fmag
thermo 10
thermo_style custom step atoms pe c_fmax
thermo_modify lost error flush yes format float %.16g
min_style fire
min_modify dmax 0.02
minimize 0.0 0.2 100 2000
write_dump all custom {OUT/'forces_minimized.dump'} id type x y z fx fy fz modify sort id
write_data {OUT/'minimized.data'}
''')
env=os.environ.copy(); env['LD_LIBRARY_PATH']='/workspace/software/openmpi-dev-5.0.7/usr/lib/x86_64-linux-gnu/openmpi/lib:'+str(LMP.parent)+':'+env.get('LD_LIBRARY_PATH',''); env['OMP_NUM_THREADS']='1'
t=time.time(); p=subprocess.run([str(LMP),'-in',str(inp),'-log','log.lammps'],cwd=OUT,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
(OUT/'stdout.log').write_text(p.stdout); (OUT/'run_metadata.txt').write_text(f'Continue FIRE from 100-step v12f checkpoint; 100 further iterations, ftol=0.2 eV/A, dmax=0.02 A; wall={time.time()-t:.3f}s; rc={p.returncode}\n')
print(p.stdout[-4000:])
if p.returncode: raise SystemExit(p.returncode)
