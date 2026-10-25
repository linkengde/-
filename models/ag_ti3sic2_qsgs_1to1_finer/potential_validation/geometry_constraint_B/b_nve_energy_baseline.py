import json, os, sys, time
from pathlib import Path
sys.path.extend(['/workspace/.venvs/lammps/lib/python3.12/site-packages'])
from mpi4py import MPI
from lammps import lammps
root=Path('/workspace/-/models/ag_ti3sic2_qsgs_1to1_finer'); b=root/'potential_validation/geometry_constraint_B'; cand=root/'potential_validation/relaxed_backbone_candidate'; os.chdir(root)
lmp=lammps(cmdargs=['-screen','none','-log','none','-sf','omp'],comm=MPI.COMM_WORLD)
commands=['units metal','atom_style atomic','boundary p p f','atom_modify sort 0 0',f'read_data {cand}/ag_ti3sic2_qsgs.data.gz',f'read_dump {b}/ag_minimized_100.dump.gz 100 x y z box no replace yes',f'include {b}/potential.geometry_wca_fixed_skeleton.inc','neighbor 2.0 bin','neigh_modify delay 0 every 1 check yes','region r_bottom block INF INF INF INF 0.0 20.0 units box','group bottom region r_bottom','group Ag type 4','group skeleton type 1 2 3','group mobile_Ag subtract Ag bottom','fix freeze_bottom bottom setforce 0 0 0','fix freeze_skeleton skeleton setforce 0 0 0','mass 4 107.8682','velocity mobile_Ag create 300.0 190427 mom yes rot no dist gaussian','timestep 0.0005','fix integrate mobile_Ag nve','thermo 1','thermo_style custom step pe ke etotal','thermo_modify flush yes','run 0']
for c in commands:lmp.command(c)
e0=float(lmp.get_thermo('etotal')); print('initial_etotal_eV',e0,flush=True)
nsteps=int(os.environ.get('AG_NVE_BASELINE_STEPS', '100'))
if nsteps <= 0:
 raise ValueError('AG_NVE_BASELINE_STEPS must be positive')
for i in range(1,nsteps+1):
 lmp.command('run 1 pre no post no')
 if i%20==0 or i==nsteps:
  pe=float(lmp.get_thermo('pe')); ke=float(lmp.get_thermo('ke')); et=float(lmp.get_thermo('etotal'))
  print(f'step={i} time_ps={i*0.0005:.5f} pe={pe:.6f} ke={ke:.6f} etotal={et:.6f} drift_eV={et-e0:.6f}',flush=True)
lmp.command('run 0 pre no post no'); ef=float(lmp.get_thermo('etotal')); print('final_after_run0_eV',ef,'final_drift_eV',ef-e0,flush=True)
lmp.close(); MPI.Finalize()
