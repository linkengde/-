#!/usr/bin/env python3
"""Short, same-geometry 64-atom cold NVE comparison; exploratory only."""
import json, os, re, subprocess
from pathlib import Path
import numpy as np
from ase.io import read, write

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'v14_stage2_p0_ag_dft_diagnostics'/'interface_cold_nve'
OUT.mkdir(parents=True,exist_ok=True)
SRC=ROOT/'fit_v13'/'validation'/'V13_VAL_INTERFACE_01.xyz'
FIT=ROOT/'fit_v14_targeted_ag_force5'
POT=Path(os.environ.get('V14_POT_DIR', '/workspace/-/models/ag_ti3sic2_qsgs_1to1_finer/potentials'))
LMP=os.environ.get('LAMMPS_BIN', '/workspace/venvs/podfit/bin/lmp')
data=OUT/'heldout_interface.data'
atoms=read(SRC,format='extxyz')
write(data,atoms,format='lammps-data',atom_style='atomic',units='metal',specorder=['Ag','C','Si','Ti'])
models={
 'V14_force5':f'''pair_style pod
pair_coeff * * {FIT/'AgTiSiC_POD_v14.pod'} {FIT/'v14_targeted_ag_f5_coefficients.pod'} Ag C Si Ti''',
 'hybrid_eam_TSCBOP_fitted_Morse':f'''pair_style hybrid eam tersoff morse 8.0
pair_coeff 1 1 eam {POT/'Ag_u3.eam'}
pair_coeff * * tersoff {POT/'TSCBOP.tersoff'} NULL C Si Ti
pair_coeff 1 2 morse 0.020000 7.804451 1.400000 8.0
pair_coeff 1 3 morse 0.020000 3.107766 1.400605 8.0
pair_coeff 1 4 morse 0.678644 1.155035 2.989514 8.0''',
 'hybrid_Morse_ZBL_short_range_guard':f'''pair_style hybrid/overlay eam tersoff morse 8.0 zbl 1.5 2.4
pair_coeff 1 1 eam {POT/'Ag_u3.eam'}
pair_coeff * * tersoff {POT/'TSCBOP.tersoff'} NULL C Si Ti
pair_coeff 1 2 morse 0.020000 7.804451 1.400000 8.0
pair_coeff 1 3 morse 0.020000 3.107766 1.400605 8.0
pair_coeff 1 4 morse 0.678644 1.155035 2.989514 8.0
pair_coeff 1 2 zbl 47 6
pair_coeff 1 3 zbl 47 14
pair_coeff 1 4 zbl 47 22'''
}
rows=[]
for name,pair in models.items():
    d=OUT/name; d.mkdir(exist_ok=True)
    inp=d/'in.nve'; log=d/'log.lammps'; dump=d/'traj.dump'
    inp.write_text(f'''clear
units metal
atom_style atomic
boundary p p p
atom_modify sort 0 0
read_data {data}
mass 1 107.8682
mass 2 12.011
mass 3 28.0855
mass 4 47.867
{pair}
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
thermo 10
thermo_style custom step time temp pe ke etotal press
velocity all create 300.0 24601 mom yes rot yes dist gaussian
fix int all nve
timestep 0.001
dump tr all custom 10 {dump} id type x y z vx vy vz fx fy fz
dump_modify tr sort id
run 100
write_data {d/'final.data'}
''')
    env=os.environ.copy();env.update({'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
    p=subprocess.run([LMP,'-in',str(inp),'-log',str(log),'-screen','none'],cwd=d,env=env,capture_output=True,text=True)
    failure=''
    if p.returncode: failure=(p.stdout+'\n'+p.stderr)[-3000:]
    text=log.read_text(); therm=[]
    for block in re.split(r'\n\s*Step\s+Time\s+Temp\s+PotEng\s+KinEng\s+TotEng\s+Press\s*\n',text)[1:]:
        for line in block.splitlines():
            vals=line.split()
            if len(vals)>=7 and vals[0].isdigit():
                try: therm.append([float(x) for x in vals[:7]])
                except ValueError: pass
    # Read the simple orthogonal/triclinic dump via ASE; preserve element map.
    frames=read(dump,index=':',format='lammps-dump-text',specorder=['Ag','C','Si','Ti'])
    fmax=max(float(np.linalg.norm(fr.get_forces(),axis=1).max()) for fr in frames)
    min_ag_ag=[]; min_ag_sk=[]; nframes=[]
    for fr in frames:
        sy=np.asarray(fr.get_chemical_symbols());dist=fr.get_all_distances(mic=True)
        ai=np.where(sy=='Ag')[0]; sk=np.where(sy!='Ag')[0]
        min_ag_ag.append(float(dist[np.ix_(ai,ai)][np.triu_indices(len(ai),1)].min()))
        min_ag_sk.append(float(dist[np.ix_(ai,sk)].min()));nframes.append(len(fr))
    arr=np.asarray(therm)
    # thermo rows are [step,time,temp,pe,ke,etotal,press].
    result={'model':name,'status':'completed' if p.returncode==0 else 'failed_during_short_run',
            'failure_tail':failure,'initial_atoms':len(atoms),'final_atoms':int(nframes[-1]),
            'duration_ps':0.1,'dt_fs':1.0,'dump_frames':len(frames),
            'T_initial_K':float(arr[0,2]) if len(arr) else None,'T_final_K':float(arr[-1,2]) if len(arr) else None,
            'T_max_K':float(arr[:,2].max()) if len(arr) else None,
            'total_energy_drift_eV_per_atom':float((arr[-1,5]-arr[0,5])/len(atoms)) if len(arr) else None,
            'max_atomic_force_eV_A_over_dump':fmax,'min_Ag_Ag_A':min(min_ag_ag),
            'min_Ag_skeleton_A':min(min_ag_sk),'note':'100-step 300 K NVE screen on one held-out 64-atom DFT interface geometry; not a phase or long-MD validation'}
    (d/'summary.json').write_text(json.dumps(result,indent=2)+'\n'); rows.append(result)
(OUT/'summary.json').write_text(json.dumps({'source':str(SRC),'results':rows},indent=2)+'\n')
print(json.dumps(rows,indent=2))
