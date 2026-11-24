#!/usr/bin/env python3
"""Minimize then short cold NVE on the 64-atom held-out ideal interface."""
import json, os, re, subprocess
from pathlib import Path
import numpy as np
from ase.io import read, write

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'v14_stage2_p0_ag_dft_diagnostics'/'interface_hybrid_fire_nve'
OUT.mkdir(parents=True,exist_ok=True)
SRC=ROOT/'fit_v13'/'validation'/'V13_VAL_INTERFACE_01.xyz'
POT=Path(os.environ.get('V14_POT_DIR', '/workspace/-/models/ag_ti3sic2_qsgs_1to1_finer/potentials'))
LMP=os.environ.get('LAMMPS_BIN', '/workspace/venvs/podfit/bin/lmp')
atoms=read(SRC,format='extxyz')
data=OUT/'heldout_interface.data'
write(data,atoms,format='lammps-data',atom_style='atomic',units='metal',specorder=['Ag','C','Si','Ti'])
log=OUT/'log.lammps'; inp=OUT/'in.fire_nve'; min_dump=OUT/'after_fire.dump'; traj=OUT/'nve.dump'
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
pair_style hybrid/overlay eam tersoff morse 8.0 zbl 1.5 2.4
pair_coeff 1 1 eam {POT/'Ag_u3.eam'}
pair_coeff * * tersoff {POT/'TSCBOP.tersoff'} NULL C Si Ti
pair_coeff 1 2 morse 0.020000 7.804451 1.400000 8.0
pair_coeff 1 3 morse 0.020000 3.107766 1.400605 8.0
pair_coeff 1 4 morse 0.678644 1.155035 2.989514 8.0
pair_coeff 1 2 zbl 47 6
pair_coeff 1 3 zbl 47 14
pair_coeff 1 4 zbl 47 22
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
variable fmag atom sqrt(fx*fx+fy*fy+fz*fz)
compute fmax all reduce max v_fmag
thermo 25
thermo_style custom step atoms pe press c_fmax
min_style fire
min_modify dmax 0.05
minimize 1.0e-8 1.0e-3 2000 20000
write_dump all custom {min_dump} id type x y z fx fy fz modify sort id
write_data {OUT/'after_fire.data'}
velocity all create 300.0 24601 mom yes rot yes dist gaussian
fix int all nve
timestep 0.001
thermo_style custom step time temp pe ke etotal press c_fmax
dump tr all custom 10 {traj} id type x y z vx vy vz fx fy fz
dump_modify tr sort id
run 100
write_data {OUT/'after_nve.data'}
''')
env=os.environ.copy();env.update({'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
p=subprocess.run([LMP,'-in',str(inp),'-log',str(log),'-screen','none'],cwd=OUT,env=env,capture_output=True,text=True)
if p.returncode: raise RuntimeError((p.stdout+'\n'+p.stderr)[-3000:])

def read_dump(path): return read(path,index=':',format='lammps-dump-text',specorder=['Ag','C','Si','Ti'])
def geometry(path):
    rows=[]
    for k,a in enumerate(read_dump(path)):
        sy=np.asarray(a.get_chemical_symbols());dist=a.get_all_distances(mic=True);ff=np.linalg.norm(a.get_forces(),axis=1)
        pair={}
        for x,y in [('Ag','Ag'),('Ag','C'),('Ag','Si'),('Ag','Ti'),('Ti','C')]:
            ix=np.where(sy==x)[0];iy=np.where(sy==y)[0];sub=dist[np.ix_(ix,iy)]
            if x==y:sub=sub[np.triu_indices(len(ix),1)]
            pair[f'min_{x}_{y}_A']=float(sub.min())
        rows.append({'frame':k,'natoms':len(a),'Fmax_eV_A':float(ff.max()),**pair})
    return rows

text=log.read_text()
def parse_thermo(text, required, ncol):
    rows=[]; active=False
    for line in text.splitlines():
        if all(x in line for x in required):
            active=True; continue
        if not active: continue
        vals=line.split()
        if len(vals)==ncol and vals[0].isdigit():
            try: rows.append([float(x) for x in vals])
            except ValueError: pass
        elif rows:
            active=False
    return rows
iters=parse_thermo(text,['Step','Atoms','PotEng','Press','c_fmax'],5)
therm=parse_thermo(text,['Step','Time','Temp','PotEng','KinEng','TotEng','Press','c_fmax'],8)
energy_match=re.search(r'Energy initial, next-to-last, final =\s*([+\-0-9.eE]+)\s+([+\-0-9.eE]+)\s+([+\-0-9.eE]+)',text)
force_match=re.search(r'Force max component initial, final =\s*([+\-0-9.eE]+)\s+([+\-0-9.eE]+)',text)
iter_match=re.search(r'Iterations, force evaluations =\s*(\d+)\s+(\d+)',text)
result={'purpose':'small exploratory check only','source':str(SRC),'potential':'Ag_u3.eam + TSCBOP.tersoff + force-fitted Ag/skeleton Morse + cross-only ZBL 1.5-2.4 A',
        'minimization':{'thermo_rows':iters,'energy_initial_nexttolast_final_eV':list(map(float,energy_match.groups())) if energy_match else None,
                        'fmax_component_initial_final_eV_A':list(map(float,force_match.groups())) if force_match else None,
                        'iterations_force_evaluations':list(map(int,iter_match.groups())) if iter_match else None,
                        'final_force_geometry':geometry(min_dump)},
        'nve':{'dt_fs':1,'steps':100,'thermo_rows':therm,
               'temperature_initial_K':therm[0][2] if therm else None,'temperature_final_K':therm[-1][2] if therm else None,
               'total_energy_drift_eV_per_atom':(therm[-1][5]-therm[0][5])/len(atoms) if therm else None,
               'geometry':geometry(traj)},
        'limitations':['one 64-atom periodic ideal interface only','Morse cross terms are empirically fitted to five interface force labels','ZBL is a short-range guard and has not been independently fitted','not evidence for the 300k-atom geometry or long heating/cooling']}
(OUT/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
