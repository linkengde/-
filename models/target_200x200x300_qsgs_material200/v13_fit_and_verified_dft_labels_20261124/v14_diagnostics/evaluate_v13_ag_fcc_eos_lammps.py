#!/usr/bin/env python3
"""Evaluate the frozen V13 POD on fixed 256-atom FCC Ag lattice constants."""
import json, os, subprocess, tempfile
from pathlib import Path
import numpy as np
from ase.build import bulk
from ase.io import write

ROOT=Path(__file__).resolve().parent
FIT=ROOT/'fit_v13'
if not FIT.exists(): FIT=ROOT/'fit_v13_weight100'
LMP=os.environ.get('LAMMPS_EXECUTABLE','/workspace/venvs/podfit/bin/lmp')
OUT=ROOT/'v14_stage2_p0_ag_dft_diagnostics'/'ag_fcc_pod_static'
LATTICE_PARAMETERS_A=[4.09,3.90,3.60,3.30]
OUT.mkdir(parents=True,exist_ok=True)
rows=[]
for a in LATTICE_PARAMETERS_A:
    tag=f'a_{a:.2f}A'; data=OUT/f'{tag}.data'; inp=OUT/f'{tag}.in'; log=OUT/f'{tag}.log'
    atoms=bulk('Ag','fcc',a=a,cubic=True).repeat((4,4,4)); atoms.pbc=True
    write(data,atoms,format='lammps-data',atom_style='atomic',units='metal',specorder=['Ag','C','Si','Ti'])
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
pair_style pod
pair_coeff * * {FIT/'AgTiSiC_POD_v13.pod'} {FIT/'v13_coefficients.pod'} Ag C Si Ti
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
variable fmag atom sqrt(fx*fx+fy*fy+fz*fz)
compute fmax all reduce max v_fmag
thermo_style custom step atoms pe press c_fmax lx
thermo 1
run 0
''')
    env=os.environ.copy(); env.update({'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
    p=subprocess.run([LMP,'-in',str(inp),'-log',str(log),'-screen','none'],cwd=OUT,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    if p.returncode: raise SystemExit(p.stdout[-2000:])
    head=None; datum=None
    for line in log.read_text(errors='replace').splitlines():
        f=line.split()
        if f and f[0]=='Step' and 'PotEng' in f: head=f
        elif head and len(f)==len(head):
            try: datum=dict(zip(head,map(float,f)))
            except ValueError: pass
    if datum is None: raise RuntimeError(f'no LAMMPS thermo data for {tag}')
    rows.append({'a_A':a,'volume_per_atom_A3':a**3/4,'nearest_neighbor_A':a/np.sqrt(2),
                 'energy_eV_per_atom':datum['PotEng']/256,'pressure_bar':datum['Press'],
                 'fmax_eV_A':datum['c_fmax'],'natoms':256,'data':str(data),'input':str(inp),'log':str(log)})
summary={'potential':'V13 weight100 POD','purpose':'static diagnostic only','rows':rows}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
