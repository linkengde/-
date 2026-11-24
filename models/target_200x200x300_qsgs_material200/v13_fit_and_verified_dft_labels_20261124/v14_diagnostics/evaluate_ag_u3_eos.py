#!/usr/bin/env python3
"""Static Ag_u3.eam FCC energy-volume reference scan for comparison only."""
import json, os, subprocess
from pathlib import Path
import numpy as np
from ase.build import bulk
from ase.io import write

ROOT=Path(__file__).resolve().parent
EAM=Path('/workspace/-/models/ag_ti3sic2_qsgs_1to1_finer/potentials/Ag_u3.eam')
OUT=ROOT/'v14_stage2_p0_ag_dft_diagnostics'/'ag_u3_reference_eos'
OUT.mkdir(parents=True,exist_ok=True)
LMP=os.environ.get('LAMMPS_EXECUTABLE','/workspace/venvs/podfit/bin/lmp')
rows=[]
for a in [4.09,3.90,3.60,3.30]:
    tag=f'a_{a:.2f}A'; atoms=bulk('Ag','fcc',a=a,cubic=True).repeat((4,4,4)); atoms.pbc=True
    data=OUT/f'{tag}.data'; inp=OUT/f'{tag}.in'; log=OUT/f'{tag}.log'
    write(data,atoms,format='lammps-data',atom_style='atomic',units='metal',specorder=['Ag'])
    inp.write_text(f'''clear
units metal
atom_style atomic
boundary p p p
atom_modify sort 0 0
read_data {data}
mass 1 107.8682
pair_style eam
pair_coeff * * {EAM}
thermo_style custom step atoms pe press lx
thermo 1
run 0
''')
    env=os.environ.copy(); env.update({'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
    p=subprocess.run([LMP,'-in',str(inp),'-log',str(log),'-screen','none'],cwd=OUT,env=env,capture_output=True,text=True)
    if p.returncode: raise SystemExit(f'{tag}: {p.stdout[-1000:]} {p.stderr[-1000:]}')
    header=None; value=None
    for line in log.read_text(errors='replace').splitlines():
        v=line.split()
        if v and v[0]=='Step' and 'PotEng' in v: header=v
        elif header and len(v)==len(header):
            try: value=dict(zip(header,map(float,v)))
            except ValueError: pass
    if value is None: raise RuntimeError(f'No thermo row for {tag}')
    rows.append({'lattice_A':a,'volume_per_atom_A3':a**3/4,'nearest_neighbor_A':a/np.sqrt(2),
      'energy_eV_per_atom':value['PotEng']/256,'pressure_bar':value['Press'],'natoms':256})
ref=rows[0]['energy_eV_per_atom']
for row in rows: row['delta_E_vs_a4p09_eV_per_atom']=row['energy_eV_per_atom']-ref
summary={'potential':'Ag_u3.eam','citation':'Foiles, Baskes, Daw, Phys. Rev. B 33, 7983 (1986)',
 'purpose':'pure-Ag reference/sampling diagnostic only; not the Ag-Ti-Si-C final potential', 'rows':rows}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
