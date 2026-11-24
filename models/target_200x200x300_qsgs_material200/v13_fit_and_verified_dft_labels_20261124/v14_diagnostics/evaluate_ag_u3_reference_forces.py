#!/usr/bin/env python3
"""Compare archived Ag_u3.eam against existing Ag DFT forces."""
import json, os, subprocess
from pathlib import Path
import numpy as np
from ase.io import read, write

ROOT=Path(__file__).resolve().parent
OLD=Path('/workspace/potential-fit/round7_rebuild/labels')
EAM=Path('/workspace/-/models/ag_ti3sic2_qsgs_1to1_finer/potentials/Ag_u3.eam')
OUT=ROOT/'v14_stage2_p0_ag_dft_diagnostics'/'ag_u3_reference_force_checks'
OUT.mkdir(parents=True,exist_ok=True)
LMP=os.environ.get('LAMMPS_EXECUTABLE','/workspace/venvs/podfit/bin/lmp')
items=[]
for src in sorted(OLD.glob('ag_cooling_*_step*_pw350.extxyz')):
    items.append((src,src.with_suffix('.json')))
candidate=ROOT/'v14_stage2_p0_ag_dft_diagnostics'/'ag_liquid_precollapse'/'label.extxyz'
if candidate.exists():
    items.append((candidate,candidate.with_name('label.json')))
results=[]
for src,meta_path in items:
    atoms=read(src,format='extxyz')
    f_ref=np.asarray(atoms.get_forces(),float)
    tag=src.stem
    data=OUT/f'{tag}.data'; inp=OUT/f'{tag}.in'; log=OUT/f'{tag}.log'; dump=OUT/f'{tag}.dump'
    write(data,atoms,format='lammps-data',atom_style='atomic',units='metal',specorder=['Ag','C','Si','Ti'])
    inp.write_text(f'''clear
units metal
atom_style atomic
boundary p p p
atom_modify sort 0 0
read_data {data}
mass 1 107.8682
pair_style eam
pair_coeff * * {EAM}
thermo_style custom step atoms pe press
thermo 1
dump d all custom 1 {dump} id type x y z fx fy fz
dump_modify d sort id
run 0
''')
    env=os.environ.copy()
    env.update({'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
    proc=subprocess.run([LMP,'-in',str(inp),'-log',str(log),'-screen','none'],cwd=OUT,env=env,capture_output=True,text=True)
    if proc.returncode:
        raise SystemExit(f'{tag}: {proc.stdout[-1000:]} {proc.stderr[-1000:]}')
    lines=dump.read_text().splitlines()
    k=lines.index('ITEM: ATOMS id type x y z fx fy fz'); n=int(lines[3])
    vals=np.array([[float(x) for x in line.split()] for line in lines[k+1:k+1+n]])
    vals=vals[np.argsort(vals[:,0].astype(int))]
    f_eam=vals[:,5:8]; err=f_eam-f_ref; errnorm=np.linalg.norm(err,axis=1)
    distances=atoms.get_all_distances(mic=True); distances[distances<1e-8]=np.inf
    meta=json.loads(meta_path.read_text()) if meta_path.exists() else {}
    row={'label_file':str(src),'natoms':len(atoms),
      'minimum_pair_distance_A':float(distances.min()),
      'dft_fmax_eV_A':float(np.linalg.norm(f_ref,axis=1).max()),
      'eam_fmax_eV_A':float(np.linalg.norm(f_eam,axis=1).max()),
      'force_rmse_component_eV_A':float(np.sqrt(np.mean(err**2))),
      'force_rmse_vector_norm_eV_A':float(np.sqrt(np.mean(np.sum(err**2,axis=1)))),
      'force_error_p95_eV_A':float(np.percentile(errnorm,95)),
      'force_error_max_eV_A':float(errnorm.max()),
      'label_protocol':meta.get('method',meta.get('dft_protocol'))}
    if 'precollapse' in str(src):
        i,j=21,26
        rij=atoms.positions[j]-atoms.positions[i]
        rij-=np.round(rij/atoms.cell.lengths())*atoms.cell.lengths()
        rhat=rij/np.linalg.norm(rij)
        row['short_pair_radial_outward_force_eV_A']={
          'DFT_on_22':float(np.dot(f_ref[i],-rhat)),
          'DFT_on_27':float(np.dot(f_ref[j],rhat)),
          'Ag_u3_EAM_on_22':float(np.dot(f_eam[i],-rhat)),
          'Ag_u3_EAM_on_27':float(np.dot(f_eam[j],rhat))}
    results.append(row)
summary={'reference_potential':str(EAM),
 'citation':'Foiles, Baskes, Daw, Phys. Rev. B 33, 7983 (1986), Ag universal-3 EAM',
 'purpose':'sampling/reference check only; not an Ag-Ti-Si-C potential or a production recommendation',
 'results':results}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
