#!/workspace/venvs/dft-labels/bin/python
"""GPAW single-point labels for a targeted Ti-C pair in full-model configurations."""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import numpy as np
from ase.io import read, write
from gpaw import FermiDirac, GPAW, Mixer, PW
from gpaw.mpi import world

HERE=Path(__file__).resolve().parent
if len(sys.argv)!=2 or sys.argv[1] not in {'initial','cold_10fs'}:
    raise SystemExit('usage: label_tic_pair.py initial|cold_10fs')
state=sys.argv[1]
source=HERE/'candidates'/f'tic_pair_{state}.extxyz'
outdir=HERE/'labels'
outdir.mkdir(exist_ok=True)
stem=f'tic_pair_{state}_pw350_gamma'
atoms=read(source)
ids=np.asarray(atoms.arrays['source_atom_id'],dtype=int)
ci=int(atoms.info['center_atom_id']); pi=int(atoms.info['partner_atom_id'])
indices={int(a):int(i) for i,a in enumerate(ids)}
if ci not in indices or pi not in indices or len(atoms)!=36:
    raise RuntimeError('Unexpected center/partner IDs or crop size')
start=time.time(); log=outdir/f'{stem}.txt'
atoms.calc=GPAW(mode=PW(350),xc='PBE',kpts=(1,1,1),symmetry='off',
    occupations=FermiDirac(0.10),mixer=Mixer(beta=0.05,nmaxold=8,weight=100),
    convergence={'density':1e-5,'eigenstates':1e-5,'energy':1e-5},maxiter=160,
    parallel={'domain':1,'band':max(1,world.size)},txt=str(log))
energy=float(atoms.get_potential_energy())
free_energy=float(atoms.get_potential_energy(force_consistent=True))
forces=np.asarray(atoms.get_forces(),dtype=float)
iterations=int(atoms.calc.get_number_of_iterations())
if not np.isfinite(forces).all() or not np.isfinite(energy):
    raise RuntimeError('DFT produced non-finite energy or forces')
center_forces={str(a):forces[indices[a]].tolist() for a in (ci,pi)}
podforces=np.asarray(atoms.arrays['POD_v12g_forces'],dtype=float)
errors={str(a):float(np.linalg.norm(forces[indices[a]]-podforces[indices[a]])) for a in (ci,pi)}
mask=np.asarray(atoms.arrays['is_DFT_force_target'],dtype=int)
atoms.arrays['REF_forces']=forces
atoms.info.update({'REF_energy':energy,'free_energy_eV':free_energy,
    'DFT_method':'GPAW PBE/PW350 Gamma, Fermi-Dirac sigma 0.10 eV',
    'DFT_center_forces_eV_A':json.dumps(center_forces),
    'DFT_POD_v12g_center_force_errors_eV_A':json.dumps(errors),
    'DFT_target_ids':json.dumps([ci,pi]),
    'label_status':'GPAW single-point; target force rows are Ti-C center pair; edge rows excluded from fitting',
    'model_data_role':'training_candidate; Si01 remains development holdout'})
atoms.calc=None
out=outdir/f'{stem}.extxyz'; write(out,atoms,format='extxyz')
logtext=log.read_text(errors='replace')
record={'status':'CONVERGED' if 'Converged in ' in logtext else 'SCF_CONVERGENCE_REVIEW_REQUIRED',
    'state':state,'center_atom_ids':[ci,pi],'elements':[atoms[indices[ci]].symbol,atoms[indices[pi]].symbol],
    'pair_distance_A':float(atoms.info['center_partner_distance_A']),'natoms':len(atoms),
    'energy_eV':energy,'free_energy_eV':free_energy,'center_DFT_forces_eV_A':center_forces,
    'center_POD_v12g_force_errors_eV_A':errors,'max_DFT_force_norm_eV_A':float(np.linalg.norm(forces,axis=1).max()),
    'DFT_target_rows':int(mask.sum()),'iterations':iterations,'elapsed_s':time.time()-start,
    'output':str(out),'log':str(log),'method':'GPAW PBE/PW350 Gamma, Fermi-Dirac sigma 0.10 eV'}
(outdir/f'{stem}.json').write_text(json.dumps(record,indent=2)+'\n')
if world.rank==0: print(json.dumps(record),flush=True)
