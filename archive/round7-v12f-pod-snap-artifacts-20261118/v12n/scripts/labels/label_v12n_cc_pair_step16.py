#!/workspace/venvs/dft-labels/bin/python
"""GPAW PBE/PW350 label for the v12m cold C-C collapse precursor."""
from __future__ import annotations
import json,time
from pathlib import Path
import numpy as np
from ase.io import read,write
from gpaw import FermiDirac,GPAW,Mixer,PW
from gpaw.mpi import world
HERE=Path(__file__).resolve().parent
SOURCE=HERE/'candidates/cc_pair_v12m_t26fs_step16.extxyz'
OUTDIR=HERE/'labels'; STEM='cc_pair_v12m_t26fs_step16_pw350_gamma'; TARGET_IDS=[128985,134363]
OUTDIR.mkdir(exist_ok=True)
atoms=read(SOURCE); ids=np.asarray(atoms.arrays['source_atom_id'],int); target=np.asarray(atoms.arrays['is_DFT_force_target'],bool); index={int(a):i for i,a in enumerate(ids)}
if len(atoms)>64 or target.sum()!=2 or any(a not in index for a in TARGET_IDS): raise RuntimeError('Unexpected C-C crop/targets')
if [atoms[index[a]].symbol for a in TARGET_IDS]!=['C','C']: raise RuntimeError('C-C target mapping mismatch')
log=OUTDIR/f'{STEM}.txt'; start=time.time()
atoms.calc=GPAW(mode=PW(350),xc='PBE',kpts=(1,1,1),symmetry='off',occupations=FermiDirac(0.10),mixer=Mixer(beta=0.05,nmaxold=8,weight=100),convergence={'density':1e-5,'eigenstates':1e-5,'energy':1e-5},maxiter=160,parallel={'domain':max(1,world.size),'band':1},txt=str(log))
energy=float(atoms.get_potential_energy()); free_energy=float(atoms.get_potential_energy(force_consistent=True)); forces=np.asarray(atoms.get_forces(),float); iterations=int(atoms.calc.get_number_of_iterations())
if not np.isfinite(forces).all() or not np.isfinite(energy): raise RuntimeError('Nonfinite DFT result')
model_forces=np.asarray(atoms.arrays['POD_v12m_forces'],float); center_forces={str(a):forces[index[a]].tolist() for a in TARGET_IDS}; center_errors={str(a):float(np.linalg.norm(forces[index[a]]-model_forces[index[a]])) for a in TARGET_IDS}
atoms.arrays['REF_forces']=forces; atoms.info.update({'REF_energy':energy,'free_energy_eV':free_energy,'DFT_method':'GPAW PBE/PW350 Gamma, Fermi-Dirac sigma 0.10 eV','DFT_center_forces_eV_A':json.dumps(center_forces),'DFT_POD_v12m_center_force_errors_eV_A':json.dumps(center_errors),'DFT_target_ids':json.dumps(TARGET_IDS),'label_status':'GPAW single-point; only two central C forces are DFT targets; edge rows excluded from DFT targets','model_data_role':'v12n training candidate; Si01 remains development holdout'}); atoms.calc=None
output=OUTDIR/f'{STEM}.extxyz'; write(output,atoms,format='extxyz'); log_text=log.read_text(errors='replace')
record={'status':'CONVERGED' if 'Converged in ' in log_text else 'SCF_CONVERGENCE_REVIEW_REQUIRED','state':'v12m_POD_cold_C-C_precursor_t26fs','center_atom_ids':TARGET_IDS,'elements':['C','C'],'center_pair_distance_A':float(atoms.info['center_pair_distance_A']),'natoms':len(atoms),'energy_eV':energy,'free_energy_eV':free_energy,'center_DFT_forces_eV_A':center_forces,'center_POD_v12m_force_errors_eV_A':center_errors,'max_DFT_force_norm_eV_A':float(np.linalg.norm(forces,axis=1).max()),'DFT_target_rows':int(target.sum()),'iterations':iterations,'elapsed_s':time.time()-start,'output':str(output),'log':str(log),'method':'GPAW PBE/PW350 Gamma, Fermi-Dirac sigma 0.10 eV'}
(OUTDIR/f'{STEM}.json').write_text(json.dumps(record,indent=2)+'\n')
if world.rank==0: print(json.dumps(record,indent=2),flush=True)
