#!/workspace/venvs/dft-labels/bin/python
"""GPAW single-point labels for SNAP v12o 10 fs short-contact crops."""
from __future__ import annotations
import json,sys,time
from pathlib import Path
import numpy as np
from ase.io import read,write
from gpaw import FermiDirac,GPAW,Mixer,PW
from gpaw.mpi import world
HERE=Path(__file__).resolve().parent
KNOWN={
 'ag_si_430045_250182':([430045,250182],['Ag','Si']),
 'ag_si_455657_339883':([455657,339883],['Ag','Si']),
 'ag_c_454089_336261':([454089,336261],['Ag','C']),
 'ag_ti_412351_186486':([412351,186486],['Ag','Ti']),
}
name=sys.argv[1] if len(sys.argv)>1 else ''
if name not in KNOWN: raise SystemExit('Choose one of: '+', '.join(KNOWN))
target_ids,elements=KNOWN[name]; stem=f'{name}_t10fs_pw350_gamma'; atoms=read(HERE/'candidates'/f'{name}_t10fs.extxyz')
ids=np.asarray(atoms.arrays['source_atom_id'],int); target=np.asarray(atoms.arrays['is_DFT_force_target'],bool); ix={int(v):i for i,v in enumerate(ids)}
if len(atoms)>96 or target.sum()!=2 or sorted(ids[target].tolist())!=sorted(target_ids): raise RuntimeError('Bad crop target mapping')
if [atoms[ix[i]].symbol for i in target_ids]!=elements: raise RuntimeError('Wrong element map')
log=HERE/'labels'/f'{stem}.txt'; start=time.time()
atoms.calc=GPAW(mode=PW(350),xc='PBE',kpts=(1,1,1),symmetry='off',occupations=FermiDirac(0.10),mixer=Mixer(beta=0.05,nmaxold=8,weight=100),convergence={'density':1e-5,'eigenstates':1e-5,'energy':1e-5},maxiter=200,parallel={'domain':max(1,world.size),'band':1},txt=str(log))
energy=float(atoms.get_potential_energy()); free=float(atoms.get_potential_energy(force_consistent=True)); forces=np.asarray(atoms.get_forces(),float); iters=int(atoms.calc.get_number_of_iterations())
if not np.isfinite(forces).all() or not np.isfinite(energy): raise RuntimeError('Nonfinite DFT output')
model=np.asarray(atoms.arrays['SNAP_v12o_forces'],float); dft={str(i):forces[ix[i]].tolist() for i in target_ids}; err={str(i):float(np.linalg.norm(forces[ix[i]]-model[ix[i]])) for i in target_ids}
atoms.arrays['REF_forces']=forces; atoms.info.update({'REF_energy':energy,'free_energy_eV':free,'DFT_method':'GPAW PBE/PW350 Gamma, Fermi-Dirac sigma 0.10 eV','DFT_center_forces_eV_A':json.dumps(dft),'DFT_SNAP_v12o_center_force_errors_eV_A':json.dumps(err),'DFT_target_ids':json.dumps(target_ids),'label_status':'GPAW single-point; central pair force rows are targets; crop-edge rows are not fit targets','model_data_role':'v12p training candidate; validation and Si01 holdout unchanged'})
atoms.calc=None; out=HERE/'labels'/f'{stem}.extxyz'; write(out,atoms,format='extxyz')
text=log.read_text(errors='replace'); status='CONVERGED' if 'Converged in ' in text else 'SCF_CONVERGENCE_REVIEW_REQUIRED'
record={'status':status,'state':'SNAP_v12o_cold_10fs','candidate':name,'center_atom_ids':target_ids,'elements':elements,'center_pair_distance_A':float(atoms.info['center_pair_distance_A']),'natoms':len(atoms),'energy_eV':energy,'free_energy_eV':free,'center_DFT_forces_eV_A':dft,'center_SNAP_v12o_force_errors_eV_A':err,'max_DFT_force_norm_eV_A':float(np.linalg.norm(forces,axis=1).max()),'DFT_target_rows':2,'iterations':iters,'elapsed_s':time.time()-start,'output':str(out),'log':str(log),'method':'GPAW PBE/PW350 Gamma, Fermi-Dirac sigma 0.10 eV'}
(HERE/'labels'/f'{stem}.json').write_text(json.dumps(record,indent=2)+'\n')
if world.rank==0: print(json.dumps(record,indent=2),flush=True)
if status!='CONVERGED': raise SystemExit(2)
