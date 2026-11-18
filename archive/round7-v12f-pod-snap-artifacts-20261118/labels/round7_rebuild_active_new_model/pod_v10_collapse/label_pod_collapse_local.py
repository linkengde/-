#!/workspace/venvs/dft-labels/bin/python
"""PBE/PW350 Gamma labels for the step-50 POD-v10 collapse crops."""
from pathlib import Path
import json
import sys
import time
import numpy as np
from ase.io import read, write
from gpaw import GPAW, FermiDirac, Mixer, PW
from gpaw.mpi import world

HERE=Path(__file__).resolve().parent
SOURCE=HERE/'candidate_geometries.extxyz'
OUTDIR=HERE/'labels'
OUTDIR.mkdir(exist_ok=True)
idx=int(sys.argv[1])
atoms=read(SOURCE,index=idx)
start=time.time()
log=OUTDIR/f'pod_collapse_local_{idx:02d}_pw350_gamma.txt'
atoms.calc=GPAW(mode=PW(350),xc='PBE',kpts=(1,1,1),symmetry='off',
    occupations=FermiDirac(0.10),mixer=Mixer(beta=0.05,nmaxold=8,weight=100),
    convergence={'density':1e-5,'eigenstates':1e-5,'energy':1e-5},maxiter=160,
    parallel={'domain':max(1,world.size)},txt=str(log))
energy=float(atoms.get_potential_energy())
free_energy=float(atoms.get_potential_energy(force_consistent=True))
forces=np.asarray(atoms.get_forces())
iterations=int(atoms.calc.get_number_of_iterations())
ids=np.asarray(atoms.arrays['source_atom_id'],dtype=int)
center_id=int(atoms.info['center_atom_id'])
partner_id=int(atoms.info['partner_atom_id'])
ci=int(np.flatnonzero(ids==center_id)[0])
pi=int(np.flatnonzero(ids==partner_id)[0])
atoms.info.update({
    'REF_energy':energy,'REF_energy_unshifted_PBE_PW350_Gamma_eV':energy,
    'free_energy_eV':free_energy,
    'label_status':'GPAW PBE/PW350 Gamma single-point labels for a POD-v10 FIRE-collapse local crop',
    'DFT_method':'GPAW PBE/PW350 Gamma, Fermi-Dirac sigma 0.10 eV',
    'DFT_cutoff_eV':350,'DFT_kpts':[1,1,1],
    'DFT_center_force_eV_A':forces[ci].tolist(),
    'DFT_center_force_norm_eV_A':float(np.linalg.norm(forces[ci])),
    'DFT_partner_force_eV_A':forces[pi].tolist(),
    'DFT_partner_force_norm_eV_A':float(np.linalg.norm(forces[pi])),
    'DFT_edge_force_note':'Finite crop; do not use edge forces as full-periodic targets without separate weighting.',
    'DFT_energy_note':'Finite-crop raw energy; do not mix directly with full-periodic training energies.',
    'model_data_role':'training-candidate only; never use as independent blind test.',
})
atoms.arrays['REF_forces']=forces
atoms.calc=None
outfile=OUTDIR/f'pod_collapse_local_{idx:02d}_pw350_gamma.extxyz'
write(outfile,atoms,format='extxyz')
record={'frame_index':idx,'source_step':int(atoms.info['source_step']),
    'center_atom_id':center_id,'center_species':atoms.info['center_species'],
    'partner_atom_id':partner_id,'partner_species':atoms.info['partner_species'],
    'pair_distance_A':float(atoms.info['center_partner_distance_A']),'formula':atoms.get_chemical_formula(),
    'natoms':len(atoms),'energy_eV':energy,'free_energy_eV':free_energy,
    'center_force_norm_eV_A':float(np.linalg.norm(forces[ci])),
    'partner_force_norm_eV_A':float(np.linalg.norm(forces[pi])),
    'max_force_norm_eV_A':float(np.linalg.norm(forces,axis=1).max()),
    'iterations':iterations,'elapsed_s':time.time()-start,
    'output':str(outfile),'log':str(log),'status':'computed; audit pending'}
(OUTDIR/f'pod_collapse_local_{idx:02d}_pw350_gamma.json').write_text(json.dumps(record,indent=2,ensure_ascii=False)+'\n')
if world.rank==0:
    print(json.dumps(record,ensure_ascii=False),flush=True)
