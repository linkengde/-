#!/usr/bin/env python3
"""Create geometry-only DFT input candidates from the archived 48-atom parent.
This script does not run DFT, assign dataset roles, or optimize structures.
"""
from pathlib import Path
import argparse, hashlib, json
import numpy as np
from ase.io import read, write

ap=argparse.ArgumentParser()
ap.add_argument('--parent',required=True)
ap.add_argument('--output-dir',required=True)
ap.add_argument('--environment-preflight', help='Optional JSON snapshot from the E runtime preflight')
args=ap.parse_args()
parent=Path(args.parent)
out=Path(args.output_dir); out.mkdir(parents=True,exist_ok=True)
raw=parent.read_bytes(); parent_sha=hashlib.sha256(raw).hexdigest()
base=read(str(parent))
expected={'Ag':0,'Ti':24,'Si':8,'C':16}
counts={s:base.get_chemical_symbols().count(s) for s in ['Ag','Ti','Si','C']}
assert len(base)==48 and counts==expected
assert base.pbc.tolist()==[True,True,True]
ids=base.arrays.get('local_new_id')
assert ids is not None and len(ids)==len(base)
si_indices=[i for i,s in enumerate(base.get_chemical_symbols()) if s=='Si']
si_idx=next(i for i in si_indices if int(ids[i])==7)

variants=[]
def make(name,transform,desc):
    a=base.copy()
    # Drop proposal annotations; preserve source atom IDs and site-kind array only.
    a.info={}
    a.positions=np.asarray(a.positions) @ transform.T
    a.set_cell(np.asarray(a.cell) @ transform.T,scale_atoms=False)
    a.info.update({
        'structure_id':name,
        'source_parent_sha256':parent_sha,
        'source_parent_id':'Ti3SiC2_48atom_bulk_baseline_PARENT',
        'transformation':desc,
        'dataset_role':'UNASSIGNED_WAITING_C_FREEZE',
        'dft_label_status':'PENDING_NO_DFT_RUN',
    })
    path=out/(name+'.extxyz')
    write(path,a,format='extxyz')
    geom={'numbers':a.numbers.astype(int).tolist(),'positions_A':np.round(a.positions,8).tolist(),'cell_A':np.round(a.cell.array,10).tolist(),'pbc':a.pbc.astype(bool).tolist()}
    geom_sha=hashlib.sha256(json.dumps(geom,separators=(',',':')).encode()).hexdigest()
    variants.append({'task_id':name,'structure_file':path.name,'structure_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'ordered_geometry_sha256_rounded_1e-8A':geom_sha,'atom_count':len(a),'composition':counts,'cell_A':a.cell.array.tolist(),'pbc':a.pbc.tolist(),'transformation':desc,'data_role':'UNASSIGNED_WAITING_C_FREEZE','label_status':'PENDING_NO_DFT_RUN'})

I=np.eye(3)
for scale,suffix in [(0.98,'m02'),(1.02,'p02')]:
    f=scale**(1/3)
    make(f'TSC48_VOL_{suffix}',f*I,f'isotropic volume factor {scale:.4f}; linear scale {f:.12f}; atoms and cell uniformly scaled; fixed-geometry single point after input freeze')
# Basal-plane engineering shear gamma = +/-0.005, simple shear x' = x + gamma*y.
for gamma,suffix in [(-0.005,'m005'),(0.005,'p005')]:
    F=np.array([[1.,gamma,0.],[0.,1.,0.],[0.,0.,1.]])
    make(f'TSC48_BASAL_XY_SHEAR_{suffix}',F,f'basal xy simple shear gamma={gamma:+.6f}; Cartesian map x_prime=x+gamma*y; fixed-geometry single point after input freeze')
for dz,suffix in [(-0.02,'m020'),(0.02,'p020')]:
    a=base.copy(); a.info={}; a.positions[si_idx,2]+=dz
    a.info.update({'structure_id':f'TSC48_SI_ID7_Z_{suffix}','source_parent_sha256':parent_sha,'source_parent_id':'Ti3SiC2_48atom_bulk_baseline_PARENT','transformation':f'local_new_id=7 (Si, spacegroup_kinds=2) Cartesian z displacement {dz:+.3f} A; cell fixed','dataset_role':'UNASSIGNED_WAITING_C_FREEZE','dft_label_status':'PENDING_NO_DFT_RUN'})
    path=out/(f'TSC48_SI_ID7_Z_{suffix}.extxyz'); write(path,a,format='extxyz')
    geom={'numbers':a.numbers.astype(int).tolist(),'positions_A':np.round(a.positions,8).tolist(),'cell_A':np.round(a.cell.array,10).tolist(),'pbc':a.pbc.astype(bool).tolist()}
    geom_sha=hashlib.sha256(json.dumps(geom,separators=(',',':')).encode()).hexdigest()
    variants.append({'task_id':f'TSC48_SI_ID7_Z_{suffix}','structure_file':path.name,'structure_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'ordered_geometry_sha256_rounded_1e-8A':geom_sha,'atom_count':len(a),'composition':counts,'cell_A':a.cell.array.tolist(),'pbc':a.pbc.tolist(),'transformation':f'local_new_id=7 (Si, spacegroup_kinds=2) Cartesian z displacement {dz:+.3f} A; cell fixed','data_role':'UNASSIGNED_WAITING_C_FREEZE','label_status':'PENDING_NO_DFT_RUN'})

manifest={
 'package_status':'INPUT_CANDIDATES_ONLY_NOT_DFT_RESULTS',
 'source_commit':'499ddf0870235594444abad780582ea4fcc89d61',
 'source_path':'research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_pure_phase_controls/inputs/Ti3SiC2_baseline_PROPOSAL.extxyz',
 'source_parent_file_sha256':parent_sha,
 'matched_existing_dft_geometry':'Parent coordinates/order/cell match Ti3SiC2_baseline_k10x10x4_v19 output (SHA256 75b083e225b8b20bbdde181c8e1847bae4b0052a854c6e9ee48885c940166415); the same baseline geometry is independently archived at k6x6x4 (output SHA256 29e9ac254984b4a4ab6a22cbac7043c2cf9a826eb498422f1860a5798e72cee4).',
 'source_geometry_check':{'atom_count':48,'composition':counts,'pbc':[True,True,True],'cell_A':base.cell.array.tolist(),'max_coord_difference_from_k10x10x4_label_A':0.0},
 'role_freeze_source':'E2 role_freeze.json; TRAIN/DEV/TEST are empty; these six points are kept in the existing 48-atom parent family and remain unassigned until C freezes a manifest',
 'proposed_calculation':{'program':'GPAW','version':'26.7.0','ASE_version':'3.29.0','xc':'PBE','basis':'plane waves','cutoff_eV':500,'PAW_data_version':'1.2.1','PAW_SHA256':{'Ti.PBE.gz':'bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f','Si.PBE.gz':'2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0','C.PBE.gz':'fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983'},'occupations':'Fermi-Dirac 0.10 eV','kpoints_proposal':[10,10,4],'kpoint_centering':'NOT_VERIFIED_IN_ARCHIVED_INPUT; C must freeze','scf_thresholds':{'energy':1e-5,'density':1e-5,'eigenstates':1e-5},'mixer':'Pulay beta=0.05, nmaxold=8, weight=100','maxiter_proposal':160,'geometry':'fixed atoms; fixed cell for each candidate; no relaxation','requested_outputs':['native total energy','force-consistent free energy','per-atom forces','stress tensor','SCF convergence and steps','input/output SHA256'],'mpi_omp':'NOT_FROZEN_BY_C'},
 'execution_authorized':False,
 'blockers':['C role and exact-input freeze not issued','no C DFT execution authorization found','no dedicated DFT work-root configured','k-point centering, MPI ranks, wall-time cap and checkpoint policy not frozen'],
 'variants':variants}
if args.environment_preflight:
    preflight_path=Path(args.environment_preflight)
    manifest['environment_preflight_file']=preflight_path.name
    manifest['environment_preflight_sha256']=hashlib.sha256(preflight_path.read_bytes()).hexdigest()
    manifest['environment_preflight']=json.loads(preflight_path.read_text())
(out/'DFT_INPUT_CANDIDATE_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps({'parent_sha256':parent_sha,'variants':len(variants),'ids':[x['task_id'] for x in variants]},indent=2))
