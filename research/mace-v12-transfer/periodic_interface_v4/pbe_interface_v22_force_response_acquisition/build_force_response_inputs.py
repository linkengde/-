#!/usr/bin/env python3
"""Build and geometry-screen six Ag-layer separation-response DFT inputs."""
from __future__ import annotations
import hashlib,json
from pathlib import Path
import numpy as np
from ase.io import read,write
from geometry_utils import compare_geometry,extxyz_geometry_frames

HERE=Path(__file__).resolve().parent
REPO=Path('/workspace/jeitschko-main-run')
V21=REPO/'research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v21_jeitschko_interface_acquisition'
TRAIN43=V21/'training_candidate_v19_train43_20261029_A/dataset/train.extxyz'
DEV3=V21/'training_candidate_v19_train43_20261029_A/dataset/valid.extxyz'
OUT=HERE/'inputs'
SHIFTS={'AgC':(-0.20,0.20),'AgSi':(-0.20,0.20),'AgTi':(-0.20,0.20)}
BASES={
 'AgC':V21/'inputs/AgC_JN1967_registryB_v21_acq01.extxyz',
 'AgSi':V21/'inputs/AgSi_JN1967_registryB_v21_acq01.extxyz',
 'AgTi':V21/'inputs/AgTi_JN1967_registryB_v21_acq01.extxyz',
}
SOURCE_CIF_SHA='26960fc12d2de584c6a13ba5fe0d15713c07f697f455dbe29791b747f07359a1'
SOURCE_PAPER_SHA='6765f99b3cc83daf3e48e4674dbeba109e5fffef8fa66cd53e266733ac4a4a63'

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
 if OUT.exists() and any(OUT.iterdir()): raise SystemExit(f'refusing to overwrite nonempty {OUT}')
 OUT.mkdir(parents=True,exist_ok=True)
 train_geom=extxyz_geometry_frames(TRAIN43)
 dev_geom=extxyz_geometry_frames(DEV3)
 base_geom={k:extxyz_geometry_frames(v)[0] for k,v in BASES.items()}
 parent_record={k:read(v) for k,v in BASES.items()}
 rows=[]
 for iface,(minus,plus) in SHIFTS.items():
  for offset in (minus,plus):
   sign='m' if offset<0 else 'p'; mag=f'{abs(offset):.2f}'.replace('.','p')
   label=f'{iface}_JN1967_registryB_gap_{sign}{mag}_sigma0p10_k6x6x1_v22'
   parent=parent_record[iface].copy()
   if parent.get_chemical_formula()!='C16Ag4Si8Ti24' or len(parent)!=52 or parent.pbc.tolist()!=[True,True,False]:
    raise SystemExit(f'unexpected parent interface geometry: {iface}')
   syms=np.asarray(parent.get_chemical_symbols()); ag=np.flatnonzero(syms=='Ag')
   if len(ag)!=4: raise SystemExit(f'expected Ag monolayer of four atoms: {iface}')
   parent.positions[ag,2]+=offset
   parent.info.update({
    'config_type':label,'candidate_role':'train_acquisition_candidate_pending_audit',
    'dataset_role':'train_acquisition_candidate_pending_audit','v22_parent_geometry':f'{iface}_JN1967_registryB_v21_acq01',
    'source_family':'Jeitschko_Nowotny_1967_single_crystal','source_doi':'10.1007/BF00899949',
    'source_cif_sha256':SOURCE_CIF_SHA,'source_paper_sha256':SOURCE_PAPER_SHA,
    'geometry_perturbation':'Rigid Ag monolayer translation along +c relative to fixed Ti3SiC2 slab',
    'Ag_layer_delta_z_A':float(offset),'fixed_geometry_unlabeled':True,'full_occupancy_model_assumption':True,
   })
   for k in ('energy','REF_energy','PW_PBE_energy_eV','PW_PBE_free_energy_eV'): parent.info.pop(k,None)
   for k in ('forces','REF_forces','PW_PBE_forces'): parent.arrays.pop(k,None)
   target=OUT/f'{label}.extxyz'; write(target,parent,format='extxyz'); atoms=read(target)
   if len(atoms)!=52 or atoms.get_chemical_formula()!='C16Ag4Si8Ti24': raise SystemExit(f'bad output formula/count: {label}')
   if any(k in atoms.info for k in ('energy','REF_energy','PW_PBE_energy_eV')) or any(k in atoms.arrays for k in ('forces','REF_forces','PW_PBE_forces')):
    raise SystemExit(f'input unexpectedly contains a label: {label}')
   d=atoms.get_all_distances(mic=True); np.fill_diagonal(d,np.inf)
   minimum=float(d.min());
   if minimum<2.0: raise SystemExit(f'implausibly short distance {minimum:.3f} A in {label}')
   central=np.flatnonzero(atoms.arrays['central_pair'])
   if len(central)!=2: raise SystemExit(f'missing marked contact pair: {label}')
   contact=float(atoms.get_distance(int(central[0]),int(central[1]),mic=True))
   screen_groups={}
   for name,refs in [('train43_geometry_only',train_geom),('dev3_geometry_only',dev_geom)]:
    comparisons=[]
    for ref in refs:
     candidate=atoms.copy(); candidate.set_cell(ref.cell,scale_atoms=True)
     cmp=compare_geometry(candidate,ref)
     if cmp is not None: comparisons.append(cmp)
    screen_groups[name]={
     'frames_read_geometry_only':len(refs),'same_composition_normalized_comparisons':len(comparisons),
     'exact_hits':sum(x['exact_duplicate'] for x in comparisons),
     'near_hits':sum(x['near_duplicate'] for x in comparisons),
     'comparisons':comparisons,
    }
   if screen_groups['dev3_geometry_only']['exact_hits'] or screen_groups['dev3_geometry_only']['near_hits']:
    raise SystemExit(f'candidate overlaps frozen dev3 geometry screen: {label}; inspect geometry_audit.json')
   if screen_groups['train43_geometry_only']['exact_hits']:
    raise SystemExit(f'exact duplicate against train43: {label}')
   rows.append({
    'label':label,'interface':iface,'input':str(target.relative_to(HERE)),'input_sha256':sha(target),
    'parent_input':str(BASES[iface].relative_to(REPO)),'parent_input_sha256':sha(BASES[iface]),
    'atoms':len(atoms),'formula':atoms.get_chemical_formula(),'cell_A':atoms.cell.array.tolist(),
    'pbc':atoms.pbc.tolist(),'Ag_layer_delta_z_A':float(offset),'marked_contact_pair_distance_A':contact,
    'minimum_periodic_pair_distance_A':minimum,'role':'train_acquisition_candidate_pending_audit',
    'launch_enabled':True,'geometry_screens':screen_groups,
   })
 # Paired perturbations must bracket each original parent geometry without changing the substrate.
 for iface in SHIFTS:
  pair=[x for x in rows if x['interface']==iface]
  if len(pair)!=2 or sorted(x['Ag_layer_delta_z_A'] for x in pair)!=[-0.2,0.2]: raise SystemExit(f'incomplete displacement pair for {iface}')
  if any(abs(x['minimum_periodic_pair_distance_A']-2.0)<1e-9 for x in pair): raise SystemExit(f'minimum-distance boundary issue for {iface}')
 method={
  'xc':'PBE','basis':'plane wave','cutoff_eV':500,'smearing_eV':0.10,'pbc':[True,True,False],
  'poissonsolver':{'dipolelayer':'xy'},'candidate_kpts':[6,6,1],'mpi_ranks':4,
  'maximum_concurrent_jobs':1,'GPAW_version':'26.7.0','ASE_version':'3.29.0','PAW_data_version':'1.2.1',
  'PAW_setup_sha256':{'Ag.PBE.gz':'bff00ac79a7fa6e8c3a71c8b34d3a9c488e6bf1c00b3e9c4c91c231551145e97','Ti.PBE.gz':'bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f','Si.PBE.gz':'2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0','C.PBE.gz':'fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983'},
  'scf_thresholds':{'energy':1e-5,'density':1e-5,'eigenstates':1e-5},'mixer':{'type':'Pulay','beta':0.05,'nmaxold':8,'weight':100},
 }
 manifest={
  'task_id':'user_directed_Jeitschko1967_AgTi3SiC2_force_response_acquisition_v22',
  'status':'geometry_screen_pass_prelaunch_resource_preflight_required',
  'source':{'family':'Jeitschko_Nowotny_1967_single_crystal','doi':'10.1007/BF00899949','source_cif_sha256':SOURCE_CIF_SHA,'source_paper_sha256':SOURCE_PAPER_SHA},
  'reference_datasets':{'train43_sha256':sha(TRAIN43),'dev3_sha256':sha(DEV3),'reference_read_scope':'geometry columns only; labels not inspected'},
  'method_candidate':method,
  'runtime_guard':{'cgroup_memory_limit_bytes':34359738368,'stop_threshold_bytes':32212254720,'minimum_free_disk_bytes':5368709120,'automatic_retries':0,'maximum_concurrent_jobs':1},
  'role_policy':'All six results are experiment-only train acquisition candidates; no official split/model change. No independent validation role assigned.',
  'geometry_perturbation_policy':'Shift the four Ag atoms together by ±0.20 Å along surface normal, keeping substrate fixed. This samples force response without using dev3 coordinates or labels.',
  'records':rows,
  'preflight_jobs':[{'label':r['label'],'geometry_label':f"{r['interface']}_JN1967_registryB_v21_acq01_Agdz{r['Ag_layer_delta_z_A']:+.2f}",'input':r['input'],'input_sha256':r['input_sha256'],'kpts':[6,6,1],'smearing_eV':0.10,'role':'train_acquisition_candidate_pending_audit','launch_enabled':True} for r in rows],
  'launch_enabled':True,'training_enabled':False,'dataset_integration_enabled':False,
  'test_or_sealed_data_opened':False,'labels_or_forces_read':False,
  'limitations':['All candidates remain one Jeitschko 1967 source family.','This is not a universal cross-archive supercell-equivalence audit.','The v21 slab mesh control supports the selected 6x6x1 mesh for the original AgSi geometry only; no new geometry-specific mesh control is included.'],
 }
 (HERE/'input_manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n')
 audit={'status':'GEOMETRY_SCREEN_PASS_RESOURCE_PREFLIGHT_REQUIRED','records':rows,'labels_or_forces_read':False,'test_or_sealed_data_opened':False}
 (HERE/'geometry_audit.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False)+'\n')
 print(json.dumps({'status':manifest['status'],'candidate_count':len(rows),'rows':[{'label':r['label'],'contact_A':r['marked_contact_pair_distance_A'],'minimum_A':r['minimum_periodic_pair_distance_A'],'dev3_near_hits':r['geometry_screens']['dev3_geometry_only']['near_hits']} for r in rows]},indent=2))

if __name__=='__main__': main()
