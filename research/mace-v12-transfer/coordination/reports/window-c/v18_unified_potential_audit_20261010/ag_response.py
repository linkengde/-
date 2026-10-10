#!/usr/bin/env python3
"""Frozen MACE energy/force response on existing, same-mesh solid Ag DFT controls."""
import csv,hashlib,importlib.util,json,os,tempfile
from pathlib import Path
import numpy as np,torch
from ase.io import read
from mace.calculators import MACECalculator
ROOT=Path(os.environ.get('CLOUD_C_REPO_ROOT',Path(__file__).resolve().parents[6])).resolve();BASE=ROOT/'research/mace-v12-transfer';OUT=Path(os.environ.get('CLOUD_C_AUDIT_OUTPUT',Path(tempfile.gettempdir())/'cloud-c-v18-audit')).resolve(); OUT.mkdir(parents=True,exist_ok=True)
torch.set_num_threads(4);torch.set_num_interop_threads(1)
purepath=BASE/'coordination/reports/window-a/v19_V18_pure_phase_force_screen/infer_pure_phase.py'
spec=importlib.util.spec_from_file_location('pure_helpers',purepath);pure=importlib.util.module_from_spec(spec);spec.loader.exec_module(pure)
folder=BASE/'periodic_interface_v4/pbe_interface_v19_pure_phase_controls/calculations'
groups={
 'displacement_k6':['Ag_baseline_k6x6x6_v19','Ag_displacement_m01_k6x6x6_v19','Ag_displacement_p01_k6x6x6_v19'],
 'shear_k8':['Ag_baseline_k8x8x8_v19','Ag_shear_m005_k8x8x8_v19','Ag_shear_p005_k8x8x8_v19'],
 'volume_k8':['Ag_baseline_k8x8x8_v19','Ag_volume_m02_k8x8x8_v19','Ag_volume_p02_k8x8x8_v19']}
refs={};audit={}
for group,labels in groups.items():
 for label in labels:
  if label in refs: continue
  f=folder/label
  s=json.loads((f/'summary.json').read_text());v=json.loads((f/'verification.json').read_text())
  assert s.get('label')==label and s.get('scf_converged') is True and v.get('status')=='PASS' and all(v['checks'].values())
  checks=[]
  for line in (f/'SHA256SUMS.txt').read_text().splitlines():
   if not line.strip():continue
   h,name=line.split(maxsplit=1);name=name.lstrip('*');p=(f/name).resolve();assert f.resolve() in p.parents and p.is_file();actual=hashlib.sha256(p.read_bytes()).hexdigest();assert h==actual;checks.append({'file':name,'sha256':actual})
  result=f/f'{label}_PW_PBE.extxyz';assert v['sha256'][result.name]==hashlib.sha256(result.read_bytes()).hexdigest()
  a=read(result,index=0);force=np.asarray(a.arrays['PW_PBE_forces'],float)
  assert len(a)==32 and a.get_chemical_formula()=='Ag32' and a.pbc.tolist()==[True]*3 and np.isfinite(force).all()
  refs[label]=(a,force,s);audit[label]={'SCF_converged':True,'verification':'PASS','sha256_files_checked':checks,'kpts':s['method']['kpts'],'native_eV':s['energy_eV_cell'],'free_eV':s['free_energy_eV_cell'],'source_sha256':s['source_sha256']}
models={}
for n,p in [('foundation',BASE/'periodic_interface_v4/mace_periodic_v12_interface_energy/foundation_models/mace-mp-0b3-medium.model'),('V18',BASE/'periodic_interface_v4/mace_periodic_v18_clean_core/training_clean_core_01/models/MACE_periodic_v18_clean_core.model')]:
 models[n]={'calc':MACECalculator(model_paths=str(p),device='cpu',default_dtype='float64'),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
scores={}
for group,labels in groups.items():
 base=labels[0];ab,fb,sb=refs[base];base_values={}
 for name,item in models.items():
  b=ab.copy();b.calc=item['calc'];base_values[name]=float(b.get_potential_energy())
 group_rows=[]
 for label in labels:
  atoms,dft_force,s=refs[label]
  assert s['method']['kpts']==sb['method']['kpts']
  item={'label':label,'group':group,'model_force_RMSE_eV_A':{},'model_relative_energy_meV_atom':{},'DFT_relative_native_meV_atom':(float(s['energy_eV_cell'])-float(sb['energy_eV_cell']))*1000/len(atoms),'DFT_relative_free_meV_atom':(float(s['free_energy_eV_cell'])-float(sb['free_energy_eV_cell']))*1000/len(atoms)}
  for name,model in models.items():
   a=atoms.copy();a.calc=model['calc'];e=float(a.get_potential_energy());pred=a.get_forces();delta=pred-dft_force
   item['model_force_RMSE_eV_A'][name]=float(np.sqrt(np.mean(np.sum(delta**2,axis=1))))
   item['model_relative_energy_meV_atom'][name]=(e-base_values[name])*1000/len(atoms)
  group_rows.append(item)
 scores[group]=group_rows
 if group=='displacement_k6':
  qminus=refs[labels[1]][0].positions; qplus=refs[labels[2]][0].positions
  moved=np.flatnonzero(np.linalg.norm(qplus-qminus,axis=1)>1e-8);assert len(moved)==1
  i=int(moved[0]);span=qplus[i,0]-qminus[i,0];assert abs(abs(span)-0.02)<1e-8
  responses={'moved_atom_index_zero_based':i,'moved_species':'Ag','signed_x_span_A':float(span),'DFT_restoring_slope_eV_A2':float(-(refs[labels[2]][1][i,0]-refs[labels[1]][1][i,0])/span)}
  for name,model in models.items(): responses[name+'_restoring_slope_eV_A2']=float(-(models[name]['calc'].get_forces(refs[labels[2]][0])[i,0]-models[name]['calc'].get_forces(refs[labels[1]][0])[i,0])/span)
  scores['symmetric_displacement_response']=responses
(OUT/'ag_response.json').write_text(json.dumps({'scope':'Eight same-composition Ag32 static DFT archives on matched k meshes; frozen CPU float64 foundation and V18 inference only; energy deltas are relative to the corresponding same-mesh baseline and native/free DFT conventions remain separate.','model_sha256':{k:v['sha256'] for k,v in models.items()},'archives':audit,'scores':scores,'limits':['Small static controls do not sample liquid Ag, melting, supercooling or solid-liquid interfaces.','Relative energy on these exact geometry families is diagnostic and not independent validation.']},indent=2)+'\n')
rows=[]
for group,values in scores.items():
 if isinstance(values,list):
  for value in values:
   row={'group':group,'label':value['label'],'DFT_relative_native_meV_atom':value['DFT_relative_native_meV_atom'],'DFT_relative_free_meV_atom':value['DFT_relative_free_meV_atom']}
   for m in models:row[f'{m}_relative_energy_meV_atom']=value['model_relative_energy_meV_atom'][m];row[f'{m}_force_RMSE_eV_A']=value['model_force_RMSE_eV_A'][m]
   rows.append(row)
with (OUT/'ag_response.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print('Ag same-mesh static response inference completed on',len(audit),'archives.')
