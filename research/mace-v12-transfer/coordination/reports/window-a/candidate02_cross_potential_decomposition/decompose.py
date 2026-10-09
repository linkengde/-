import json
from pathlib import Path
import numpy as np
from ase.io import read
from ase.neighborlist import neighbor_list
ROOT=Path(__file__).resolve().parent;REPO=ROOT.parents[5];bundle=REPO/'research/mace-v12-transfer/exports/cross_potential_dft_reference_bundle_20261009'
full=json.loads((ROOT.parent/'candidate02_cross_potential_dft_diagnosis/predictions.json').read_text());internal=json.loads((ROOT/'predictions.json').read_text());manifest=json.loads((bundle/'bundle_manifest.json').read_text());rows=[]
params={'Ti':(.94814194,2.,2.33708864),'Si':(.93944948,2.,2.40407837),'C':(1.35077230,2.,2.15255924)}
for e,p,q in zip(manifest['entries'],full,internal):
 assert e['label']==p['label']==q['label'];a=read(bundle/e['input']);ref=read(bundle/e['result_extxyz']).arrays['PW_PBE_forces'];ff=np.array(p['forces_eV_A']);fi=np.array(q['forces_eV_A']);fc=np.zeros_like(fi);ec=0.;symbols=np.array(a.get_chemical_symbols());i,j,d,v=neighbor_list('ijdD',a,8.)
 for ii,jj,rr,vec in zip(i,j,d,v):
  if (symbols[ii]=='Ag')==(symbols[jj]=='Ag'):continue
  x=symbols[jj] if symbols[ii]=='Ag' else symbols[ii];D,alpha,r0=params[x];z=np.exp(-alpha*(rr-r0));ec+=.5*D*(z*z-2*z);fc[ii]+=2*alpha*D*(z-z*z)*vec/rr
 assert np.allclose(ff-fi,fc,rtol=1e-10,atol=1e-9),e['label']
 assert abs(p['energy_eV']-q['energy_eV']-ec)<1e-8,e['label']
 # Necessary lower bound: arbitrary central cross-pair amplitudes cannot change
 # any force component orthogonal to all Ag-X neighbour directions at an atom.
 residual=[]
 for atom in range(len(a)):
  directions=[vec/rr for ii,jj,rr,vec in zip(i,j,d,v) if ii==atom and ((symbols[ii]=='Ag')!=(symbols[jj]=='Ag'))]
  target=ref[atom]-fi[atom]
  if directions:
   A=np.array(directions).T;projection=A@np.linalg.lstsq(A,target,rcond=1e-12)[0];residual.append(target-projection)
  else:residual.append(target)
 norm=lambda f:float(np.sqrt(np.mean(np.sum(f*f,axis=1))))
 rows.append({'label':e['label'],'role':e['role'],'full_vector_rmse_eV_A':norm(ff-ref),'internal_only_vector_rmse_eV_A':norm(fi-ref),'cross_force_vector_rms_eV_A':norm(fc),'cross_energy_eV':ec,'maximum_cross_reconstruction_error_eV_A':float(np.abs(ff-fi-fc).max()),'central_cross_only_necessary_lower_bound_eV_A':norm(np.array(residual)),'by_species':{s:{'full_rmse':norm((ff-ref)[symbols==s]),'internal_only_rmse':norm((fi-ref)[symbols==s]),'cross_rms':norm(fc[symbols==s])} for s in sorted(set(symbols))}})
report={'status':'STATIC_DECOMPOSITION_RECONSTRUCTION_PASS','definition':'Internal EAM/Tersoff unchanged, Morse D=0 for internal calculation. Independent analytical unshifted cutoff8 Morse reconstructs full-minus-internal energy/forces including periodic neighbours.','cases':rows,'limits':['Internal-only error is not pure-phase validation: reference still includes interface interactions','Negative cross potential energy alone does not prove dynamic heat release','Necessary geometric lower bound permits arbitrary per-neighbour central amplitudes: passing it does not prove any shared Morse parameters can fit','No fitting or parameter adoption; original files unchanged; not Stage69-confirmed or independent validation']}
(ROOT/'decomposition.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'diagnostics':[r for r in rows if r['role']=='reused_training_family_local_reference']},indent=2))
