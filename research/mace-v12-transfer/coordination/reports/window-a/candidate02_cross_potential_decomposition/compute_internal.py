import csv,hashlib,json,tempfile
from pathlib import Path
import numpy as np
from ase.io import read,write
from ase.geometry import find_mic
from lammps import lammps
ROOT=Path(__file__).resolve().parent;REPO=ROOT.parents[5]
BUNDLE=REPO/'research/mace-v12-transfer/exports/cross_potential_dft_reference_bundle_20261009'
POT=REPO/'research/mace-v12-transfer/coordination/reports/window-b/current_heating_potential_handoff_20261009'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for line in (POT/'SHA256SUMS').read_text().splitlines():
 digest,name=line.split('  ',1);assert sha(POT/name)==digest
manifest=json.loads((BUNDLE/'bundle_manifest.json').read_text());rows=[];pred=[]
for e in manifest['entries']:
 a=read(BUNDLE/e['input']);ref=read(BUNDLE/e['result_extxyz'])
 assert sha(BUNDLE/e['input'])==e['input_sha256']
 assert np.allclose(a.positions,ref.positions,rtol=0,atol=1e-12) and np.allclose(a.cell,ref.cell,rtol=0,atol=1e-12)
 assert a.get_chemical_symbols()==ref.get_chemical_symbols() and np.array_equal(a.pbc,ref.pbc)
 assert np.array_equal(a.arrays['lammps_id'],ref.arrays['lammps_id'])
 with tempfile.TemporaryDirectory(prefix='candidate02-static-') as tmp:
  data=Path(tmp)/'input.data';write(data,a,format='lammps-data',atom_style='atomic',specorder=['Ti','Si','C','Ag'],masses=True)
  l=lammps(cmdargs=['-log','none','-screen','none'])
  try:
   for c in ['units metal','atom_style atomic','boundary '+' '.join('p' if v else 'f' for v in a.pbc),'read_data '+str(data)]:l.command(c)
   for c in (POT/'potential.interface_calibration_v0.inc').read_text().splitlines():
    if c.strip() and not c.lstrip().startswith('#'):
     if ' morse ' in c and c.startswith('pair_coeff'):
      fields=c.split();fields[4]='0.0';c=' '.join(fields)
     l.command(c.replace('potentials/',str(POT/'potentials')+'/'))
   l.command('run 0');energy=float(l.get_thermo('pe'));version=l.version()
   ids=l.numpy.extract_atom('id').copy()[:len(a)];order=np.argsort(ids)
   assert np.array_equal(ids[order],np.arange(1,len(a)+1))
   xyz=l.numpy.extract_atom('x').copy()[:len(a)][order];f=l.numpy.extract_atom('f').copy()[:len(a)][order]
   delta,lengths=find_mic(xyz-a.positions,a.cell,a.pbc)
   assert np.max(lengths)<1e-10 and np.isfinite(f).all() and np.isfinite(energy), e["label"]
  finally:l.close()
 df=f-ref.arrays['PW_PBE_forces'];pair=np.flatnonzero(a.arrays['central_pair']);assert len(pair)==2
 i,j=pair
 if a[i].symbol!='Ag':i,j=j,i
 v=a.get_distance(int(i),int(j),mic=True,vector=True);unit=v/np.linalg.norm(v)
 modelpair=float(np.dot(f[j]-f[i],unit));refpair=float(np.dot(ref.arrays['PW_PBE_forces'][j]-ref.arrays['PW_PBE_forces'][i],unit))
 rms=float(np.sqrt(np.mean(np.sum(df*df,axis=1))));sep=abs(modelpair-refpair)
 rows.append(dict(label=e['label'],role=e['role'],contact='Ag-'+a[j].symbol,distance_A=float(np.linalg.norm(v)),force_vector_rmse_eV_A=rms,max_atom_error_eV_A=float(np.linalg.norm(df,axis=1).max()),reference_pair_eV_A=refpair,model_pair_eV_A=modelpair,separation_error_eV_A=sep,vector_pass=rms<=.05,separation_pass=sep<=.10,model_energy_eV=energy,native_DFT_eV=e['native_energy_eV'],free_DFT_eV=e['free_energy_eV'],raw_energy_gate='NOT_SCORED_DIFFERENT_ZEROS'))
 pred.append({'label':e['label'],'input_sha256':e['input_sha256'],'original_atom_ids':a.arrays['lammps_id'].tolist(),'symbols':a.get_chemical_symbols(),'energy_eV':energy,'forces_eV_A':f.tolist(),'per_species_vector_rmse_eV_A':{s:float(np.sqrt(np.mean(np.sum(df[np.array(a.get_chemical_symbols())==s]**2,axis=1)))) for s in sorted(set(a.get_chemical_symbols()))}})
with (ROOT/'metrics.csv').open('w') as out:
 w=csv.DictWriter(out,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
(ROOT/'predictions.json').write_text(json.dumps(pred,indent=2)+'\n')
# Same-parent signed diagnostics only: relative differences cancel constant energy zeros.
responses=[]
for contact in ['Ag-C','Ag-Si','Ag-Ti']:
 rr=[r for r in rows if r['contact']==contact and r['role']=='reused_training_family_local_reference'];assert len(rr)==2
 minus=next(r for r in rr if '_m_' in r['label']);plus=next(r for r in rr if '_p_' in r['label'])
 responses.append({'contact':contact,'labels':[minus['label'],plus['label']],'model_delta_eV':plus['model_energy_eV']-minus['model_energy_eV'],'DFT_native_delta_eV':plus['native_DFT_eV']-minus['native_DFT_eV'],'DFT_free_delta_eV':plus['free_DFT_eV']-minus['free_DFT_eV'],'note':'Local signed residual perturbation, not a general normal-distance scan; native/free kept distinct.'})
report={'status':'INTERNAL_ONLY_COUNTERFACTUAL_DIAGNOSIS','lammps_version':version,'potential_hashes':{str(p.relative_to(POT)):sha(p) for p in sorted(POT.rglob('*')) if p.is_file()},'cases':len(rows),'all_force_gates_pass':all(r['vector_pass'] and r['separation_pass'] for r in rows),'local_signed_responses':responses,'constraints':['Not confirmed Stage69 configuration','Not independent validation: reused training parents and same-geometry numerical controls','No absolute-energy gate: energy zeros unmatched and no separated-phase reference','Counterfactual Morse D=0 only in memory; source files unchanged; no adopted parameters/dynamics/heating/sealed access']}
(ROOT/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'cases':len(rows),'force_gates_pass':report['all_force_gates_pass'],'diagnostics':[r for r in rows if r['role']=='reused_training_family_local_reference'],'responses':responses},indent=2))
