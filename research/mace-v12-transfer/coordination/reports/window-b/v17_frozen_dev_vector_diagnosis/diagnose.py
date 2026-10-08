import json,csv,hashlib
from pathlib import Path
import numpy as np
import torch
from ase.io import read
from mace.calculators import MACECalculator
ROOT=Path(__file__).resolve().parents[6];OUT=Path(__file__).resolve().parent;BASE=ROOT/'research/mace-v12-transfer/periodic_interface_v4';entry=BASE/'mace_periodic_v17_interface_energy';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
s=json.loads((entry/'training_provisional_03/selection_record.json').read_text());model=ROOT/s['model_path'];assert sha(model)==s['model_sha256']=='559b12c024fe8e2155e26ffc7ca48cf9fa1c2682ea37ff9c5367af14a5cc43d3';assert sha(entry/'data/valid.extxyz')==s['development_validation_input_sha256'];ev=json.loads((entry/'evaluation_dev_01/evaluation.json').read_text());assert sha(entry/'training_provisional_03/selection_record.json')==ev['selection_record_sha256']
p16=BASE/'pbe_interface_v17_main_targeted_acquisition/development_validation_evaluation_v16/per_atom.csv';old=list(csv.DictReader(p16.open()));torch.set_num_threads(1);torch.set_default_dtype(torch.float64);calc=MACECalculator(models=[torch.load(model,map_location='cpu',weights_only=False).double().eval()],device='cpu',default_dtype='float64')
atomrows=[];results=[]
for a,target in zip(read(entry/'data/valid.extxyz',':'),ev['per_structure']):
 label=a.info['config_type'];assert label==target['label'];refpath=BASE/'pbe_interface_v17_main_targeted_acquisition/development_validation_archives'/label/(label+'_PW_PBE.extxyz');assert sha(refpath)==target['reference_sha256']
 b=a.copy();b.calc=calc;energy=b.get_potential_energy();fp=b.get_forces();fr=a.arrays['REF_forces'];err=fp-fr;central=np.flatnonzero(a.arrays['central_pair']);i=next(k for k in central if a[k].symbol=='Ag');j=next(k for k in central if a[k].symbol!='Ag');axis=a.positions[j]-a.positions[i];axis/=np.linalg.norm(axis);dist=a.get_all_distances(mic=True);shell=np.minimum(dist[i],dist[j])<=4.5
 previous=[r for r in old if r['label']==label];assert len(previous)==len(a)
 e16=np.array([[float(r['error_'+d+'_eV_A']) for d in ['fx','fy','fz']] for r in previous]);assert all(int(r['atom_id'])==int(a.arrays['lammps_id'][k]) for k,r in enumerate(previous))
 longitudinal=err@axis;transverse=err-longitudinal[:,None]*axis
 rmse=lambda x:float(np.sqrt(np.mean(np.sum(x*x,axis=1))))
 metrics={'force_vector_RMSE_eV_A':rmse(err),'separating_force_abs_error_eV_A':float(abs(longitudinal[j]-longitudinal[i])),'energy_abs_error_meV_atom':float(abs(energy-a.info['REF_energy'])*1000/len(a))}
 for key,value in metrics.items():assert abs(value-target[key])<(1e-4 if key=='energy_abs_error_meV_atom' else 1e-6),(label,key,value,target[key])
 versions={}
 for name,e in [('V16',e16),('V17',err)]:
  proj=e@axis;perp=e-proj[:,None]*axis;versions[name]={'vector_RMSE':rmse(e),'longitudinal_RMS':float(np.sqrt(np.mean(proj**2))),'transverse_vector_RMS':rmse(perp),'pair_Ag_signed_error':float(proj[i]),'pair_X_signed_error':float(proj[j]),'pair_separation_signed_error':float(proj[j]-proj[i]),'shell_RMSE':rmse(e[shell]),'outside_RMSE':rmse(e[~shell]),'species':{sp:rmse(e[np.array(a.get_chemical_symbols())==sp]) for sp in set(a.get_chemical_symbols())}}
 norms=np.linalg.norm(err,axis=1);top=np.argsort(norms)[-5:][::-1];maxima=[]
 for k in top:
  neighbors=np.argsort(dist[k])[1:7];maxima.append({'id':int(a.arrays['lammps_id'][k]),'species':a[k].symbol,'error_vector':err[k].tolist(),'norm':float(norms[k]),'longitudinal':float(longitudinal[k]),'transverse':float(np.linalg.norm(transverse[k])),'neighbors':[{'id':int(a.arrays['lammps_id'][n]),'species':a[n].symbol,'distance_A':float(dist[k,n])} for n in neighbors]})
 results.append({'label':label,'reproduced_metrics':metrics,'versions':versions,'maxima_and_species_neighbors':maxima,'marked_pair_ids':[int(a.arrays['lammps_id'][i]),int(a.arrays['lammps_id'][j])]})
 for k in range(len(a)):
  row={'label':label,'index':k,'id':int(a.arrays['lammps_id'][k]),'species':a[k].symbol,'shell':bool(shell[k]),'marked':bool(a.arrays['central_pair'][k]),'longitudinal_error':float(longitudinal[k]),'transverse_error':float(np.linalg.norm(transverse[k]))}
  for prefix,vector in [('position',a.positions[k]),('reference',fr[k]),('prediction',fp[k]),('error',err[k])]:
   for d,v in zip('xyz',vector):row[prefix+'_'+d]=float(v)
  atomrows.append(row)
with (OUT/'per_atom.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(atomrows[0]));w.writeheader();w.writerows(atomrows)
(OUT/'diagnosis.json').write_text(json.dumps({'model_sha256':sha(model),'V16_per_atom_sha256':sha(p16),'results':results,'metric_tolerance':{'energy_meV_atom':1e-4,'force_eV_A':1e-6},'shell_radius_A':4.5,'axis':'direct marked Ag-to-X matching published evaluator'},indent=2)+'\n')
lines=['# Frozen V17 same-development vector diagnosis','','Published V17 energy/vector/separation metrics reproduced within1e-4 meV/atom energy and1e-6 eV/A force (CPU/runtime arithmetic tolerance). Only approved frozen epoch79 model and three scored development references were used; no unchanged V16 inference, fitting, DFT or withheld access. Signed separation error is projected X error minus projected Ag error, so improvements in bulk RMS need not improve this two-atom difference.','']
for r in results:
 lines.append('## '+r['label'])
 for name,v in r['versions'].items():lines.append(f"{name}: vector={v['vector_RMSE']:.6f}, shell={v['shell_RMSE']:.6f}, outside={v['outside_RMSE']:.6f}, transverse={v['transverse_vector_RMS']:.6f}; Ag signed={v['pair_Ag_signed_error']:.6f}, X signed={v['pair_X_signed_error']:.6f}, separation signed={v['pair_separation_signed_error']:.6f} eV/A.")
 lines.append('Largest V17 residual IDs: '+str([(m['id'],m['species'],round(m['norm'],5)) for m in r['maxima_and_species_neighbors']]))
lines+=['','## Acquisition recommendations for A review','First rank contact-shell neighbor residuals from this report and preserve pair contributions separately. For AgTi, target marked Ti plus nearest framework-C response and transverse Ti/Ag registry; propose opposite-sign small displacements along the observed residual directions, with geometry/contact checks and new-role isolation before any DFT approval. AgC prioritize large Ti transverse shell residual; AgSi prioritize marked Si/Ag differential projection and neighboring framework response. These directions are diagnostic hypotheses, not evidence of a sole cause or authorization to calculate.','The existing eight pairs probe particular parent environments; residual maxima on moved development registries can require additional neighbor coverage. A fixed-recipe clean30 comparison is useful but combines provenance filtering, coverage removal and mesh-changed frame0; separate these confounds. Keep force/separation gates unchanged and sealed tests unavailable for tuning.']
(OUT/'report.md').write_text('\n'.join(lines)+'\n');(OUT/'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='SHA256SUMS.txt'));print('frozen inference and metrics PASS')
