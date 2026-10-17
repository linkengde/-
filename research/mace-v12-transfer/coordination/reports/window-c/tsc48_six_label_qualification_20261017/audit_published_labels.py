"""Offline read-only audit. No calculator execution. Run with ASE/numpy installed.
Arguments: exported E result root, previous C manifest, output directory.
"""
import sys,json,hashlib,re,csv
from pathlib import Path
import numpy as np
from ase.io import read
root,prior,out=map(Path,sys.argv[1:]);out.mkdir(parents=True,exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
prior=json.loads(prior.read_text());old={r['structure_id']:r for r in prior['six_input_records']}
manpath=root/'tscbop_run0_sixpoint/analysis/nonholdout_backbone_next_step_20261017/P0_SIX_POINT_PUBLICATION_MANIFEST.json'
man=json.loads(manpath.read_text()); records=[];checks=0;sources=[]
for sums in root.rglob('SHA256SUMS.txt'):
 for line in sums.read_text().splitlines():
  if not line.strip():continue
  h,n=line.split(maxsplit=1);p=(sums.parent/n.lstrip('*')).resolve()
  assert p.is_relative_to(root.resolve()),str(p)
  assert p.is_file() and sha(p)==h,str(p)
  checks+=1
for rec in man['records']:
 d=root/Path(rec['path']).relative_to('research/mace-v12-transfer/coordination/reports/window-e/dft_supplement_results_20261016')
 r=json.loads((d/'RESULT.json').read_text());status=json.loads((d/'RUN_STATUS.json').read_text());command=json.loads((d/'COMMAND.json').read_text());log=(d/'gpaw.log').read_text()
 for fn,key in [('input.extxyz','input_sha256'),('labelled.extxyz','labelled_extxyz_sha256'),('RESULT.json','result_json_sha256'),('RUN_STATUS.json','run_status_sha256'),('SHA256SUMS.txt','run_sha256_manifest_sha256')]:assert sha(d/fn)==rec[key]
 assert sha(d/'input.extxyz')==old[r['task_id']]['input_sha256']==r['input_sha256']
 assert sha(d/'labelled.extxyz')==r['output_extxyz_sha256']
 a=read(d/'input.extxyz');b=read(d/'labelled.extxyz');assert len(a)==len(b)==48
 assert a.get_chemical_symbols()==b.get_chemical_symbols()
 assert {s:a.get_chemical_symbols().count(s) for s in ['Ti','Si','C']}=={'Ti':24,'Si':8,'C':16}
 assert np.array_equal(a.positions,b.positions) and np.array_equal(a.cell.array,b.cell.array) and np.array_equal(a.pbc,b.pbc) and a.pbc.all()
 for name in ['local_new_id','spacegroup_kinds']:assert np.array_equal(a.arrays[name],b.arrays[name])
 f=np.array(r['forces_eV_per_A']);s=np.array(r['stress_eV_per_A3']);assert f.shape==(48,3) and s.shape==(3,3)
 assert np.isfinite(f).all() and np.isfinite(s).all();assert np.allclose(s,s.T,atol=1e-12)
 assert np.max(np.abs(f-b.calc.results['forces']))<=5.1e-9
 assert abs(b.calc.results['energy']-r['energy_eV'])<1e-10 and abs(b.calc.results['free_energy']-r['free_energy_eV'])<1e-10
 assert np.allclose(np.array(r['stress_GPa']),s*160.21766208,atol=2e-6,rtol=0)
 assert np.allclose(b.calc.results['stress'],[s[0,0],s[1,1],s[2,2],s[1,2],s[0,2],s[0,1]],atol=1e-12,rtol=0)
 assert r['scf_converged'] and r['status']=='SUCCESS' and status.get('scf_converged', True)
 assert f"Converged in {r['scf_iterations']} steps" in log
 for label,key in [('Extrapolated','energy_eV'),('Free energy','free_energy_eV')]:
  value=float(re.search(r'  '+label+r':\s+(-?[\d.]+)',log).group(1));assert abs(value-r[key])<=5.1e-7
 logforces=np.array([[float(x) for x in row] for row in re.findall(r'\|\s*(-?[\d.]+)\s*\|\s*(-?[\d.]+)\s*\|\s*(-?[\d.]+)\s*\|\s*# (?:Ti|Si|C)\s+\d+',log)])
 assert logforces.shape==(48,3) and np.max(np.abs(logforces-f))<5.1e-6
 stresspart=log.split('Stress tensor  # eV/Ang^3')[-1]
 logs=np.array([[float(x) for x in row] for row in re.findall(r'\|\s*(-?[\d.]+)\s*\|\s*(-?[\d.]+)\s*\|\s*(-?[\d.]+)\s*\|',stresspart)[:3]])
 assert logs.shape==(3,3) and np.max(np.abs(logs-s))<=5.1e-7
 for token in ["PW(ecut=500.0)","size=array([6, 6, 2])","width=0.1","xc=XC(name='PBE')","'energy': 1e-05","'density': 1e-05","'eigenstates': 1e-05","mixer=MSR1"]:assert token in log,token
 assert r['kpoints']==[6,6,2] and r['parameter_lock_sha256']==man['parameter_lock_sha256']==command['parameter_lock_sha256']
 for p in d.iterdir():
  if p.is_file():sources.append({'path':str(p.relative_to(root)),'sha256':sha(p)})
 records.append({'structure_id':r['task_id'],'source_family':'bulk_parent_48','input_sha256':sha(d/'input.extxyz'),'result_sha256':sha(d/'RESULT.json'),'label_sha256':sha(d/'labelled.extxyz'),'scf_iterations':r['scf_iterations'],'energy_eV':r['energy_eV'],'free_energy_eV':r['free_energy_eV'],'max_force_eV_A':float(np.linalg.norm(f,axis=1).max()),'launcher_exit_code':rec['launcher_exit_code'],'label_qualification':'PASS_WITH_LIMITATIONS_FOR_DIAGNOSTICS','data_role':'DIAGNOSTIC_ONLY','training_use':'HOLD_PENDING_NUMERICAL_AND_PROVENANCE_REVIEW_AND_EXPLICIT_USE_APPROVAL','stale_extxyz_label_status':b.info.get('dft_label_status'),'PAW_log_ids':re.findall(r'    id:\s+(\w+)',log),'parameter_lock_sha256':r['parameter_lock_sha256'],'runner_sha256':command['runner_sha256']})
assert len(records)==6
result={'schema_version':'2.0','e_commit':'f9f2e88e34b4c1e18e6a42ed16fd7c949f7c22ef','family':'bulk_parent_48','checks_passed':{'sha_entries':checks,'records':6,'geometry_order_ids':True,'result_vs_xyz_and_raw_log':True},'TRAIN':[],'DEV':[],'TEST':[],'role_change_authorized':False,'fit_authorized':False,'new_calculation_authorized':False,'records':records,'limitations':['Five missing supervisor exit codes remain unknown; scientific completion independently supported by raw SCF/forces/stress outputs.','Actual 6x6x2 mesh and MSR1 differ from original 10x10x4/Pulay proposal; no new geometry k-point/stress convergence demonstrated.','Referenced execution lock and runner source not present at their stated paths in reviewed commit; PAW archive SHA runtime binding pending.','labelled.extxyz carries stale PENDING_NO_DFT_RUN metadata; original bytes retained and qualified result interpreted through this manifest.']}
(out/'FAMILY_USAGE_MANIFEST_v2.json').write_text(json.dumps(result,indent=2)+'\n')
(out/'SOURCE_SHA256.json').write_text(json.dumps(sources,indent=2)+'\n')
with (out/'SIX_POINT_LABEL_REVIEW.csv').open('w') as f:
 fields=['structure_id','scf_iterations','energy_eV','free_energy_eV','max_force_eV_A','launcher_exit_code','label_qualification','data_role'];w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore',lineterminator='\n');w.writeheader();w.writerows(records)
print(json.dumps(result['checks_passed']));print([(r['structure_id'],r['scf_iterations'],r['PAW_log_ids']) for r in records])
