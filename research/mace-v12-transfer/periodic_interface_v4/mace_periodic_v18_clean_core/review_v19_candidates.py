import hashlib,json
from pathlib import Path
import numpy as np
from ase.io import read
from ase.geometry import find_mic
root=Path.cwd(); src=root/'research/mace-v12-transfer/coordination/reports/window-b/v19_residual_direction_geometry_design'; manifest=json.loads((src/'manifest.json').read_text())
parent_file=root/'research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v17_interface_energy/data/valid.extxyz';parents={a.info['config_type'].split('_v17_dev_validation')[0]:a for a in read(parent_file,':')}
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
rows=[]
for r in manifest['records']:
 if r['candidate_id'] not in manifest['recommended']:continue
 p=root/r['path'];assert sha(p)==r['sha256'];a=read(p);b=parents[r['interface']]
 assert r['screening_pass'] and not r['exact_hits'] and not r['contact_pair_flags'] and not r['parent_relative_contact_flags']
 assert np.array_equal(a.numbers,b.numbers) and np.array_equal(a.cell.array,b.cell.array) and np.array_equal(a.pbc,b.pbc)
 for key in ('lammps_id','central_pair'):assert np.array_equal(a.arrays[key],b.arrays[key])
 assert a.calc is None and not any('force' in k.lower() or 'energy' in k.lower() for k in list(a.info)+list(a.arrays))
 assert np.isfinite(a.positions).all()
 delta=find_mic(a.positions-b.positions,a.cell,a.pbc)[0];expected=np.zeros_like(delta)
 for moved in r['moved_atoms']:
  idx=np.flatnonzero(a.arrays['lammps_id']==moved['atom_id']);assert len(idx)==1 and a[idx[0]].symbol==moved['species'];expected[idx[0]]=moved['delta_xyz_A']
 assert np.allclose(delta,expected,rtol=0,atol=3e-8)
 assert abs(np.linalg.norm(delta,axis=1).max()-r['amplitude_A'])<3e-8
 distances=a.get_all_distances(mic=True);parent_distances=b.get_all_distances(mic=True)
 minima={}
 for i in range(len(a)):
  for j in range(i):
   pair='-'.join(sorted([a[i].symbol,a[j].symbol]));minima[pair]=min(minima.get(pair,float('inf')),float(distances[i,j]))
 for pair,m in r['species_minima_A'].items():assert abs(minima[pair]-m['distance_A'])<3e-8
 rows.append({'candidate_id':r['candidate_id'],'input_path':r['path'],'sha256':sha(p),'amplitude_A':r['amplitude_A'],'moved_atoms':r['moved_atoms'],'species_minima_A':minima,'identity_and_displacement_verified':True,'proposal_role':'future_training_acquisition_reused_development_family'})
assert len(rows)==6
result={'status':'GEOMETRY_REVIEW_PASS; DFT_NOT_STARTED','candidate_manifest_sha256':sha(src/'manifest.json'),'parent_file_sha256':sha(parent_file),'rows':rows,'coordinate_serialization_tolerance_A':3e-8,'limitations':['Intentional near-family correlation; not independent validation.','Parent-relative0.05A screening is not proof of physical representativeness.','Sealed near-overlap remains unknown; labels not inspected.'],'planned_owner_split':'A positive signs/B negative signs,3 each; draft only, actual DFT authorization requires reviewed execution pack.'}
output=root/'research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v18_clean_core/review_v19_candidates.json';output.write_text(json.dumps(result,indent=2)+'\n');print(result['status'],len(rows),'candidates')
