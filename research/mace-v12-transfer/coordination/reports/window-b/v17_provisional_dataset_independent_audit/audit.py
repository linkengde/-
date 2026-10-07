import json,hashlib,importlib.util
from pathlib import Path
from collections import Counter
import numpy as np
from ase.io import read
ROOT=Path(__file__).resolve().parents[6];OUT=Path(__file__).resolve().parent;WB=OUT.parent;BASE=ROOT/'research/mace-v12-transfer/periodic_interface_v4';data=BASE/'mace_periodic_v17_interface_energy/data'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
spec=importlib.util.spec_from_file_location('geom',WB/'v17_split_geometry_design/final_set/generate_split_candidates.py');g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)
spec=importlib.util.spec_from_file_location('rank',BASE/'mace_periodic_v17_interface_energy/preflight_v17.py');pf=importlib.util.module_from_spec(spec);spec.loader.exec_module(pf)
m=json.loads((data/'dataset_manifest.json').read_text());assert sha(data/'train.extxyz')==m['train_input_sha256'];train=read(data/'train.extxyz',':');old=read(BASE/'mace_periodic_v16_interface_energy/data/train.extxyz',':');assert len(train)==43 and m['roles']=={'train':43,'development_validation':3,'test':0}
ledger=[]
for i,a in enumerate(train):
 assert np.isfinite(a.info['REF_energy']) and a.arrays['REF_forces'].shape==(len(a),3) and np.isfinite(a.arrays['REF_forces']).all()
 if i<35:
  b=old[i];assert np.array_equal(a.numbers,b.numbers) and np.array_equal(a.positions,b.positions) and np.array_equal(a.cell,b.cell) and np.array_equal(a.pbc,b.pbc)
  assert a.info['REF_energy']==b.info['REF_energy'] and np.array_equal(a.arrays['REF_forces'],b.arrays['REF_forces'])
 else:
  label=m['verified_V17_pair_labels'][i-35];directory='pbe_interface_v17_main_targeted_acquisition' if '_p_' in label else 'pbe_interface_v17_parallel_targeted_acquisition';p=BASE/directory/'calculations'/label/(label+'_PW_PBE.extxyz');b=read(p)
  assert np.array_equal(a.numbers,b.numbers) and np.array_equal(a.positions,b.positions) and np.array_equal(a.cell,b.cell) and np.array_equal(a.pbc,b.pbc)
  assert a.info['REF_energy']==b.info['PW_PBE_energy_eV'] and np.array_equal(a.arrays['REF_forces'],b.arrays['PW_PBE_forces'])
 ledger.append({'frame':i,'config':a.info.get('config_type'),'formula':a.get_chemical_formula(),'role':a.info.get('v17_role'),'provenance':a.info.get('v17_provenance_status'),'source_copy':'PASS','method_declaration':a.info.get('dft_method',a.info.get('source_method'))})
N=np.array([[a.get_chemical_symbols().count(s) for s in m['composition_matrix_elements']] for a in train]);assert pf.exact_rank(N)==4 and N.tolist()==m['composition_matrix']
assert dict(Counter(x['provenance'] for x in ledger))==m['train_provenance_status_counts']
hashes=[]
for rel,h in m['source_files_sha256'].items():
 if any(x in rel for x in ['development_validation_archives','development_validation_evaluation','withheld']):
  hashes.append({'path':rel,'status':'not_opened_restricted_source'});continue
 p=ROOT/rel;actual=sha(p) if p.is_file() else None;hashes.append({'path':rel,'expected':h,'actual':actual,'status':'PASS' if actual==h else 'STALE_OR_MISSING'})
# Development geometry only; withheld exclusion via frozen role canonical hashes.
dev=g.extxyz_geometry_frames(data/'valid.extxyz');assert len(dev)==3
roles=json.loads((WB/'v17_split_geometry_design/final_set/role_manifest.json').read_text());withheld_hashes={r['geometry_sha256'] for r in roles['candidates'] if r['proposed_role']=='withheld_test'}
overlap=[]
for i,a in enumerate(train):
 assert g.geometry_sha(a) not in withheld_hashes
 for j,b in enumerate(dev):
  match=g.compare_geometry(a,b)
  if match and (match['exact_duplicate'] or match['near_duplicate']):overlap.append({'train':i,'dev':j,'match':match})
assert not overlap
result={'status':'PROVISIONAL_TRAIN_AUDIT_PASS_WITH_SOURCE_SNAPSHOT_NOTES','train_count':43,'composition_rank_exact':4,'provenance_counts':dict(Counter(x['provenance'] for x in ledger)),'frame0_new_label_integrated':False,'train_rows':ledger,'source_hash_checks':hashes,'development_geometry_overlap':overlap,'withheld_geometry_hash_overlap':False,'restricted_outputs_opened':False}
(OUT/'audit.json').write_text(json.dumps(result,indent=2)+'\n');stale=[x['path'] for x in hashes if x['status']=='STALE_OR_MISSING']
(OUT/'report.md').write_text('# Independent provisional V17 dataset audit\n\n43 training frames verified:35 inherited labels exactly preserved and8 signed diagnostic labels exactly copied from source output. Exact composition rank4;29 directly supported and14 partial-source rows. Newly computed frame0 replacement is not integrated: source frame0 remains unchanged. Finite training energies/forces, ordered species, coordinates, cell/PBC and manifest provenance counts pass. Three development geometries were parsed without label values; no exact/near train overlap under published thresholds. Frozen withheld geometry hashes do not match training; withheld inputs and outputs were not opened.\n\nSource snapshot changes: '+json.dumps(stale)+'\n\nAny stale partial aggregate manifest reflects subsequent B aggregate reconciliation rather than a changed label. Per-source details in audit.json retain expected/current hashes; other unexplained changes require A review. Do not silently rewrite the production manifest. Original inherited method/convention uncertainties remain: this pool is PROVISIONAL, not fully PW-PBE certified. Training-source checks and full rank do not certify model force quality. Builder/preflight were inspected; preflight was not run because it reads development reference values. No production dataset changes or model/DFT calculations during audit.\n')
(OUT/'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='SHA256SUMS.txt'))
print('train43 rank4; source snapshot discrepancies:',stale)
