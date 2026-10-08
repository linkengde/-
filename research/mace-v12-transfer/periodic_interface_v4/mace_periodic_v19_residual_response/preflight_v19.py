"""Read-only training gate; requires approved independent development and B audit."""
import hashlib,importlib.util,json,subprocess,sys
from pathlib import Path
import numpy as np
from ase.io import read
ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[3]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
subprocess.run([sys.executable,str(ROOT/'build_training_candidate.py'),'--check'],check=True)
data=ROOT/'data'
if not (data/'dataset_manifest.json').exists():raise SystemExit('BLOCKED: independent V19 development labels and official dataset not integrated; no training')
m=json.loads((data/'dataset_manifest.json').read_text())
assert m['status']=='A_APPROVED_V19_DATASET' and m['roles']=={'train':36,'development_validation':3,'test':0}
assert (data/'train.extxyz').read_bytes()==(ROOT/'training_candidate/train.extxyz').read_bytes()
assert sha(data/'train.extxyz')==m['train_sha256'] and sha(data/'valid.extxyz')==m['development_sha256']
assert m['withheld_labels_opened'] is False and not (data/'test.extxyz').exists()
for p,h in m['source_files_sha256'].items():
 f=(REPO/p).resolve();assert f.is_relative_to(REPO) and sha(f)==h
train=read(data/'train.extxyz',':');dev=read(data/'valid.extxyz',':');assert len(train)==36 and len(dev)==3
for a in train+dev:assert np.isfinite(a.info['REF_energy']) and a.arrays['REF_forces'].shape==(len(a),3) and np.isfinite(a.arrays['REF_forces']).all()
assert len(m['development_rows'])==3 and len(set(r['interface'] for r in m['development_rows']))==3
for r in m['development_rows']:assert r['role']=='independent_development' and r['A_lineage_review']=='APPROVED_NOT_REUSED_V17_DEVELOPMENT_FAMILY'
spec=importlib.util.spec_from_file_location('geometry',REPO/'research/mace-v12-transfer/coordination/reports/window-b/v17_split_geometry_design/final_set/generate_split_candidates.py');geom=importlib.util.module_from_spec(spec);spec.loader.exec_module(geom)
parents=read(ROOT.parent/'mace_periodic_v17_interface_energy/data/valid.extxyz',':')
for a in dev:
 for b in train+parents:
  result=geom.compare_geometry(a,b);assert not result or not(result['exact_duplicate'] or result['near_duplicate']),'Development/acquisition-family geometric overlap'
approval=ROOT/'entry_review_approval.json'
if not approval.exists():raise SystemExit('BLOCKED: B official entry audit and A approval not recorded; no training')
review=json.loads(approval.read_text());assert review['status']=='A_APPROVED_AFTER_B_AUDIT' and review['dataset_manifest_sha256']==sha(data/'dataset_manifest.json')
for p,h in review['review_files_sha256'].items():assert sha(REPO/p)==h
assert m['foundation_model_sha256']==sha(ROOT.parent/'mace_periodic_v12_interface_energy/foundation_models/mace-mp-0b3-medium.model')
print('V19_PREFLIGHT_PASS: train36/dev3/test0, roles/source/isolation/audit approvals verified')
