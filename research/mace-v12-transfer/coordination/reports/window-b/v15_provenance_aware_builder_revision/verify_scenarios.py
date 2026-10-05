from pathlib import Path
import json,hashlib,subprocess
import numpy as np
from ase.io import read
root=Path('/workspace/-'); rel='research/mace-v12-transfer/coordination/reports/window-b/v15_provenance_aware_builder_revision'; p=root/rel
old=root/'research/mace-v12-transfer/coordination/reports/window-b/v15_builder_draft/generated_data'
scenarios={}
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def equal(a,b):
 assert np.array_equal(a.numbers,b.numbers) and np.array_equal(a.pbc,b.pbc)
 assert np.allclose(a.positions,b.positions,atol=1e-8,rtol=0) and np.allclose(a.cell,b.cell,atol=1e-8,rtol=0)
 assert set(a.arrays)==set(b.arrays)
 for k in a.arrays:assert np.allclose(a.arrays[k],b.arrays[k],atol=1e-8,rtol=0),k
 assert set(a.info)==set(b.info)
 for k in a.info:assert np.array_equal(a.info[k],b.info[k]),k
 assert np.isfinite(b.arrays['REF_forces']).all() and np.isfinite(b.info['REF_energy'])
reference=read(old/'train.extxyz',index=':')
for name,count in [('baseline_30',30),('exclude_unresolved_29',29)]:
 d=p/name; m=json.loads((d/'dataset_manifest.json').read_text()); v=json.loads((d/'verification_report.json').read_text())
 frames=read(d/'train.extxyz',index=':'); expected=reference if count==30 else reference[:9]+reference[10:]
 assert len(frames)==len(expected)==count
 for a,b in zip(expected,frames):equal(a,b)
 for s in ['valid','test']:assert (d/(s+'.extxyz')).read_bytes()==(old/(s+'.extxyz')).read_bytes()
 assert m['energy_composition_matrix']['rank']==4 and all(not hits for hits in v['frozen_V15_geometry_overlap_by_split'].values())
 assert m['evidence_tier_counts_train']=={'partial_recovery':14,**({'unresolved':1} if count==30 else {}),'inherited_declaration_not_reverified':12,'verified_residual_archive':3}
 for path,h in m['sources_sha256'].items():assert sha(root/path)==h
 assert len([e for e in m['excluded_frame_ledger'] if e['config_type']=='AgTi_2p3871A_periodic_PW_PBE_force_only'])==(count==29)
 for r in m['per_frame_provenance']:assert 'source_method_as_stored' in r and 'dft_method_as_stored' in r and 'original_run_evidence' in r
 scenarios[name]={'split_sizes':m['split_sizes'],'composition_rank':4,'evidence_tier_counts_train':m['evidence_tier_counts_train'],'all_retained_geometry_labels_metadata_equal_prior_draft':True,'valid_test_byte_identical':True,'finite_labels':True,'source_checksums_match':True,'blind_input_overlap_by_split':v['frozen_V15_geometry_overlap_by_split'],'historical_overlap_counts':{s:len(h) for s,h in v['historical_V14_training_split_overlaps'].items()},'final_overlap_counts':{s:len(h) for s,h in v['draft_V15_training_split_overlaps'].items()},'excluded_frame_ledger':m['excluded_frame_ledger'],'train_sha256':sha(d/'train.extxyz')}
script=p/'build_v15_draft.py'; python='/workspace/.venvs/gpaw-mpi/bin/python'
# User-requested validation: existing runs and production paths must be refused.
refusals=[]
for path in [p/'baseline_30',root/'research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v15_interface_energy/data']:
 result=subprocess.run([python,'-B',str(script),'--repo-root',str(root),'--output-dir',str(path)],cwd='/tmp',text=True,capture_output=True)
 assert result.returncode!=0
 refusals.append({'output_path':str(path.relative_to(root)),'exit_code':result.returncode,'diagnostic':result.stderr.strip()})
result={'scenarios':scenarios,'refusal_checks':refusals,'no_blind_labels_read':True,'no_calculations':True,'comparison_tolerances':{'position_cell_A':1e-8,'force_eV_A':1e-8,'energy_metadata':'exact stored equality'},'recommendation':'29-frame scenario for limited numerical screening only; A decides integration. 14 partial and 12 explicitly declared inherited frames still lack original-run revalidation.'}
(p/'comparison.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({name:{k:v for k,v in x.items() if k not in ['excluded_frame_ledger']} for name,x in scenarios.items()},indent=2))
