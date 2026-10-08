import json,hashlib,ast,subprocess,importlib.util
from pathlib import Path
import numpy as np
from ase.io import read
ROOT=Path(__file__).resolve().parents[6];OUT=Path(__file__).resolve().parent;BASE=ROOT/'research/mace-v12-transfer/periodic_interface_v4';entry=BASE/'mace_periodic_v18_clean_core';WB=OUT.parent;sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
m=json.loads((entry/'data/dataset_manifest.json').read_text());assert m['roles']=={'train':30,'development_validation':3,'test':0}
assert sha(entry/'data/train.extxyz')==sha(WB/'v17_clean_core_builder_draft/train_candidate.extxyz')==m['train_input_sha256'];assert sha(entry/'data/valid.extxyz')==sha(BASE/'mace_periodic_v17_interface_energy/data/valid.extxyz')==m['development_validation_input_sha256']
for p,h in m['source_files_sha256'].items():assert sha(ROOT/p)==h,(p,'hash mismatch')
a=read(entry/'data/train.extxyz',':');assert len(a)==30
for x in a:assert np.isfinite(x.info['REF_energy']) and x.arrays['REF_forces'].shape==(len(x),3) and np.isfinite(x.arrays['REF_forces']).all()
spec=importlib.util.spec_from_file_location('pf',BASE/'mace_periodic_v17_interface_energy/preflight_v17.py');pf=importlib.util.module_from_spec(spec);spec.loader.exec_module(pf);matrix=np.array([[x.get_chemical_symbols().count(s) for s in ['Ag','C','Si','Ti']] for x in a]);assert pf.exact_rank(matrix)==4
spec=importlib.util.spec_from_file_location('g',WB/'v17_split_geometry_design/final_set/generate_split_candidates.py');g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g);dev=g.extxyz_geometry_frames(entry/'data/valid.extxyz');assert len(dev)==3
for x in a:
 for y in dev:
  match=g.compare_geometry(x,y);assert not match or not(match['exact_duplicate'] or match['near_duplicate'])
assert len([r for r in m['rows'] if 'signed_design' in r])==8
for name in ['build_v18_dataset.py','finalize_v18.py','evaluate_v18_dev.py']:ast.parse((entry/name).read_text())
subprocess.run(['bash','-n',str(entry/'run_v18_training.sh')],check=True)
script=(entry/'run_v18_training.sh').read_text()
for token in ['--seed 45','--batch_size 1','--max_num_epochs 80','--lr 0.0001','--weight_decay 5e-7','--energy_weight 100','--forces_weight 1000','--save_all_checkpoints','--keep_checkpoints','epoch 79','mace-mp-0b3-medium']:assert token in script
assert '--test_file' not in script
notes=['Default threads4, allowed runtime override1–4: A should record actual threads to preserve controlled comparison.','Finalizer is safe in intended set-e runner flow; called standalone it sets training_exit_code0 from log evidence without independent launcher receipt.','Evaluator verifies current dataset against manifest but does not bind selection dataset_manifest/train/dev hashes back to current inputs; standalone selection/data drift guard could be strengthened.','Pair projection uses unwrapped Cartesian displacement, safe for current marked pairs but not a generic MIC implementation.']
result={'status':'ENTRY_AUDIT_PASS_WITH_LIMITATIONS','counts':m['roles'],'rank_exact':4,'train_candidate_byte_identity':True,'dev_byte_identity':True,'sources_hashes_pass':True,'signed_pairs':8,'frame0_caveat':m['frame0_kpoint_convergence'],'notes':notes,'script_hashes':{n:sha(entry/n) for n in ['build_v18_dataset.py','run_v18_training.sh','finalize_v18.py','evaluate_v18_dev.py']}}
(OUT/'audit.json').write_text(json.dumps(result,indent=2)+'\n');(OUT/'report.md').write_text('# V18 independent entry audit\n\nTrain30/dev3/test0, exact rank4, candidate/dev byte identity, all pinned source hashes, finite training labels and eight intended signed diagnostics pass. Train/dev geometry isolation passes at published exact/near thresholds. Frame0 source/order/cell preserved; original1x4x2 versus newGamma kpoint-convergence caveat remains. No withheld structure or label was opened.\n\nDefault recipe matches seed45/medium/lr1e-4/wd5e-7/batch1/80epochs/energy100/force1000, four threads, no test_file. Runner requires final epoch79 export; finalizer checks all80 epochs and completion markers, freezes selection before development scoring; evaluator exports per-atom predictions and fixed gates. Syntax/AST checks passed without training/inference.\n\nLimitations/recommendations:\n'+ '\n'.join('- '+x for x in notes)+'\n\nNo production edits or model calculations. This entry check does not prove numerical/physical convergence, force accuracy or independent generalization.\n');(OUT/'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='SHA256SUMS.txt'));print(result['status'])
