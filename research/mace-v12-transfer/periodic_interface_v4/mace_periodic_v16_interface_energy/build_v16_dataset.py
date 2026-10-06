#!/usr/bin/env python3
"""V16 controlled six-label augmentation; fresh holdout references never read."""
import hashlib, importlib.util, json, shutil
from pathlib import Path
import numpy as np
from ase.io import read, write
ROOT=Path(__file__).resolve().parent
PROJECT=ROOT.parent
REPO=ROOT.parents[3]
DATA=ROOT/'data'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def require(ok,msg):
    if not ok:raise ValueError(msg)
require(not DATA.exists() or not any(DATA.iterdir()),'Existing V16 data must be preserved')
parent=PROJECT/'mace_periodic_v15_interface_energy/data'
m=json.loads((parent/'dataset_manifest.json').read_text())
for split in ('train','valid','test'):require(sha(parent/(split+'.extxyz'))==m['output_sha256'][split],'Parent split hash mismatch')
train=read(parent/'train.extxyz',index=':')
bp=REPO/'research/mace-v12-transfer/coordination/reports/window-b/v15_provenance_aware_builder_revision/build_v15_draft.py'
require(sha(bp)==m['reviewed_builder_sha256'],'Reviewed geometry implementation changed')
spec=importlib.util.spec_from_file_location('reviewed_geometry',bp);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
sources={str(parent.relative_to(REPO)/'dataset_manifest.json'):sha(parent/'dataset_manifest.json')}
additions=[]
for name in ('pbe_interface_v16_main_acquisition','pbe_interface_v16_parallel_acquisition'):
    folder=PROJECT/name;manifest=json.loads((folder/'input_manifest.json').read_text())
    sources[str((folder/'input_manifest.json').relative_to(REPO))]=sha(folder/'input_manifest.json')
    for rec in manifest['records']:
        result=folder/'calculations'/rec['label'];verification=json.loads((result/'verification.json').read_text());summary=json.loads((result/'summary.json').read_text());progress=json.loads((result/'progress.json').read_text())
        require(verification['status']=='PASS' and all(verification['checks'].values()) and summary['scf_converged'] is True and progress['status']=='complete','Incomplete acquisition')
        for member,h in verification['sha256'].items():
            require(Path(member).name==member and sha(result/member)==h,'Archive hash mismatch')
            sources[str((result/member).relative_to(REPO))]=h
        inp=folder/rec['input'];require(sha(inp)==rec['input_sha256']==summary['source_sha256'],'Input hash mismatch')
        sources[str(inp.relative_to(REPO))]=sha(inp)
        original=read(inp);a=read(result/(rec['label']+'_PW_PBE.extxyz'))
        require(np.array_equal(original.numbers,a.numbers) and np.array_equal(original.arrays['lammps_id'],a.arrays['lammps_id']) and np.array_equal(original.arrays['central_pair'],a.arrays['central_pair']),'Identity mismatch')
        require(np.array_equal(original.pbc,a.pbc) and np.all(a.pbc) and np.allclose(original.cell,a.cell,atol=1e-12,rtol=0) and np.allclose(original.positions,a.positions,atol=1e-12,rtol=0),'Geometry mismatch')
        require(not mod.overlaps(train,[a]),'New training acquisition overlaps retained training geometry')
        e=float(a.info['PW_PBE_energy_eV']);f=np.asarray(a.arrays['PW_PBE_forces']);require(np.isfinite(e) and f.shape==(len(a),3) and np.isfinite(f).all() and abs(e-summary['energy_eV_cell'])<1e-10,'Invalid acquisition labels')
        require(summary['method']['xc']=='PBE' and summary['method']['cutoff_eV']==500 and summary['method']['kpts']==[1,1,1],'Method mismatch')
        a.calc=None;a.info['REF_energy']=e;a.arrays['REF_forces']=f.copy();a.info.pop('PW_PBE_energy_eV',None);a.arrays.pop('PW_PBE_forces',None);a.info['config_type']=rec['label']+'_PW_PBE_v16_train';a.info['dataset_role']='training_acquisition'
        train.append(a);additions.append({'label':rec['label'],'source':str(result.relative_to(REPO)),'source_sha256':sha(result/(rec['label']+'_PW_PBE.extxyz')),'energy_convention':summary.get('energy_convention'),'provenance_status':'verified_V16_PW_PBE_archive'})
require(len(train)==35,'Expected 29+6=35 training frames')
blind=PROJECT/'pbe_interface_v16_blind_holdouts';blind_manifest=json.loads((blind/'input_manifest.json').read_text());held=[]
for rec in blind_manifest['records']:
    inp=blind/rec['input'];require(sha(inp)==rec['input_sha256'],'Frozen blind input changed');a=read(inp)
    require(a.calc is None and 'REF_energy' not in a.info and 'REF_forces' not in a.arrays,'Blind input labeled')
    for split,frames in [('train',train),('valid',read(parent/'valid.extxyz',index=':')),('test',read(parent/'test.extxyz',index=':'))]:require(not mod.overlaps(frames,[a]),'Fresh blind input overlaps '+split)
    held.append({'label':rec['label'],'input_sha256':rec['input_sha256'],'role':'fresh_withheld_local_registry; excluded from training/validation/model selection'})
matrix=np.array([[a.get_chemical_symbols().count(e) for e in ['Ag','C','Si','Ti']] for a in train]);require(np.linalg.matrix_rank(matrix)==4,'Composition rank not4')
DATA.mkdir(exist_ok=True)
write(DATA/'train.extxyz',train,format='extxyz')
for split in ('valid','test'):shutil.copy2(parent/(split+'.extxyz'),DATA/(split+'.extxyz'))
for a in read(DATA/'train.extxyz',index=':'):require(np.isfinite(a.info['REF_energy']) and np.isfinite(a.arrays['REF_forces']).all(),'Serialized labels invalid')
m.update({'version':'V16 controlled six-acquisition numerical screening','split_sizes':{'train':35,'valid':2,'test':3},'energy_label_counts':{'train':35,'valid':2,'test':3},'force_label_counts':{'train':35,'valid':2,'test':3},'output_sha256':{s:sha(DATA/(s+'.extxyz')) for s in ('train','valid','test')},'energy_composition_matrix':{'elements':['Ag','C','Si','Ti'],'rank':4,'matrix':matrix.tolist()},'parent_V15_manifest_sha256':sha(parent/'dataset_manifest.json'),'V16_additions':additions,'V16_sources_sha256':sources,'frozen_holdout_INPUTS_only':held,'blind_label_content_read':False,'inherited_provenance_ledger_scope':'per_frame_provenance inherited records refer to parent V15; V16_additions describes six new rows'})
(DATA/'dataset_manifest.json').write_text(json.dumps(m,indent=2)+'\n')
(DATA/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+p.name+'\n' for p in sorted(DATA.iterdir()) if p.name!='SHA256SUMS.txt'))
print('V16 build PASS:35/2/3, rank4, six verified additions, fresh blind isolation')
