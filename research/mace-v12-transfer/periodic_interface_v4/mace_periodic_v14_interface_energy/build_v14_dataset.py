#!/usr/bin/env python3
"""Build v14 only from completed, independently checked acquisition labels."""
from hashlib import sha256
from io import StringIO
import json
from pathlib import Path
import shutil
import numpy as np
from ase.io import read, write

ROOT=Path(__file__).resolve().parent
PROJECT=ROOT.parent
PARENT=PROJECT/'mace_periodic_v13_interface_energy'
DATA=ROOT/'data'
if DATA.exists() and any(DATA.iterdir()):
    raise SystemExit('v14 data already exists; preserve it and inspect its manifest.')
if not (PARENT/'results/v13_independent_holdout_comparison.json').is_file():
    raise SystemExit('Score all v13 acquisition labels before assigning them to v14 training.')
sources=[
 (PROJECT/'pbe_interface_v13_holdouts', 'AgC_lateral_registry_holdout_v13'),
 (PROJECT/'pbe_interface_v13_holdouts', 'AgSi_lateral_registry_holdout_v13'),
 (PROJECT/'pbe_interface_v14_parallel_acquisition', 'AgTi_registry_probe_v13_01'),
 (PROJECT/'pbe_interface_v14_parallel_acquisition', 'AgC_registry_strain_acq_v14_01'),
 (PROJECT/'pbe_interface_v14_parallel_acquisition', 'AgSi_registry_strain_acq_v14_01'),
]
parent_train=(PARENT/'data/train.extxyz').read_bytes()
parent_valid=PARENT/'data/valid.extxyz'; parent_test=PARENT/'data/test.extxyz'
valid=read(parent_valid,index=':'); test=read(parent_test,index=':')
frames=[]; records=[]
for folder,label in sources:
    result_dir=folder/'calculations'/label
    verification=json.loads((result_dir/'verification.json').read_text())
    summary=json.loads((result_dir/'summary.json').read_text())
    progress=json.loads((result_dir/'progress.json').read_text())
    if verification.get('status')!='PASS' or summary.get('scf_converged') is not True or progress.get('status')!='complete':
        raise SystemExit('Unverified or unconverged DFT source: '+label)
    output=next(result_dir.glob('*_PW_PBE.extxyz'))
    expected=verification['sha256'][output.name]
    if sha256(output.read_bytes()).hexdigest()!=expected: raise SystemExit('Output hash mismatch: '+label)
    atoms=read(output)
    energy=float(atoms.info['PW_PBE_energy_eV'])
    forces=np.asarray(atoms.arrays['PW_PBE_forces'],dtype=float)
    if not np.all(atoms.pbc) or not np.isfinite(energy) or forces.shape!=(len(atoms),3) or not np.isfinite(forces).all():
        raise SystemExit('Bad labels or nonperiodic geometry: '+label)
    for old in read(PARENT/'data/train.extxyz',index=':'):
        if (len(old)==len(atoms) and old.get_chemical_symbols()==atoms.get_chemical_symbols()
                and np.allclose(old.cell.array,atoms.cell.array,atol=1e-12,rtol=0)
                and np.allclose(old.positions,atoms.positions,atol=1e-10,rtol=0)):
            raise SystemExit('A newly assigned training frame duplicates a v13 training geometry: '+label)
    atoms.info['REF_energy']=energy; atoms.arrays['REF_forces']=forces.copy()
    atoms.info['config_type']=label+'_periodic_PW_PBE_energy_force_v14_train'
    atoms.info['source_method']='GPAW 26.7.0 PW-PBE 500 eV Gamma single point'
    atoms.info['previous_role']='v13-independent geometry check; scored before entering v14 training'
    atoms.info.pop('PW_PBE_energy_eV',None); atoms.arrays.pop('PW_PBE_forces',None); atoms.calc=None
    frames.append(atoms)
    records.append({'label':label,'source_directory':str(result_dir.relative_to(PROJECT)),
                    'output_sha256':expected,'summary_sha256':sha256((result_dir/'summary.json').read_bytes()).hexdigest(),
                    'atoms':len(atoms),'energy_eV_cell':energy,'scf_iterations':summary['scf_iterations'],
                    'central_pair':summary['central_pair']})
buffer=StringIO(); write(buffer,frames,format='extxyz')
DATA.mkdir(parents=True,exist_ok=True)
prefix=parent_train if parent_train.endswith(b'\n') else parent_train+b'\n'
(DATA/'train.extxyz').write_bytes(prefix+buffer.getvalue().encode())
shutil.copy2(parent_valid,DATA/'valid.extxyz'); shutil.copy2(parent_test,DATA/'test.extxyz')
train=read(DATA/'train.extxyz',index=':')
if len(train)!=31 or len(valid)!=2 or len(test)!=3: raise SystemExit('Unexpected v14 split sizes.')
if (DATA/'valid.extxyz').read_bytes()!=parent_valid.read_bytes() or (DATA/'test.extxyz').read_bytes()!=parent_test.read_bytes():
    raise SystemExit('Frozen v13 valid/test must be preserved byte-for-byte.')
for atoms in train[-len(frames):]:
    if not ('REF_energy' in atoms.info and 'REF_forces' in atoms.arrays): raise SystemExit('New frame lacks both DFT labels.')
energy_frames=[a for a in train if 'REF_energy' in a.info]
elements=['Ag','C','Si','Ti']
matrix=np.asarray([[a.get_chemical_symbols().count(e) for e in elements] for a in energy_frames],dtype=float)
rank=int(np.linalg.matrix_rank(matrix))
if rank!=4: raise SystemExit('Energy composition rank is '+str(rank)+'/4.')
holdout_records=[]
holdout_sources=[(PROJECT/'pbe_interface_v14_main_holdouts',x) for x in ('AgC_registry_holdout_v14_01','AgSi_registry_holdout_v14_01')]
holdout_sources += [(PROJECT/'pbe_interface_v14_parallel_acquisition','AgTi_registry_holdout_v14_01')]
for folder,label in holdout_sources:
    result_dir=folder/'calculations'/label
    v=json.loads((result_dir/'verification.json').read_text())
    if v.get('status')!='PASS': raise SystemExit('v14 holdout is not yet archived and verified: '+label)
    f=next(result_dir.glob('*_PW_PBE.extxyz'))
    if sha256(f.read_bytes()).hexdigest()!=v['sha256'][f.name]: raise SystemExit('v14 holdout hash mismatch: '+label)
    held=read(f)
    for used in train:
        if (len(held)==len(used) and held.get_chemical_symbols()==used.get_chemical_symbols()
                and np.allclose(held.cell.array,used.cell.array,atol=1e-12,rtol=0)
                and np.allclose(held.positions,used.positions,atol=1e-10,rtol=0)):
            raise SystemExit('A reserved v14 holdout duplicates training geometry: '+label)
    holdout_records.append({'label':label,'path':str(f.relative_to(PROJECT)),
                            'sha256':sha256(f.read_bytes()).hexdigest(),
                            'role':'blind_v14_holdout; exclude from training and checkpoint selection'})
manifest={'model_version':'MACE v14 Ag/Ti/Si/C; five checked v13 acquisition labels added',
          'parent_v13_model_sha256':sha256((PARENT/'checkpoints/MACE_periodic_v13_interface_energy_run-45.model').read_bytes()).hexdigest(),
          'parent_train_sha256':sha256((PARENT/'data/train.extxyz').read_bytes()).hexdigest(),
          'split_sizes':{'train':len(train),'valid':len(valid),'test':len(test)},
          'energy_label_counts':{'train':len(energy_frames),'valid':sum('REF_energy' in a.info for a in valid),'test':sum('REF_energy' in a.info for a in test)},
          'force_label_counts':{'train':sum('REF_forces' in a.arrays for a in train),'valid':sum('REF_forces' in a.arrays for a in valid),'test':sum('REF_forces' in a.arrays for a in test)},
          'energy_composition_matrix':{'elements':elements,'rank':rank},
          'added_training_labels':records,
          'reserved_blind_holdouts':holdout_records,
          'frozen_v13_validation_and_test_byte_identical':True,
          'output_sha256':{name:sha256((DATA/f'{name}.extxyz').read_bytes()).hexdigest() for name in ('train','valid','test')},
          'scope_note':'Five registry/local-noise acquisition checks are derived from previous cluster motifs; three additional labels are held out for v14 scoring. This does not by itself represent extended bulk interfaces, liquid Ag or production-temperature structures.'}
(DATA/'dataset_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
(DATA/'README.md').write_text('# MACE v14 data\n\nFive converged, hash-checked labels previously scored against v13 were added to its 26-frame training set. The three reserved `v14_holdout` labels remain outside training, validation and checkpoint selection. The v13 frozen validation/test files are copied byte-for-byte. See `dataset_manifest.json`.\n')
print(json.dumps(manifest,indent=2))
