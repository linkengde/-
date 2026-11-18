#!/usr/bin/env python3
"""Stage an auditable, non-official v19+v25 MTP/ACE candidate dataset."""
from __future__ import annotations
import hashlib, json
from pathlib import Path
import numpy as np
from ase.io import read, write

ROOT = Path('/workspace/-/research/mace-v12-transfer/periodic_interface_v4')
V19 = ROOT/'mace_periodic_v19_residual_response/data'
V25 = ROOT/'pbe_interface_v25_tisic_skeleton_acquisition/dataset'
OUT = Path('/workspace/potential-comparison/mace_fit_data_20261118')
SYMBOLS = ['Ag', 'C', 'Si', 'Ti']
TYPE = {s:i for i,s in enumerate(SYMBOLS)}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def geom_sha(a):
    h=hashlib.sha256()
    h.update(' '.join(a.get_chemical_symbols()).encode())
    h.update(np.asarray(a.cell, dtype='<f8').round(8).tobytes())
    h.update(np.asarray(a.positions, dtype='<f8').round(8).tobytes())
    h.update(np.asarray(a.pbc, dtype=np.uint8).tobytes())
    return h.hexdigest()

def assert_data(a, source):
    assert 'REF_energy' in a.info, f'missing REF_energy: {source}'
    assert 'REF_forces' in a.arrays, f'missing REF_forces: {source}'
    assert np.asarray(a.arrays['REF_forces']).shape == (len(a),3), source
    assert np.isfinite(float(a.info['REF_energy']))
    assert np.isfinite(a.arrays['REF_forces']).all()
    assert set(a.get_chemical_symbols()) <= set(SYMBOLS), source
    assert np.linalg.det(np.asarray(a.cell)) > 0, source

# Input datasets are hash-checked before any transformation.
for data_dir in [V19, V25]:
    import subprocess
    subprocess.run(['sha256sum','-c','SHA256SUMS.txt'],cwd=data_dir,check=True,stdout=subprocess.PIPE,text=True)

v19_train=read(V19/'train.extxyz',':')
v19_dev=read(V19/'valid.extxyz',':')
v25_train=read(V25/'train.extxyz',':')
v25_dev=read(V25/'valid.extxyz',':')
assert (len(v19_train),len(v19_dev),len(v25_train),len(v25_dev)) == (36,3,19,4)

native_free_differences=[]
for i,a in enumerate(v25_train+v25_dev):
    # Preserve all raw labels in metadata; use native extrapolated energy to
    # match the A-approved v19 dataset's declared energy convention.
    native=a.info.get('PW_PBE_energy_eV')
    free=a.info.get('PW_PBE_free_energy_eV')
    assert native is not None and free is not None, f'v25 energy keys missing at frame {i}'
    native_free_differences.append(float(native)-float(free))
    a.info['REF_energy_original_free_eV']=float(a.info['REF_energy'])
    a.info['REF_energy']=float(native)
    a.info['energy_target_convention']='GPAW native extrapolated total energy (force_consistent=False); forces from sigma=0.10 eV calculation'
    a.info['mtp_ace_candidate_energy_rekey']='PW_PBE_energy_eV'
    assert_data(a,f'v25 frame {i}')
for i,a in enumerate(v19_train+v19_dev):
    assert_data(a,f'v19 frame {i}')

all_sets={'v19_train':v19_train,'v19_development':v19_dev,'v25_train':v25_train,'v25_development':v25_dev}
hashes={k:[geom_sha(a) for a in v] for k,v in all_sets.items()}
intersections={}
keys=list(hashes)
for i,k1 in enumerate(keys):
  for k2 in keys[i+1:]:
    overlap=sorted(set(hashes[k1]) & set(hashes[k2]))
    if overlap: intersections[f'{k1}__{k2}']=overlap

train=v19_train+v25_train
dev=v19_dev+v25_dev
write(OUT/'train.extxyz',train,format='extxyz')
write(OUT/'dev.extxyz',dev,format='extxyz')

def write_cfg(path, frames):
    with path.open('w') as f:
      for a in frames:
        f.write('BEGIN_CFG\n Size\n    %d\n Supercell\n' % len(a))
        for row in np.asarray(a.cell): f.write('    % .15f % .15f % .15f\n' % tuple(row))
        f.write(' AtomData:  id type cartes_x cartes_y cartes_z fx fy fz\n')
        forces=np.asarray(a.arrays['REF_forces'])
        for i,(sym,pos,force) in enumerate(zip(a.get_chemical_symbols(),a.positions,forces),1):
            f.write('    %d %d % .15f % .15f % .15f % .15f % .15f % .15f\n' % (i,TYPE[sym],*pos,*force))
        f.write(' Energy\n    % .15f\n' % float(a.info['REF_energy']))
        label=str(a.info.get('config_type','unnamed')).replace('\n',' ')
        f.write(' Feature   config_type %s\nEND_CFG\n\n' % label)
write_cfg(OUT/'train.cfg',train)
write_cfg(OUT/'dev.cfg',dev)

# Summary of role, geometry, labels, DFT convention, and force range.
def summary(name, frames):
  formulas={}
  for a in frames: formulas[a.get_chemical_formula(mode='hill')]=formulas.get(a.get_chemical_formula(mode='hill'),0)+1
  maxf=[float(np.linalg.norm(a.arrays['REF_forces'],axis=1).max()) for a in frames]
  mins=[]
  for a in frames:
    d=a.get_all_distances(mic=True)
    d[d<1e-8]=np.inf
    mins.append(float(d.min()))
  return {'frames':len(frames),'atom_count_total':sum(len(a) for a in frames),'atom_count_range':[min(len(a) for a in frames),max(len(a) for a in frames)],'symbols':sorted({s for a in frames for s in a.get_chemical_symbols()}),'formulas':formulas,'pbc_values':sorted({tuple(bool(x) for x in a.pbc) for a in frames}),'all_have_energy_and_forces':True,'max_force_eV_A_range':[min(maxf),max(maxf)],'minimum_pair_distance_A_range':[min(mins),max(mins)]}
manifest={
 'status':'EXPLORATORY_MTP_ACE_CANDIDATE_NOT_OFFICIAL',
 'target_energy_convention':'GPAW native extrapolated total energy (force_consistent=False); forces from Fermi-Dirac sigma=0.10 eV',
 'source_datasets':{
   'v19':{'path':str(V19),'status':'A_APPROVED_V19_DATASET','train_sha256':sha(V19/'train.extxyz'),'dev_sha256':sha(V19/'valid.extxyz'),'energy_convention':'native extrapolated total energy'},
   'v25':{'path':str(V25),'status':'non_official_skeleton_train_and_development_candidate','train_sha256':sha(V25/'train.extxyz'),'dev_sha256':sha(V25/'valid.extxyz'),'original_energy_convention':'finite-width free energy at sigma=0.10 eV','staged_energy_key':'PW_PBE_energy_eV (native extrapolated total energy)','native_minus_free_eV_range':[min(native_free_differences),max(native_free_differences)]}},
 'roles':{'train':summary('train',train),'development':summary('dev',dev),'components':{k:summary(k,v) for k,v in all_sets.items()}},
 'geometry_hash_intersections':intersections,
 'limitations':['No blind/test set exists.','The three v19 development frames share one COD 9009647 source family.','All v25 train/development frames descend from a single Ti3SiC2 crystallographic source family; v25 development is correlated with its training family.','v25 native extrapolated energies are used to match v19, while forces remain the sigma=0.10 eV GPAW forces; this thermodynamic energy/force convention is not exactly variationally consistent.','MTP/ACE fit gates and Si-hotspot center errors must be evaluated on the fixed development set; no production or MD qualification is implied.'],
 'files_sha256':{},
}
for name in ['train.extxyz','dev.extxyz','train.cfg','dev.cfg']:
    manifest['files_sha256'][name]=sha(OUT/name)
(OUT/'candidate_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps({k:manifest[k] for k in ['status','target_energy_convention','roles','geometry_hash_intersections','limitations','files_sha256']},indent=2))
