"""Shared static and launch checks; no calculator construction."""
import hashlib,json,os,shutil,sys
from pathlib import Path
import numpy as np
from ase.io import read
ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[5]
ENERGY='GPAW native extrapolated energy, force_consistent=False; free_energy recorded separately'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def inventory():
 seen=set()
 for line in (ROOT/'SHA256SUMS.txt').read_text().splitlines():
  digest,name=line.split('  ',1);p=ROOT/name
  assert p.resolve().is_relative_to(ROOT.resolve()) and name not in seen
  assert p.is_file() and sha(p)==digest,'Entry inventory mismatch: '+name
  seen.add(name)
 actual={str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file() and p.name!='SHA256SUMS.txt' and '__pycache__' not in p.parts and 'calculations' not in p.relative_to(ROOT).parts}
 assert actual==seen,'Entry inventory incomplete/extraneous'
def validate(manifest,record):
 assert record['pbc']==[True,True,False],'Slab must be nonperiodic z'
 assert record['poissonsolver'] in ({},{'dipolelayer':'xy'}),'Wrong Poisson setting'
 assert record['method']==dict(manifest['method'],pbc=record['pbc'],poissonsolver=record['poissonsolver'])
 assert record['method']['kpts']==[6,6,1] and record['method']['smearing_eV']==.1
 assert record['method']['cutoff_eV']==500 and record['method']['xc']=='PBE'
 assert manifest['energy_convention']==ENERGY
 src=ROOT/record['input'];assert src.resolve().is_relative_to(ROOT.resolve()) and sha(src)==record['input_sha256'],'Input hash mismatch'
 a=read(src);assert a.pbc.tolist()==record['pbc']
 assert np.isfinite(a.positions).all() and np.isfinite(a.cell.array).all() and np.linalg.det(a.cell)>0
 assert np.allclose(a.cell.array[2,:2],0,rtol=0,atol=1e-12) and np.allclose(a.cell.array[:2,2],0,rtol=0,atol=1e-12)
 assert a.get_chemical_symbols()==record['ordered_symbols']
 assert a.arrays['lammps_id'].tolist()==record['ordered_ids'] and a.arrays['central_pair'].tolist()==record['pair_markings']
 assert not any(k in a.info for k in ('REF_energy','energy','PW_PBE_energy_eV')) and not any(k in a.arrays for k in ('forces','REF_forces','PW_PBE_forces'))
 source=REPO/manifest['source_geometry_path'];assert sha(source)==manifest['source_geometry_sha256']
 parent=read(source);assert np.array_equal(a.positions,parent.positions) and np.array_equal(a.cell.array,parent.cell.array)
 assert a.get_chemical_symbols()==parent.get_chemical_symbols()
 return a
def load(label,launch=False):
 m=json.loads((ROOT/'input_manifest.json').read_text());r=next(x for x in m['records'] if x['label']==label)
 if launch:
  assert m['launch_enabled'] is True and r['launch_enabled'] is True,'Launch disabled: A integration/registration required'
 inventory();validate(m,r)
 if launch:
  sys.path.insert(0,str(REPO/'research/mace-v12-transfer/coordination'));import sync_tasks as sync
  task=sync.load_task(r['owner_task']);sync.ensure_owner(r['owner_task'],task)
  assert r['owner_instance']==sync.identity() and label in task['labels']
  quota=Path('/sys/fs/cgroup/cpu.max').read_text().split();cores=len(os.sched_getaffinity(0))
  if quota[0]!='max':cores=min(cores,int(quota[0])//int(quota[1]))
  assert cores>=4 and shutil.disk_usage(ROOT).free>=m['disk_budget']['minimum_start_bytes']
  import importlib.metadata as meta, gpaw_data
  assert meta.version('gpaw')==r['method']['GPAW_version'] and meta.version('gpaw-data')==r['method']['PAW_data_version']
  for name,digest in r['method']['PAW_setup_sha256'].items():assert sha(Path(gpaw_data.datapath())/name)==digest
 return m,r
