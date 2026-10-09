"""Registered, fixed-geometry Ag-layer separation GPAW controls."""
import hashlib,json,os,shutil,sys
from pathlib import Path
import numpy as np
from ase.io import read
from ase.geometry import find_mic
ROOT=Path(__file__).resolve().parent;REPO=ROOT.parents[3]
ENERGY='GPAW native extrapolated energy, force_consistent=False; free_energy recorded separately'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def inventory():
 seen=set()
 for line in (ROOT/'SHA256SUMS.txt').read_text().splitlines():
  digest,name=line.split('  ',1);p=ROOT/name
  assert p.resolve().is_relative_to(ROOT.resolve()) and name not in seen and p.is_file() and sha(p)==digest,'Entry inventory mismatch: '+name
  seen.add(name)
 actual={str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file() and p.name!='SHA256SUMS.txt' and '__pycache__' not in p.parts and 'calculations' not in p.relative_to(ROOT).parts}
 assert actual==seen,'Entry inventory incomplete/extraneous'
def validate(m,r):
 assert r['pbc']==[True,True,False] and r['poissonsolver']=={'dipolelayer':'xy'}
 assert r['method']==dict(m['method'],pbc=r['pbc'],poissonsolver=r['poissonsolver'])
 assert r['method']['kpts']==[6,6,1] and r['method']['smearing_eV']==.1 and r['method']['cutoff_eV']==500 and r['method']['xc']=='PBE'
 assert m['energy_convention']==ENERGY
 source=REPO/r['parent_geometry_path'];assert source.resolve().is_relative_to(REPO.resolve()) and sha(source)==r['parent_geometry_sha256']
 src=read(source);dest=ROOT/r['input'];assert dest.resolve().is_relative_to(ROOT.resolve()) and sha(dest)==r['input_sha256'];a=read(dest)
 ag=np.array(src.get_chemical_symbols())=='Ag';assert ag.sum()==4 and src.arrays['central_pair'].sum()==2
 expected=src.copy();expected.cell[2,2]=44.2058765;expected.positions[:,2]+=5.;expected.positions[ag,2]+=r['ag_layer_shift_z_A'];expected.set_pbc(r['pbc'])
 assert a.get_chemical_symbols()==expected.get_chemical_symbols()==r['ordered_symbols']
 assert np.array_equal(a.arrays['lammps_id'],expected.arrays['lammps_id']) and a.arrays['lammps_id'].tolist()==r['ordered_ids']
 assert np.array_equal(a.arrays['central_pair'],expected.arrays['central_pair']) and a.arrays['central_pair'].tolist()==r['pair_markings']
 assert np.allclose(a.positions,expected.positions,rtol=0,atol=1e-12) and np.allclose(a.cell.array,expected.cell.array,rtol=0,atol=1e-12) and a.pbc.tolist()==r['pbc']
 assert not any(k in a.info for k in ('REF_energy','energy','PW_PBE_energy_eV')) and not any(k in a.arrays for k in ('forces','REF_forces','PW_PBE_forces'))
 fw=np.array(src.get_chemical_symbols())!='Ag';dist=[]
 for i in np.flatnonzero(ag):
  for j in np.flatnonzero(fw):
   v,_=find_mic(a.positions[j]-a.positions[i],a.cell,pbc=[True,True,False]);dist.append(float(np.linalg.norm(v)))
 assert min(dist)>0.8,'Unsafe Ag-framework overlap'
 return a
def load(label,launch=False):
 m=json.loads((ROOT/'input_manifest.json').read_text());r=next(x for x in m['records'] if x['label']==label)
 if launch:assert m['launch_enabled'] is True and r['launch_enabled'] is True,'Launch disabled pending A target review'
 inventory();validate(m,r)
 if launch:
  assert r.get('owner_task') in ('window-a','window-b') and r.get('owner_instance')
  sys.path.insert(0,str(REPO/'research/mace-v12-transfer/coordination'));import sync_tasks as sync
  task=sync.load_task(r['owner_task']);sync.ensure_owner(r['owner_task'],task);assert r['owner_instance']==sync.identity() and label in task['labels']
  quota=Path('/sys/fs/cgroup/cpu.max').read_text().split();cores=len(os.sched_getaffinity(0))
  if int(os.environ.get('OMPI_COMM_WORLD_SIZE','1'))>1:
   assert int(os.environ['OMPI_COMM_WORLD_SIZE'])==4;cpus=Path('/sys/fs/cgroup/cpuset.cpus.effective').read_text().strip();cores=sum((int(y.split('-')[1])-int(y.split('-')[0])+1) if '-' in y else 1 for y in cpus.split(','))
  if quota[0]!='max':cores=min(cores,int(quota[0])//int(quota[1]))
  assert cores>=4 and shutil.disk_usage(ROOT).free>=m['disk_budget']['minimum_start_bytes']
  import importlib.metadata as meta,gpaw_data
  assert meta.version('gpaw')==r['method']['GPAW_version'] and meta.version('gpaw-data')==r['method']['PAW_data_version'] and meta.version('ase')==r['method']['ASE_version']
  for name,digest in r['method']['PAW_setup_sha256'].items():assert sha(Path(gpaw_data.datapath())/name)==digest
 return m,r
