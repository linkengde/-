import argparse, hashlib, itertools, json
from pathlib import Path
import numpy as np
from ase.io import read

p=argparse.ArgumentParser()
p.add_argument('--left',required=True)
p.add_argument('--right',required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args()
repo=Path('/workspace/-')
base=repo/'research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_pure_phase_controls'
manifest=json.loads((base/'input_manifest.json').read_text())
records={r['label']:r for r in manifest['records']}
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def load(label):
 r=records[label]; folder=base/'calculations'/label
 checked={}
 for line in (folder/'SHA256SUMS.txt').read_text().splitlines():
  digest,name=line.split('  ',1); f=folder/name
  if sha(f)!=digest: raise ValueError(f'SHA mismatch: {f}')
  checked[name]=digest
 v=json.loads((folder/'verification.json').read_text())
 s=json.loads((folder/'summary.json').read_text())
 if v.get('status')!='PASS' or not all(v.get('checks',{}).values()) or not s.get('scf_converged'):
  raise ValueError(f'Archive not converged/verified: {label}')
 if s.get('source_sha256')!=r['input_sha256'] or v.get('input_sha256')!=r['input_sha256']:
  raise ValueError(f'Input hash mismatch: {label}')
 if s['method']!=r['method'] or list(s['method']['kpts'])!=r['kpts']:
  raise ValueError(f'Method/mesh mismatch: {label}')
 atoms=read(folder/f'{label}_PW_PBE.extxyz',format='extxyz')
 if len(atoms)!=48 or atoms.get_chemical_formula()!='C16Si8Ti24': raise ValueError('Unexpected composition')
 return r,s,atoms,checked
ra,sa,aa,ha=load(a.left); rb,sb,ab,hb=load(a.right)
ma={k:v for k,v in sa['method'].items() if k!='kpts'}; mb={k:v for k,v in sb['method'].items() if k!='kpts'}
for name,x,y in [('method',ma,mb),('source hash',sa['source_sha256'],sb['source_sha256']),('symbols',aa.numbers,ab.numbers),('positions',aa.positions,ab.positions),('cell',aa.cell.array,ab.cell.array),('pbc',aa.pbc,ab.pbc),('IDs',aa.arrays['local_new_id'],ab.arrays['local_new_id'])]:
 equal=(x==y) if isinstance(x,dict) or isinstance(x,str) else np.array_equal(x,y)
 if not equal: raise ValueError(f'Identity mismatch: {name}')
df=np.asarray(ab.arrays['PW_PBE_forces'])-np.asarray(aa.arrays['PW_PBE_forces'])
vectors=np.linalg.norm(df,axis=1)
force={'vector_rms_eV_A':float(np.sqrt(np.mean(vectors**2))), 'maximum_atom_vector_difference_eV_A':float(vectors.max()), 'component_rms_xyz_eV_A':np.sqrt(np.mean(df**2,axis=0)).tolist(), 'species':{}}
for symbol in sorted(set(aa.get_chemical_symbols())):
 mask=np.array(aa.get_chemical_symbols())==symbol
 force['species'][symbol]={'n':int(mask.sum()),'vector_rms_eV_A':float(np.sqrt(np.mean(vectors[mask]**2))),'maximum_atom_vector_difference_eV_A':float(vectors[mask].max())}
n=len(aa)
delta_native=(sb['energy_eV_cell']-sa['energy_eV_cell'])*1000/n
delta_free=(sb['free_energy_eV_cell']-sa['free_energy_eV_cell'])*1000/n
report={'left':a.left,'right':a.right,'source_hash':sa['source_sha256'],'same_geometry_method_except_kpts':True,'archive_verification':{'left':'PASS','right':'PASS','left_files':ha,'right_files':hb},'native_energy_delta_meV_atom':delta_native,'free_energy_delta_meV_atom':delta_free,'force':force,'thresholds':{'energy_abs_meV_atom':2.0,'force_vector_rms_eV_A':0.01},'pass':{'native_energy':abs(delta_native)<=2.0,'free_energy':abs(delta_free)<=2.0,'force_vector_rms':force['vector_rms_eV_A']<=0.01}}
a.output.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
