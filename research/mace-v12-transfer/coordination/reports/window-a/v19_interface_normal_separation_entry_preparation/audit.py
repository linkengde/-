import hashlib,json
from pathlib import Path
import numpy as np
from ase.io import read
from ase.geometry import find_mic
ROOT=Path(__file__).resolve().parent;REPO=ROOT.parents[5]
m=json.loads((ROOT/'manifest.json').read_text());assert len(m['cases'])==9
for r in m['cases']:
 src=REPO/r['source'];out=ROOT/r['proposal'];assert hashlib.sha256(src.read_bytes()).hexdigest()==r['source_sha256'];assert hashlib.sha256(out.read_bytes()).hexdigest()==r['proposal_sha256']
 a=read(src);b=read(out);assert a.get_chemical_symbols()==b.get_chemical_symbols() and np.array_equal(a.arrays['lammps_id'],b.arrays['lammps_id']) and np.array_equal(a.arrays['central_pair'],b.arrays['central_pair']) and np.array_equal(a.pbc,b.pbc) and np.allclose(a.cell,b.cell,rtol=0,atol=1e-12)
 ag=np.array(a.get_chemical_symbols())=='Ag';delta=b.positions-a.positions
 assert np.allclose(delta[~ag],0,atol=1e-12) and np.allclose(delta[ag,0:2],0,atol=1e-12) and np.allclose(delta[ag,2],r['Ag_layer_shift_z_A'],atol=1e-12)
 assert r['launch_enabled'] is False and r['owner'] is None
 d=[]
 for i in np.flatnonzero(ag):
  for j in np.flatnonzero(~ag):
   v,_=find_mic(b.positions[j]-b.positions[i],b.cell,pbc=[True,True,False]);d.append(np.linalg.norm(v))
 assert abs(min(d)-r['minimum_interfacial_distance_A'])<1e-8 and min(d)>r['overlap_floor_A']
print('PASS: 9 source/proposal records; hashes, species/order/IDs/cell/PBC, rigid layer movement and distance screens')
