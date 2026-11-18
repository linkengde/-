#!/workspace/venvs/dft-labels/bin/python
"""Final geometry/force/short-contact audit for the v12g 458k FIRE screen."""
from pathlib import Path
import json
import numpy as np
from scipy.spatial import cKDTree
from collections import Counter

HERE=Path(__file__).resolve().parent
OUT=HERE/'fire_v12g_458k'
MODEL={'1':'Ag','2':'C','3':'Si','4':'Ti'}
SRC={'1':'Ti','2':'Si','3':'C','4':'Ag'}

def read_atoms_dump(path):
    lines=Path(path).read_text().splitlines()
    h=next(i for i,s in enumerate(lines) if s.startswith('ITEM: ATOMS'))
    cols=lines[h].split()[2:]
    rows=np.loadtxt(lines[h+1:],dtype=float)
    return cols,rows

def read_data(path,n):
    lines=Path(path).read_text().splitlines()
    a=next(i for i,s in enumerate(lines) if s.startswith('Atoms'))+2
    return np.asarray([[float(q) for q in s.split()[:5]] for s in lines[a:a+n]],dtype=float)

cols,r=read_atoms_dump(OUT/'forces_minimized.dump')
ix={k:cols.index(k) for k in ('id','type','x','y','z','fx','fy','fz')}
ids=r[:,ix['id']].astype(np.int64); types=r[:,ix['type']].astype(np.int8)
pos=r[:,[ix['x'],ix['y'],ix['z']]]; forces=r[:,[ix['fx'],ix['fy'],ix['fz']]]
n=len(ids); norms=np.linalg.norm(forces,axis=1)
coords=pos.copy();coords[:,0]%=200;coords[:,1]%=200;coords[:,2]+=100
D,NN=cKDTree(coords,boxsize=[200,200,1000]).query(coords,k=2)
near=D[:,1]; nn=NN[:,1]; short=near<1.8
pair_counts=Counter(); atom_counts=Counter(); unique=set()
for i in np.where(short)[0]:
    a,b=sorted((MODEL[str(int(types[i]))],MODEL[str(int(types[nn[i]]))]))
    atom_counts[f'{a}-{b}']+=1
    unique.add(tuple(sorted((int(ids[i]),int(ids[nn[i]])))))
for a,b in unique:
    ia=np.searchsorted(ids,a); ib=np.searchsorted(ids,b)
    x,y=sorted((MODEL[str(int(types[ia]))],MODEL[str(int(types[ib]))]))
    pair_counts[f'{x}-{y}']+=1
src=read_data(HERE/'Ag_Ti3SiC2_bicontinuous_200x200x300.data',n)
src_ids=src[:,0].astype(int); src_types=src[:,1].astype(int); init=src[:,2:5]
order=np.argsort(ids); sorted_ids=ids[order]
if not np.array_equal(sorted_ids,src_ids):
    raise RuntimeError('Atom ID set differs from source geometry')
rows_for_src=order[np.searchsorted(sorted_ids,src_ids)]
final_in_src=pos[rows_for_src]
delta=final_in_src-init;delta[:,0]-=200*np.rint(delta[:,0]/200);delta[:,1]-=200*np.rint(delta[:,1]/200)
disp=np.linalg.norm(delta,axis=1)
fixed=(src_types!=4)&(init[:,2]>=0)&(init[:,2]<=20)
fixed_disp=disp[fixed]
maxi=int(np.argmax(norms)); mini=int(np.argmin(near))
result={
 'status':'AUDITED','step':100,'natoms':int(n),'unique_ids':int(len(np.unique(ids))),
 'nonfinite_positions':int((~np.isfinite(pos)).any(axis=1).sum()),'nonfinite_forces':int((~np.isfinite(forces)).any(axis=1).sum()),
 'force_max_eV_A':float(norms[maxi]),'force_p99_eV_A':float(np.percentile(norms,99)),'force_p95_eV_A':float(np.percentile(norms,95)),
 'atoms_force_above_500_eV_A':int((norms>500).sum()),'max_force_atom_id':int(ids[maxi]),'max_force_model_type':MODEL[str(int(types[maxi]))],
 'min_neighbor_A':float(near[mini]),'min_neighbor_pair_ids':[int(ids[mini]),int(ids[nn[mini]])],
 'nearest_neighbor_pairs_below_1p8_A':int(len(unique)),'atom_counts_nearest_below_1p8_A':int(short.sum()),
 'short_pair_types_by_atom':dict(atom_counts),'short_pair_types_unique':dict(pair_counts),
 'max_displacement_A':float(disp.max()),'p99_displacement_A':float(np.percentile(disp,99)),
 'atoms_displaced_more_than_1A':int((disp>1).sum()),'atoms_displaced_more_than_2A':int((disp>2).sum()),
 'fixed_bottom_skeleton_atoms':int(fixed.sum()),'fixed_bottom_skeleton_max_displacement_A':float(fixed_disp.max()),
 'fixed_bottom_skeleton_p99_displacement_A':float(np.percentile(fixed_disp,99)),
 'xyz_min_A':pos.min(axis=0).tolist(),'xyz_max_A':pos.max(axis=0).tolist(),
 'atoms_z_above_200_A':int((pos[:,2]>200).sum()),'atoms_z_below_0_A':int((pos[:,2]<0).sum()),
 'interpretation':'FIRE hit its 100-iteration limit. Do not use minimized coordinates for MD: extensive Ti-C and Ag-skeleton short contacts indicate collapse.'}
(OUT/'coordinate_force_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
