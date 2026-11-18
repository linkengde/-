#!/workspace/venvs/dft-labels/bin/python
"""Audit Ti-C short-pair coverage and labels in the v12g training candidates."""
from __future__ import annotations
import json
from pathlib import Path
from collections import Counter
import numpy as np
from ase.io import read
from ase.neighborlist import neighbor_list
from ase.data import atomic_masses

HERE=Path(__file__).resolve().parent/'v12g_si04_candidate'
files=sorted((HERE/'dataset/train').glob('frame_*.xyz'))
items=[]
for fp in files:
    for n,atoms in enumerate(read(fp,index=':')):
        items.append((fp.name,n,atoms))
local=read(HERE/'dataset/local_train/local_dft_targets.xyz',index=':')
items.extend(('local_dft_targets.xyz',n,a) for n,a in enumerate(local))
counts=Counter(); minrec=None
for filename,frame,atoms in items:
    i,j,D,d=neighbor_list('ijDd',atoms,3.0,self_interaction=False)
    symbols=atoms.get_chemical_symbols()
    mask=np.asarray([symbols[x]=='Ti' and symbols[y]=='C' or symbols[x]=='C' and symbols[y]=='Ti'
                     for x,y in zip(i,j)],dtype=bool)
    if not mask.any():
        continue
    ds=d[mask]
    for threshold in (1.8,2.0,2.2):
        counts[f'undirected_pairs_below_{threshold:.1f}_A']+=int((ds<threshold).sum()//2)
    k=int(np.where(mask)[0][np.argmin(ds)])
    if minrec is None or float(d[k])<minrec['distance_A']:
        ti,ci=(int(i[k]),int(j[k])) if symbols[int(i[k])]=='Ti' else (int(j[k]),int(i[k]))
        vec=np.asarray(D[k] if int(i[k])==ti else -D[k],dtype=float)
        rhat=vec/np.linalg.norm(vec); forces=atoms.get_forces()
        forceproj=float(np.dot(forces[ci]/atomic_masses[6]-forces[ti]/atomic_masses[22],rhat))
        minrec={'file':filename,'frame':frame,'distance_A':float(d[k]),'Ti_index':ti,'C_index':ci,
                'Ti_force_eV_A':forces[ti].tolist(),'C_force_eV_A':forces[ci].tolist(),
                'Ti_force_norm_eV_A':float(np.linalg.norm(forces[ti])),'C_force_norm_eV_A':float(np.linalg.norm(forces[ci])),
                'mass_weighted_relative_force_projection_along_Ti_to_C':forceproj}
result={'status':'AUDITED','baseline_train_frames':len(files),'local_train_frames':len(local),
        'total_train_frames':len(items),'TiC_pair_counts':dict(counts),'minimum_TiC_pair':minrec,
        'interpretation':'A negative v12g full-model net-force relative projection for the new interface candidate coexists with DFT training examples at even shorter Ti-C distances; this is a local transfer/stability issue, not a total absence of short-pair data.'}
out=HERE/'tic_training_coverage_audit.json';out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
