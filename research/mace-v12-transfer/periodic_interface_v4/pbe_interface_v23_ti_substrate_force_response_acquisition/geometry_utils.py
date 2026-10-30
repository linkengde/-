"""Geometry-only extxyz reader and duplicate screen used for v22 inputs.

The reader consumes symbols, positions, cell, PBC and optional persistent IDs
only.  Its near-duplicate thresholds match the project geometry screen.
"""
from __future__ import annotations
import re
from collections import Counter
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.geometry import find_mic

NEAR_AG_RMS_A = 0.15
NEAR_FRAMEWORK_RMS_A = 0.15
EXACT_MAX_A = 1e-5
CELL_TOL_A = 1e-5

def extxyz_geometry_frames(path: Path):
    result=[]
    with Path(path).open('r',encoding='utf-8') as f:
        frame=0
        while True:
            first=f.readline()
            if not first: break
            if not first.strip(): continue
            n=int(first.strip()); header=f.readline().rstrip('\n')
            match=re.search(r'Properties=([^ ]+)',header)
            if match is None: raise ValueError(f'missing extxyz Properties in {path} frame {frame}')
            schema=match.group(1).split(':'); fields={}; offset=0
            for i in range(0,len(schema),3):
                name,kind,width=schema[i],schema[i+1],int(schema[i+2]); fields[name]=(offset,width,kind); offset+=width
            if 'species' not in fields or 'pos' not in fields: raise ValueError(f'missing geometry columns in {path} frame {frame}')
            lattice=re.search(r'Lattice="([^"]+)"',header)
            if lattice is None: raise ValueError(f'missing cell in {path} frame {frame}')
            cell=np.asarray([float(x) for x in lattice.group(1).split()],dtype=float).reshape(3,3)
            pbcm=re.search(r'pbc="([^"]+)"',header)
            pbc=[x.upper() in {'T','TRUE','1'} for x in pbcm.group(1).split()] if pbcm else [True,True,True]
            symbols=[]; positions=[]; ids=[]
            for _ in range(n):
                row=f.readline().split(); s0,_,_=fields['species']; p0,_,_=fields['pos']
                symbols.append(row[s0]); positions.append([float(x) for x in row[p0:p0+3]])
                if 'lammps_id' in fields:
                    q,_,_=fields['lammps_id']; ids.append(int(row[q]))
            atoms=Atoms(symbols=symbols,positions=positions,cell=cell,pbc=pbc)
            if ids: atoms.new_array('lammps_id',np.asarray(ids,dtype=int))
            result.append(atoms); frame+=1
    return result

def compare_geometry(a,b):
    if len(a)!=len(b) or Counter(a.get_chemical_symbols())!=Counter(b.get_chemical_symbols()) or not np.array_equal(a.pbc,b.pbc): return None
    cell_delta=float(np.max(np.abs(a.cell.array-b.cell.array)))
    if cell_delta>CELL_TOL_A: return None
    aid=a.arrays.get('lammps_id'); bid=b.arrays.get('lammps_id')
    if aid is not None and bid is not None and set(map(int,aid))==set(map(int,bid)):
        amap={int(x):i for i,x in enumerate(aid)}; bmap={int(x):i for i,x in enumerate(bid)}; ids=sorted(amap)
        ia=np.asarray([amap[x] for x in ids]); ib=np.asarray([bmap[x] for x in ids]); method='persistent_lammps_id'
    elif a.get_chemical_symbols()==b.get_chemical_symbols():
        ia=np.arange(len(a)); ib=np.arange(len(b)); method='same_element_order'
    else:
        return None
    delta,_=find_mic(a.positions[ia]-b.positions[ib],a.cell,pbc=a.pbc)
    symbols=np.asarray(a.get_chemical_symbols())[ia]; framework=symbols!='Ag'; ag=symbols=='Ag'
    common=np.median(delta[framework],axis=0) if framework.any() else np.median(delta,axis=0)
    magnitudes=np.linalg.norm(delta-common,axis=1)
    ag_rms=float(np.sqrt(np.mean(magnitudes[ag]**2))) if ag.any() else 0.0
    fw_rms=float(np.sqrt(np.mean(magnitudes[framework]**2))) if framework.any() else 0.0
    return {'matching_method':method,'cell_max_delta_A':cell_delta,'Ag_registry_RMS_A':ag_rms,
            'framework_RMS_A':fw_rms,'all_atom_RMS_A':float(np.sqrt(np.mean(magnitudes**2))),
            'max_atom_displacement_A':float(np.max(magnitudes)),
            'exact_duplicate':bool(cell_delta<=CELL_TOL_A and float(np.max(magnitudes))<=EXACT_MAX_A),
            'near_duplicate':bool(ag_rms<=NEAR_AG_RMS_A and fw_rms<=NEAR_FRAMEWORK_RMS_A)}
