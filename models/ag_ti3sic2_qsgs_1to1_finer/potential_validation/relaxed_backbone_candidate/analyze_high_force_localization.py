#!/usr/bin/env python3
"""Locate high-force atoms and compare interface versus same-phase force components."""
from __future__ import annotations
import gzip, json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

HERE=Path(__file__).resolve().parent
ELEMENT={1:'Ti',2:'Si',3:'C',4:'Ag'}
GEOM=json.loads((HERE/'geometry_report.json').read_text())
BOX=np.asarray(GEOM['simulation_box_A'],float)
MATERIAL_Z=float(GEOM['material_dimensions_A'][2])
TREE_BOX=[BOX[0],BOX[1],1000.]

def load_dump(name):
    with gzip.open(HERE/name,'rt') as stream: return np.loadtxt(stream,skiprows=9)
def stats(values):
    if not len(values): return None
    return {'n':int(len(values)),'p50':float(np.percentile(values,50)),'p95':float(np.percentile(values,95)),'p99':float(np.percentile(values,99)),'max':float(np.max(values))}

full=load_dump('forces_after300_FULL.dump.gz'); no_cross=load_dump('forces_after300_NO_CROSS.dump.gz')
if not np.array_equal(full[:,:2],no_cross[:,:2]): raise RuntimeError('ID/type mismatch')
if not np.allclose(full[:,2:5],no_cross[:,2:5],atol=1e-10,rtol=0): raise RuntimeError('Coordinate mismatch')
type_id=full[:,1].astype(int); pos=full[:,2:5].copy(); pos[:,0]%=BOX[0]; pos[:,1]%=BOX[1]
fvec=full[:,5:8]; nvec=no_cross[:,5:8]; cvec=fvec-nvec
fmag=np.linalg.norm(fvec,axis=1); nmag=np.linalg.norm(nvec,axis=1); cmag=np.linalg.norm(cvec,axis=1)
ag_indices=np.flatnonzero(type_id==4); sk_indices=np.flatnonzero(type_id!=4)
ag_tree=cKDTree(pos[ag_indices],boxsize=TREE_BOX); sk_tree=cKDTree(pos[sk_indices],boxsize=TREE_BOX); all_tree=cKDTree(pos,boxsize=TREE_BOX)
near_s_dist,near_s_local=sk_tree.query(pos[ag_indices],k=1); near_a_dist,near_a_local=ag_tree.query(pos[sk_indices],k=1)
d_opposite=np.empty(len(full)); opposite_index=np.empty(len(full),dtype=int)
d_opposite[ag_indices]=near_s_dist; opposite_index[ag_indices]=sk_indices[near_s_local]
d_opposite[sk_indices]=near_a_dist; opposite_index[sk_indices]=ag_indices[near_a_local]

def neighbor_counts(i):
    local=np.asarray([j for j in all_tree.query_ball_point(pos[i],r=8.0) if j!=i],dtype=int)
    delta=pos[local]-pos[i]; delta[:,0]-=BOX[0]*np.rint(delta[:,0]/BOX[0]); delta[:,1]-=BOX[1]*np.rint(delta[:,1]/BOX[1]); dist=np.linalg.norm(delta,axis=1)
    out={}
    for cutoff in (3.0,3.25,3.5,4.0,8.0):
        inside=dist<=cutoff
        out[f'within_{cutoff:g}_A']={ELEMENT[t]:int(np.sum(type_id[local[inside]]==t)) for t in (1,2,3,4)}
    return out

def top_atoms(indices,count=10):
    order=indices[np.argsort(fmag[indices])[-count:][::-1]]; result=[]
    for i in order:
        j=int(opposite_index[i]); pool=ag_indices if type_id[i]==4 else np.flatnonzero(type_id==type_id[i]); tree=ag_tree if type_id[i]==4 else cKDTree(pos[pool],boxsize=TREE_BOX)
        ds,js=tree.query(pos[i],k=min(2,len(pool)))
        if np.ndim(ds)==0: same=None
        else:
            k=1 if int(pool[js[0]])==i else 0; same={'element':ELEMENT[int(type_id[pool[js[k]]])],'id':int(full[pool[js[k]],0]),'distance_A':float(ds[k])}
        result.append({'id':int(full[i,0]),'element':ELEMENT[int(type_id[i])],'position_A':[float(x) for x in full[i,2:5]],'distance_to_top_surface_plane_A':float(MATERIAL_Z-full[i,4]),'distance_to_bottom_plane_A':float(full[i,4]),'force_full_eV_A':float(fmag[i]),'force_without_cross_eV_A':float(nmag[i]),'cross_force_component_eV_A':float(cmag[i]),'opposite_phase_nearest':{'element':ELEMENT[int(type_id[j])],'id':int(full[j,0]),'distance_A':float(d_opposite[i])},'nearest_same_element_neighbor':same,'neighbor_counts_by_element':neighbor_counts(int(i))})
    return result

report={'status':'RUN0_LOCAL_FORCE_DIAGNOSTIC_ONLY','atoms':int(len(full)),'force_definition':'Full force vector at current calibration-v0 coordinates; no-cross vector at identical coordinates after setting Ag-Ti/Si/C Morse D0=0; vector difference estimates cross-Morse contribution.','top_surface_plane_z_A':MATERIAL_Z,'by_element':{}}
for t,name in ELEMENT.items():
    sel=np.flatnonzero(type_id==t); bins={'interface_le_4A':d_opposite[sel]<=4.0,'interface_4_to_8A':(d_opposite[sel]>4.0)&(d_opposite[sel]<=8.0),'interior_gt_8A':d_opposite[sel]>8.0}
    report['by_element'][name]={'count':int(len(sel)),'nearest_opposite_phase_distance_A':stats(d_opposite[sel]),'force_by_opposite_phase_distance':{label:{'n':int(mask.sum()),'full_force_eV_A':stats(fmag[sel][mask]),'no_cross_force_eV_A':stats(nmag[sel][mask]),'cross_component_eV_A':stats(cmag[sel][mask])} for label,mask in bins.items()},'top_10_full_force_atoms':top_atoms(sel,10),'max_full_force_eV_A':float(np.max(fmag[sel])),'max_no_cross_force_eV_A':float(np.max(nmag[sel])),'max_cross_component_eV_A':float(np.max(cmag[sel]))}
# The dynamics/minimization tests held all atoms below z=20 A fixed and kept the
# entire ceramic skeleton rigid. Report those populations separately from mobile Ag.
report['by_population']={}
for label,sel in {
    'mobile_Ag_above_20A':(type_id==4)&(full[:,4]>=20.0),
    'fixed_bottom_Ag_below_20A':(type_id==4)&(full[:,4]<20.0),
}.items():
    ids=np.flatnonzero(sel)
    bins={'interface_le_4A':d_opposite[ids]<=4.0,'interface_4_to_8A':(d_opposite[ids]>4.0)&(d_opposite[ids]<=8.0),'interior_gt_8A':d_opposite[ids]>8.0}
    report['by_population'][label]={
        'count':int(len(ids)),
        'full_force_eV_A':stats(fmag[ids]),
        'no_cross_force_eV_A':stats(nmag[ids]),
        'cross_component_eV_A':stats(cmag[ids]),
        'force_by_opposite_phase_distance':{k:{'n':int(m.sum()),'full_force_eV_A':stats(fmag[ids][m]),'no_cross_force_eV_A':stats(nmag[ids][m]),'cross_component_eV_A':stats(cmag[ids][m])} for k,m in bins.items()},
        'top_10_full_force_atoms':top_atoms(ids,10),
    }
(HERE/'high_force_localization_after300.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
