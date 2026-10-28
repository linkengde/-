#!/usr/bin/env python3
"""Audit composition, periodic winding, open-z spanning, and voxel-scale features."""
from __future__ import annotations
import hashlib, json
from pathlib import Path
import numpy as np
from scipy import ndimage

HERE=Path(__file__).resolve().parent
GEOMETRY=HERE/'nano_geometry_ag607wt_ti3sic2_393wt.npz'

def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()

def periodic_labels(mask, periodic_axes):
 labels,count=ndimage.label(mask)
 parent=np.arange(count+1,dtype=np.int64)
 def root(a):
  a=int(a)
  while parent[a]!=a:
   parent[a]=parent[parent[a]]
   a=int(parent[a])
  return a
 for axis in periodic_axes:
  low=np.take(labels,0,axis=axis).ravel()
  high=np.take(labels,-1,axis=axis).ravel()
  pairs=np.unique(np.stack((low,high),axis=1),axis=0)
  for a,b in pairs:
   if a and b:
    ra,rb=root(a),root(b)
    if ra!=rb: parent[rb]=ra
 roots=np.fromiter((root(i) for i in range(count+1)),dtype=np.int64,count=count+1)
 return roots[labels]

def spans_opposite_faces(mask, axis, other_periodic_axes):
 labels=periodic_labels(mask,other_periodic_axes)
 low=np.unique(np.take(labels,0,axis=axis))
 high=np.unique(np.take(labels,-1,axis=axis))
 return bool(np.intersect1d(low[low>0],high[high>0],assume_unique=True).size)

def chords(mask,axis,periodic,spacing):
 rows=np.moveaxis(mask,axis,-1).reshape(-1,mask.shape[axis])
 sizes=[]
 for row in rows:
  if not row.any(): continue
  if row.all(): sizes.append(row.size); continue
  if periodic:
   start=int(np.flatnonzero(row!=np.roll(row,1))[0])
   row=np.roll(row,-start)
  pad=np.concatenate(([False],row,[False]))
  trans=np.flatnonzero(pad[1:]!=pad[:-1])
  sizes.extend((trans[1::2]-trans[::2]).tolist())
 return np.asarray(sizes,dtype=float)*spacing

def describe(values):
 if not values.size: return {'count':0}
 return {'count':int(values.size),'p10_A':float(np.percentile(values,10)),
         'median_A':float(np.median(values)),'p90_A':float(np.percentile(values,90)),
         'max_A':float(values.max())}

with np.load(GEOMETRY,allow_pickle=False) as z:
 skeleton=np.asarray(z['skeleton'],dtype=bool)
 size=np.asarray(z['size_A'],dtype=float)
 parent_hash=str(z['parent_geometry_sha256'].item())
 rank_hash=str(z['growth_rank_sha256'].item())
spacing=size/np.asarray(skeleton.shape)
phases={'Ag':~skeleton,'Ti3SiC2':skeleton}
metrics={}
for name,mask in phases.items():
 lbl=periodic_labels(mask,(0,1))
 count=int(np.unique(lbl[lbl>0]).size)
 zlow=np.unique(lbl[:,:,0]); zhigh=np.unique(lbl[:,:,-1])
 zspan=bool(np.intersect1d(zlow[zlow>0],zhigh[zhigh>0],assume_unique=True).size)
 tiled=np.tile(mask,(3,3,1))
 dx=ndimage.distance_transform_edt(tiled,sampling=spacing)[mask.shape[0]:2*mask.shape[0],mask.shape[1]:2*mask.shape[1],:]
 diameter=2.0*dx[mask]
 metrics[name]={
  'voxel_fraction':float(mask.mean()),
  'periodic_xy_open_z_components':count,
  'largest_component_fraction':float(np.bincount(lbl.ravel())[1:].max()/mask.sum()),
  'periodic_x_winding':spans_opposite_faces(mask,0,(1,)),
  'periodic_y_winding':spans_opposite_faces(mask,1,(0,)),
  'open_z_spans_both_surfaces':zspan,
  'z_surface_touch_fraction':{'bottom':float(mask[:,:,0].mean()),'top':float(mask[:,:,-1].mean())},
  'interface_distance_diameter_A':{
   'minimum':float(diameter.min()),'p01':float(np.percentile(diameter,1)),
   'p05':float(np.percentile(diameter,5)),'median':float(np.median(diameter)),
   'p95':float(np.percentile(diameter,95)),'maximum':float(diameter.max()),
   'fraction_below_6A':float(np.mean(diameter<6.0)),
   'fraction_below_10A':float(np.mean(diameter<10.0)),
   'method':'2x voxel-center Euclidean distance to non-phase; x/y are periodic via 3x tiling and z is open'
  },
  'lineal_chord_A':{a:describe(chords(mask,axis,axis<2,float(spacing[axis]))) for axis,a in enumerate(('x','y','z'))}
 }
expected={'Ag':0.3995857976962546,'Ti3SiC2':0.6004142023037454}
for name in phases:
 m=metrics[name]
 if abs(m['voxel_fraction']-expected[name])>2e-6:
  raise AssertionError(f'{name} fraction outside the 1-voxel discretization: {m["voxel_fraction"]}')
 if m['periodic_xy_open_z_components']!=1 or m['largest_component_fraction']<1.0:
  raise AssertionError(f'{name} phase is not one component under periodic x/y and open z')
 if not (m['periodic_x_winding'] and m['periodic_y_winding'] and m['open_z_spans_both_surfaces']):
  raise AssertionError(f'{name} phase does not percolate in all required directions')
report={
 'status':'AG607WT_GEOMETRY_AUDIT_PASSED',
 'geometry':str(GEOMETRY.relative_to(HERE.parent)),
 'geometry_sha256':digest(GEOMETRY),
 'shape_xyz':list(skeleton.shape),
 'domain_size_A':size.tolist(),
 'voxel_spacing_A':spacing.tolist(),
 'boundary_convention':'periodic x/y, open z',
 'phase_metrics':metrics,
 'source_lineage':{'parent_qsgs_geometry_sha256':parent_hash,'growth_rank_sha256':rank_hash},
 'integrity':{'binary_partition_exact':bool(np.all((~skeleton)^skeleton)),
              'voxel_count':int(skeleton.size),'phase_counts':{'Ag':int((~skeleton).sum()),'Ti3SiC2':int(skeleton.sum())}},
 'interpretation':'Interface-distance diameters and chord lengths are raster-resolved morphology descriptors, not atomistic minimum necks. Features below the ~3 A voxel size are unresolved.'
}
(HERE/'geometry_audit.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
