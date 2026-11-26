#!/usr/bin/env python3
"""Generate a fresh 1-A voxel QSGS field and screen for thin balanced branches."""
from pathlib import Path
import importlib.util,itertools,json
import numpy as np
from scipy import ndimage
from ase.build import bulk

HERE=Path(__file__).resolve().parent
HELPER=HERE/'support'/'build_microtwig_candidate.py'
spec=importlib.util.spec_from_file_location('qsgs_helpers',HELPER)
g=importlib.util.module_from_spec(spec); spec.loader.exec_module(g)
g.SIZE=np.array([200.,200.,200.]);g.GRID=np.array([200,200,200]);g.VACUUM_PADDING=1000.
SEED=20261141;GRID=g.GRID;SIZE=g.SIZE
P_SEED=.0016;P_FORWARD=.070;P_BRANCH=.030;P_BACKWARD=.00060;TARGET=.84
OFF=np.asarray([p for p in itertools.product((-1,0,1),repeat=3) if any(p)],np.int8)
UNIT=OFF.astype(float);UNIT/=np.linalg.norm(UNIT,axis=1)[:,None]
RNG=np.random.default_rng(SEED)

def shift_open_z(a,o):
 b=np.roll(a,tuple(o),axis=(0,1,2));dz=int(o[2])
 if dz>0:b[:,:,:dz]=False
 elif dz<0:b[:,:,dz:]=False
 return b

def grow_rank():
 shape=tuple(GRID);nvox=int(np.prod(GRID));occupied=RNG.random(shape)<P_SEED
 heading=np.full(shape,-1,np.int8);seeds=np.flatnonzero(occupied)
 heading.flat[seeds]=RNG.integers(0,len(OFF),len(seeds),dtype=np.int8)
 rank=np.full(shape,nvox,np.int32);RNG.shuffle(seeds);rank.flat[seeds]=np.arange(len(seeds),dtype=np.int32)
 count=len(seeds);generation=0
 while count<int(TARGET*nvox):
  added=np.zeros(shape,bool);next_heading=np.full(shape,-1,np.int8)
  for k in RNG.permutation(len(OFF)):
   o=OFF[k];src=shift_open_z(occupied,o);sh=np.roll(heading,tuple(o),axis=(0,1,2));dz=int(o[2])
   if dz>0:sh[:,:,:dz]=-1
   elif dz<0:sh[:,:,dz:]=-1
   valid=src&(sh>=0)&~occupied&~added;ids=np.flatnonzero(valid)
   if not len(ids):continue
   cosine=UNIT[sh.ravel()[ids]]@UNIT[k]
   prob=np.where(cosine>.5,P_FORWARD,np.where(cosine>=-.5,P_BRANCH,P_BACKWARD))
   take=ids[RNG.random(len(ids))<prob]
   if len(take):added.flat[take]=True;next_heading.flat[take]=k
  ids=np.flatnonzero(added)
  if not len(ids):raise RuntimeError('QSGS growth stopped before target occupancy')
  RNG.shuffle(ids);rank.flat[ids]=np.arange(count,count+len(ids),dtype=np.int32)
  heading[added]=next_heading[added];occupied|=added;count+=len(ids);generation+=1
  if generation%5==0:print('growth',generation,f'{count/nvox:.1%}',flush=True)
 print('rank complete',generation,f'{count/nvox:.1%}',flush=True)
 return rank,generation

def largest_mask(rank,f):
 labels,n=ndimage.label(rank<int(f*rank.size));counts=np.bincount(labels.ravel());counts[0]=0
 return labels==int(np.argmax(counts))

def main():
 rank,generations=grow_rank();np.save(HERE/'qsgs_rank_1A_seed20261141.npy',rank)
 print('building TSC/Ag crystal-site pools',flush=True)
 tp,ts=g.build_tsc_lattice();aglat=bulk('Ag','fcc',a=g.AG_A0,cubic=True).repeat((49,49,49));aglat.set_cell(np.diag(SIZE),scale_atoms=True);ap=aglat.positions.copy()
 results=[]
 for f in [.48,.50,.52,.54,.56,.58]:
  print('screen threshold',f,flush=True)
  mask=largest_mask(rank,f)
  tsc=tp[g.inside_mask(mask,tp)];syms=ts[g.inside_mask(mask,tp)]
  sk,ss,fu,hist,tiinfo=g.clean_atomic_stoichiometry(tsc,syms)
  candidates=ap[~g.inside_mask(mask,ap)];d,_=g.tree(sk).query(candidates,workers=4)
  candidates=candidates[d>=2.0]
  keep,aginfo=g.connected_component(candidates,3.15);ag=candidates[keep]
  massratio=len(ag)*g.MASS['Ag']/(fu*g.FORMULA_MASS)
  dist=ndimage.distance_transform_edt(mask,sampling=(1,1,1))[mask]
  bins=mask.reshape(10,20,10,20,10,20).mean(axis=(1,3,5))
  top=mask[:,:,180:200].reshape(10,20,10,20,20).mean(axis=(1,3,4))
  row={'rank_threshold':f,'mask_fraction':float(mask.mean()),'TSC_formula_units':int(fu),
       'TSC_atoms':int(len(sk)),'Ag_main_atoms':int(len(ag)),'Ag_TSC_mass_ratio':float(massratio),
       'TSC_components_discarded':hist,'Ag_component_info':aginfo,
       'TSC_branch_inside_distance_A':{'p50':float(np.percentile(dist,50)),'p90':float(np.percentile(dist,90)),
         'p99':float(np.percentile(dist,99)),'max':float(dist.max())},
       '20A_bins':{'p10':float(np.percentile(bins,10)),'median':float(np.median(bins)),
         'p90':float(np.percentile(bins,90)),'max':float(bins.max()),'over90':int(np.count_nonzero(bins>.9))},
       'top_20A_bins':{'mean':float(top.mean()),'max':float(top.max()),'over90':int(np.count_nonzero(top>.9))},
       'TSC_topology':g.winding_rank(sk,3.25),'Ag_topology':g.winding_rank(ag,3.15)}
  results.append(row);print(json.dumps(row),flush=True)
 (HERE/'qsgs_1A_screen_seed41_results.json').write_text(json.dumps({'seed':SEED,'voxel_A':1.0,
   'parameters':{'P_seed':P_SEED,'P_forward':P_FORWARD,'P_branch':P_BRANCH,'P_backward':P_BACKWARD,
     'growth_target':TARGET,'generations':generations},'results':results},indent=2)+'\n')

if __name__=='__main__':main()
