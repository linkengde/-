#!/workspace/venvs/dft-labels/bin/python
"""Select Round7 local crops for a locked, never-trained-on DFT blind set."""
from __future__ import annotations
import csv, gzip, json, math
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.io import write
from scipy.spatial import cKDTree

HERE=Path(__file__).resolve().parent
RUN=Path('/workspace/-/models/ag_ti3sic2_qsgs_200x200x150_md_finer_skeleton/runs/round7_fiveframe_joint_analysis_to_1p6ps_20261109')
CORRECTED=Path('/workspace/-/models/ag_ti3sic2_qsgs_200x200x150_md_finer_skeleton/runs/round7_corrected_radial_core_analysis_1p6ps_20261109')
TRAJ=RUN/'round7_5frames_full_system_0p6_to_1p6ps.dump.gz'
CLUSTER=CORRECTED/'round7_1p6ps_largest_candidate_cluster.dump'
CANDIDATE_SUMMARY=CORRECTED/'round7_1p6ps_joint_screen.json'
RADIUS=5.5; CELL=22.0; LX=200.0; LY=200.0
SYM={1:'Ti',2:'Si',3:'C',4:'Ag'}
WANTED={6000,11000,16000}

def read_frames(path):
    out={}
    with gzip.open(path,'rt') as f:
        while True:
            line=f.readline()
            if not line: break
            if not line.startswith('ITEM: TIMESTEP'): raise ValueError(line)
            step=int(f.readline())
            if not f.readline().startswith('ITEM: NUMBER OF ATOMS'): raise ValueError('atom count header')
            n=int(f.readline())
            bh=f.readline().split()
            bounds=np.array([list(map(float,f.readline().split()[:2])) for _ in range(3)])
            ah=f.readline().split()[2:]
            ix={v:i for i,v in enumerate(ah)}
            raw=np.loadtxt((f.readline() for _ in range(n)),dtype=float,ndmin=2)
            if step in WANTED:
                o=np.argsort(raw[:,ix['id']]); raw=raw[o]
                out[step]={'step':step,'ids':raw[:,ix['id']].astype(np.int64),
                           'types':raw[:,ix['type']].astype(np.int8),
                           'xyz':raw[:,[ix['x'],ix['y'],ix['z']]],
                           'bounds':bounds,'box_header':' '.join(bh)}
    missing=WANTED-set(out)
    if missing: raise RuntimeError(f'missing snapshots {sorted(missing)}')
    return out

def read_ids(path):
    lines=path.read_text().splitlines()
    i=next(i for i,x in enumerate(lines) if x.startswith('ITEM: ATOMS'))
    cols=lines[i].split()[2:]; ix={k:j for j,k in enumerate(cols)}
    a=np.loadtxt(lines[i+1:],dtype=float,ndmin=2)
    return a[:,ix['id']].astype(np.int64)

def circular_center(x):
    z=np.exp(2j*np.pi*x/LX).mean()
    return float((np.angle(z)%(2*np.pi))*LX/(2*np.pi))

def pbc_delta(xyz,center):
    d=xyz-center
    d[:,0]-=LX*np.rint(d[:,0]/LX)
    d[:,1]-=LY*np.rint(d[:,1]/LY)
    return d

def choose_nearest_target(frame,target,species=4,exclude_ids=()):
    mask=frame['types']==species
    if exclude_ids: mask &= ~np.isin(frame['ids'],np.asarray(list(exclude_ids),dtype=np.int64))
    inds=np.flatnonzero(mask)
    d=pbc_delta(frame['xyz'][inds],np.asarray(target,dtype=float))
    dist=np.linalg.norm(d,axis=1)
    return int(frame['ids'][inds[np.argmin(dist)]]),float(dist.min())

def crop(frame,center_id,reason,partner_id=None):
    i=int(np.searchsorted(frame['ids'],center_id))
    if frame['ids'][i]!=center_id: raise KeyError(center_id)
    d=pbc_delta(frame['xyz'],frame['xyz'][i]); dist=np.linalg.norm(d,axis=1)
    keep=np.flatnonzero(dist<=RADIUS)
    local=Atoms([SYM[int(t)] for t in frame['types'][keep]],positions=d[keep]+CELL/2,
                cell=np.eye(3)*CELL,pbc=True)
    local.new_array('source_atom_id',frame['ids'][keep].astype(np.int32))
    if partner_id is None:
        others=keep[frame['ids'][keep]!=center_id]
        partner_id=int(frame['ids'][others[np.argmin(dist[others])]])
    pi=int(np.searchsorted(frame['ids'],partner_id))
    local.info.update({'config_type':f'blind_round7_step{frame["step"]}_{reason}',
       'geometry_source':str(TRAJ),'source_step':int(frame['step']),
       'center_atom_id':int(center_id),'center_species':SYM[int(frame['types'][i])],
       'partner_atom_id':int(partner_id),'partner_species':SYM[int(frame['types'][pi])],
       'center_partner_distance_A':float(dist[pi]),'crop_radius_A':RADIUS,'crop_cell_A':CELL,
       'source_pbc':'x/y periodic, z open; crop embedded in a cubic vacuum box',
       'label_role':'locked independent blind-test candidate; do not use in fitting or model selection',
       'label_status':'geometry only; DFT pending',
       'boundary_note':'Finite local crop; evaluate only declared center and partner forces; do not mix raw energy with periodic training energies.'})
    return local

frames=read_frames(TRAJ)
cluster_ids=read_ids(CLUSTER)
end=frames[16000]
idmap={int(v):i for i,v in enumerate(end['ids'])}
idx=np.array([idmap[int(v)] for v in cluster_ids if int(v) in idmap],dtype=int)
cluster_xyz=end['xyz'][idx]
cx=circular_center(cluster_xyz[:,0]); cy=circular_center(cluster_xyz[:,1])
cz=float(np.median(cluster_xyz[:,2]))
center_vec=np.array([cx,cy,cz])
dd=pbc_delta(cluster_xyz,center_vec); core_id=int(cluster_ids[np.argmin(np.linalg.norm(dd,axis=1))])

# Pick an Ag from the largest candidate component that is closest to the ceramic.
skel=np.flatnonzero(end['types']!=4)
skel_xyz=end['xyz'][skel]; skel_ids=end['ids'][skel]
tiles=[]; tile_ids=[]
for sx in (-LX,0,LX):
 for sy in (-LY,0,LY):
  tiles.append(skel_xyz+np.array([sx,sy,0.])); tile_ids.append(skel_ids)
tree=cKDTree(np.concatenate(tiles))
tiles_ids=np.concatenate(tile_ids)
candidate_points=end['xyz'][idx]
dist,near=tree.query(candidate_points,k=1,workers=1)
interface_loc=int(np.argmin(dist)); interface_ag_id=int(cluster_ids[interface_loc])
interface_skel_id=int(tiles_ids[near[interface_loc]])
interface_skel_i=idmap[interface_skel_id]

# Two earlier, spatially separated source-region Ag sites provide a temporal baseline.
early_id,early_dist=choose_nearest_target(frames[6000],[100.,100.,140.],species=4)
mid_id,mid_dist=choose_nearest_target(frames[11000],[125.,100.,140.],species=4,exclude_ids={early_id})
mid_partner=None

selected=[
 ('round7_early_source_ag',6000,early_id,None,'Ag near source-axis position at 0.60 ps'),
 ('round7_mid_source_ag',11000,mid_id,None,'distinct Ag near a displaced source-axis position at 1.10 ps'),
 ('round7_main_pool_core_ag',16000,core_id,None,'Ag nearest the periodic centroid of the corrected 1.60 ps largest candidate component'),
 ('round7_main_pool_skeleton_interface',16000,interface_skel_id,interface_ag_id,'skeleton atom paired to the largest-component Ag nearest the skeleton'),
]
out=[]; manifest=[]
for name,step,cid,pid,reason in selected:
 atoms=crop(frames[step],cid,name,pid)
 out.append(atoms)
 manifest.append({'name':name,'timestep':step,'time_ps':step*0.0001,'center_atom_id':cid,
     'center_species':atoms.info['center_species'],'partner_atom_id':atoms.info['partner_atom_id'],
     'partner_species':atoms.info['partner_species'],'center_partner_distance_A':atoms.info['center_partner_distance_A'],
     'crop_atoms':len(atoms),'formula':atoms.get_chemical_formula(),'reason':reason,
     'training_or_tuning_use':'forbidden','status':'geometry selected; DFT pending'})
write(HERE/'blind_geometries.extxyz',out,format='extxyz')
(HERE/'blind_geometry_manifest.json').write_text(json.dumps({'dataset_role':'locked blind test; never include in v11 fitting or model selection',
 'source_trajectory':str(TRAJ),'source_potential':'Round7 MTP energy-force pilot; configurations are model-generated and DFT-unlabeled before this selection',
 'corrected_core_analysis':str(CANDIDATE_SUMMARY),'candidate_cluster_dump':str(CLUSTER),
 'crop_radius_A':RADIUS,'crop_cell_A':CELL,'frames':manifest},indent=2,ensure_ascii=False)+'\n')
print(json.dumps(manifest,indent=2,ensure_ascii=False))
