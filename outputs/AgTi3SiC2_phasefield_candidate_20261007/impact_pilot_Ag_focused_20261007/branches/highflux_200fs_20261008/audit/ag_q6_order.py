#!/usr/bin/env python3
"""Exploratory Ag-only l=6 bond-orientational order audit."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scipy.special import sph_harm_y
import audit_crater as ac

ROOT=Path(__file__).resolve().parents[1]
CENTER=np.array([9.9573898076,24.2578814220])
CUTOFF=3.5
MATRIX_ID_MAX=4305

def q6_by_id(path):
    a=ac.read_data(path)
    r=ac.min_image_xy(a,CENTER)
    keep=(a['typ']==4)&(a['ids']<=MATRIX_ID_MAX)
    local=keep&(r<8.0)&(a['xyz'][:,2]>=12)&(a['xyz'][:,2]<55)
    matrix=np.flatnonzero(keep)
    ids=a['ids'][matrix]
    xyz=a['xyz'][matrix]
    avec=np.array([a['lx'],0.0,0.0]); bvec=np.array([a['tilt'][0],a['ly'],0.0])
    points=np.concatenate([xyz+i*avec+j*bvec for i in (-1,0,1) for j in (-1,0,1)],axis=0)
    labels=np.tile(ids,9)
    tree=cKDTree(points)
    selected=np.flatnonzero(local)
    ids_out=[]; qvals=[]; nneigh=[]
    for ix in selected:
        center= a['xyz'][ix]
        js=tree.query_ball_point(center,CUTOFF)
        vectors=points[js]-center
        nids=labels[js]
        d=np.linalg.norm(vectors,axis=1)
        use=(nids!=a['ids'][ix])&(d>1e-7)&(d<=CUTOFF)
        v=vectors[use]; n=len(v)
        if n<3: continue
        rr=np.linalg.norm(v,axis=1)
        theta=np.arccos(np.clip(v[:,2]/rr,-1,1))
        phi=np.mod(np.arctan2(v[:,1],v[:,0]),2*np.pi)
        vals=np.array([np.mean(sph_harm_y(6,m,theta,phi)) for m in range(-6,7)])
        q=float(np.sqrt(4*np.pi/13*np.sum(np.abs(vals)**2)))
        ids_out.append(int(a['ids'][ix])); qvals.append(q); nneigh.append(n)
    return {'path':str(path),'local_Ag_matrix':int(local.sum()),'valid_centers':len(qvals),
            'ids':ids_out,'q6':qvals,'nneigh':nneigh}

def summarize(v):
    q=np.asarray(v['q6']); n=np.asarray(v['nneigh'])
    return {'local_Ag_matrix':v['local_Ag_matrix'],'valid_centers':len(q),
            'Ag_Ag_neighbors_mean':float(n.mean()) if len(n) else None,
            'Q6_mean':float(q.mean()) if len(q) else None,
            'Q6_median':float(np.median(q)) if len(q) else None,
            'Q6_p25':float(np.percentile(q,25)) if len(q) else None,
            'Q6_p75':float(np.percentile(q,75)) if len(q) else None,
            'fraction_Q6_gt_0p35':float(np.mean(q>0.35)) if len(q) else None,
            'fraction_with_10plus_Ag_neighbors':float(np.mean(n>=10)) if len(n) else None}

def common_change(a,b):
    qa=dict(zip(a['ids'],a['q6'])); qb=dict(zip(b['ids'],b['q6']))
    ids=sorted(set(qa)&set(qb)); delta=np.array([qb[i]-qa[i] for i in ids])
    return {'common_matrix_Ag_ids':len(ids),
            'mean_Q6_change':float(delta.mean()) if len(delta) else None,
            'median_Q6_change':float(np.median(delta)) if len(delta) else None,
            'fraction_Q6_drop_gt_0p1':float(np.mean(delta<-.1)) if len(delta) else None}

files={
 'pristine_base':'structures/pristine_base_4305.data',
 'after40_cool':'structures/start_after40.data',
 'after45_cool':'structures/cooled_after45_0p5ps.data'}
raw={k:q6_by_id(ROOT/v) for k,v in files.items()}
summary={k:summarize(v) for k,v in raw.items()}
summary['common_id_changes']={'pristine_to_after45':common_change(raw['pristine_base'],raw['after45_cool']),
                              'after40_to_after45':common_change(raw['after40_cool'],raw['after45_cool'])}
summary['method']='Ag-only Steinhardt Q6 from Ag-Ag neighbors within 3.5 A; local shell r<8 A, 12<=z<55 A. Exploratory ordering indicator; threshold 0.35 is descriptive, not a liquid/solid classifier.'
(ROOT/'audit/ag_q6_order.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
