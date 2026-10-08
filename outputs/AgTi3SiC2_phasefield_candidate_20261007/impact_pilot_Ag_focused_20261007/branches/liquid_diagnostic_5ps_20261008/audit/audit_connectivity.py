from pathlib import Path
import sys
import numpy as np
from scipy.spatial import cKDTree
from audit_crater import read_data


def largest_component(coords, ids, cutoff, a, b):
    n=len(coords)
    points=np.concatenate([coords+i*a+j*b for i in (-1,0,1) for j in (-1,0,1)],axis=0)
    labels=np.tile(np.arange(n),9)
    tree=cKDTree(points)
    parent=np.arange(n)
    size=np.ones(n,dtype=int)
    def find(x):
        while parent[x]!=x:
            parent[x]=parent[parent[x]]
            x=parent[x]
        return x
    def union(x,y):
        rx,ry=find(x),find(y)
        if rx==ry: return
        if size[rx]<size[ry]: rx,ry=ry,rx
        parent[ry]=rx; size[rx]+=size[ry]
    for i,p in enumerate(coords):
        for k in tree.query_ball_point(p,cutoff):
            j=int(labels[k])
            if j!=i: union(i,j)
    counts={}
    for i in range(n):
        r=find(i); counts[r]=counts.get(r,0)+1
    largest=max(counts.values()) if counts else 0
    return largest, len(counts), largest/n if n else 0


def run(path):
    d=read_data(path)
    a=np.array([d['lx'],0.0,0.0]); b=np.array([d['tilt'][0],d['ly'],0.0])
    print(f'[{path}] total={len(d["ids"])}')
    for name,mask,cutoff in [
      ('Ti3SiC2',d['typ']!=4,3.1),
      ('Ag',d['typ']==4,3.35),
    ]:
      largest,ncomp,fraction=largest_component(d['xyz'][mask],d['ids'][mask],cutoff,a,b)
      print(f'  {name}: cutoff={cutoff:.2f} A; largest={largest}/{mask.sum()} ({fraction:.4f}); clusters={ncomp}')

if __name__=='__main__':
    run(sys.argv[1] if len(sys.argv)>1 else 'structures/base_active_local_relaxed.data')
    run(sys.argv[2] if len(sys.argv)>2 else 'structures/cool_04_quench_0p5ps.data')
