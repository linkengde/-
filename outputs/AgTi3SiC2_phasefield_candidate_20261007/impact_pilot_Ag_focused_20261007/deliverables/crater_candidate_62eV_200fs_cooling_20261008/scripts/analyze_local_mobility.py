#!/usr/bin/env python3
"""Compare local Ag and Ti3SiC2 displacements over the packaged no-ion hold."""
from pathlib import Path
import gzip
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import audit_crater as ac
start=ac.read_data(ROOT/'structures/model_after_487th_62eV_impact.data')
r=ac.min_image_xy(start);z=start['xyz'][:,2];typ=start['typ']
sel=(r<4)&(z>=10)&(z<35);ids=start['ids'][sel];types=start['typ'][sel]
frames=[]
with gzip.open(ROOT/'trajectory/cool_after_62eV_0p5ps.dump.gz','rt') as f:
 while True:
  line=f.readline()
  if not line: break
  if line.strip()!='ITEM: TIMESTEP': continue
  step=int(f.readline());f.readline();n=int(f.readline());f.readline()
  for _ in range(3):f.readline()
  cols=f.readline().split()[2:];ix=cols.index('id');pos=[cols.index(c) for c in ('xu','yu','zu')]
  arr=np.array([[float(x) for x in f.readline().split()] for _ in range(n)])
  if step in (0,20000):frames.append((step,arr[:,ix].astype(int),arr[:,pos]))
if len(frames)!=2:raise SystemExit('Expected start and final frames')
(_,id0,p0),(_,id1,p1)=frames
p0={int(i):v for i,v in zip(id0,p0)};p1={int(i):v for i,v in zip(id1,p1)}
for name,mask in [('Ag',types==4),('Ti3SiC2',types!=4)]:
 selected=ids[mask];disp=np.array([np.linalg.norm(p1[int(i)]-p0[int(i)]) for i in selected if int(i) in p1])
 print(f'{name}: n={len(disp)} RMSD={np.sqrt(np.mean(disp**2)):.5f} A median={np.median(disp):.5f} A fraction>1A={np.mean(disp>1):.5f} fraction>2A={np.mean(disp>2):.5f}')
