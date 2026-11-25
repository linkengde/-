#!/usr/bin/env python3
"""Count OVITO PTM structure classes in each Ag-only LAMMPS dump frame."""
import csv, sys
import numpy as np
from ovito.io import import_file
from ovito.modifiers import PolyhedralTemplateMatchingModifier

if len(sys.argv) != 3:
    raise SystemExit('usage: analyze_ptm.py trajectory.dump output.csv')
source, out = sys.argv[1:]
pipeline = import_file(source, multiple_frames=True)
ptm = PolyhedralTemplateMatchingModifier()
ptm.output_rmsd = True
pipeline.modifiers.append(ptm)
rows=[]
seen=set()
for frame in range(pipeline.source.num_frames):
    data=pipeline.compute(frame)
    step=int(data.attributes.get('Timestep', frame))
    if step in seen:
        continue
    seen.add(step)
    st=np.asarray(data.particles['Structure Type'],dtype=np.int32)
    counts=np.bincount(st,minlength=9)
    n=len(st)
    rows.append({
        'step':step,'time_ps':step*0.001,'n_Ag':n,
        'PTM_other':int(counts[0]),'PTM_FCC':int(counts[1]),
        'PTM_HCP':int(counts[2]),'PTM_BCC':int(counts[3]),
        'PTM_ICO':int(counts[4]),'PTM_FCC_HCP_fraction':float((counts[1]+counts[2])/n),
        'PTM_other_fraction':float(counts[0]/n),
    })
with open(out,'w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=rows[0].keys())
    writer.writeheader(); writer.writerows(rows)
print(f'Wrote {len(rows)} unique frames to {out}')
