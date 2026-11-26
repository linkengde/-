#!/usr/bin/env python3
"""Summarize Q6 and OVITO PTM over the seed41 Ag heating trajectory."""
import csv, json
from pathlib import Path
import numpy as np
from ovito.io import import_file
from ovito.modifiers import PolyhedralTemplateMatchingModifier

root=Path(__file__).resolve().parent
traj=root/'seed41_heat_q6_msd_ag.dump'
pipeline=import_file(str(traj))
pipeline.modifiers.append(PolyhedralTemplateMatchingModifier())
summary=[]
for frame in range(pipeline.source.num_frames):
    data=pipeline.compute(frame)
    attrs=data.attributes
    step=int(attrs['Timestep'])
    # The first saved frame is the end of the 0.1 ps cold equilibration;
    # subsequent frames use the heating timestep of 0.0005 ps.
    time_ps=0.1+(step-1000)*0.0005
    n=len(data.particles['Particle Identifier'])
    row={'frame':frame,'step':step,'time_ps':time_ps,'n_Ag':n}
    for key in ('FCC','HCP','BCC','ICO','SC','OTHER'):
        row[key.lower()]=int(attrs.get(f'PolyhedralTemplateMatching.counts.{key}',0))
    row['ordered_fcc_hcp_pct']=100*(row['fcc']+row['hcp'])/n
    q6=np.asarray(data.particles['c_q6ag[1]'])
    row['q6_mean']=float(np.mean(q6))
    row['q6_gt_0p30_pct']=100*float(np.mean(q6>0.30))
    summary.append(row)
    print(json.dumps(row),flush=True)

(root/'seed41_ptm_q6_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
with (root/'seed41_ptm_q6_summary.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(summary[0]))
    w.writeheader();w.writerows(summary)
