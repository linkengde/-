#!/usr/bin/env python3
"""Run OVITO PTM on representative frames of the 1024-atom hybrid cycle."""
import json
from pathlib import Path
import numpy as np
from ovito.io import import_file
from ovito.modifiers import PolyhedralTemplateMatchingModifier

ROOT=Path(__file__).resolve().parent
RUN=ROOT/'v14_stage2_p0_ag_dft_diagnostics'/'hybrid_1024_ag_melt_cool_2200K'
DUMP=RUN/'cycle.dump'
OUT=RUN/'ptm_selected_frames.json'
pipe=import_file(str(DUMP))
pipe.modifiers.append(PolyhedralTemplateMatchingModifier())
nframes=pipe.source.num_frames
selected=sorted(set([0,*range(0,nframes,20),nframes-1]))
labels={0:'OTHER',1:'FCC',2:'HCP',3:'BCC',4:'ICO',5:'SC',6:'CUBIC_DIAMOND',7:'HEX_DIAMOND',8:'GRAPHENE'}
rows=[]
for frame in selected:
    data=pipe.compute(frame)
    types=np.asarray(data.particles['Particle Type'])
    st=np.asarray(data.particles['Structure Type'])
    ag=types==1
    counts={labels.get(int(t),str(int(t))):int(np.sum(st[ag]==t)) for t in np.unique(st[ag])}
    total={labels.get(int(t),str(int(t))):int(np.sum(st==t)) for t in np.unique(st)}
    rows.append({'frame':int(frame),'step':int(data.attributes['Timestep']),
                 'time_ps':float(data.attributes['Timestep'])*0.001,
                 'Ag_count':int(ag.sum()),'PTM_Ag_counts':counts,'PTM_all_counts':total})
result={'method':'OVITO PolyhedralTemplateMatchingModifier default settings',
        'note':'The Ag film is only a few atomic layers in this small periodic interface; PTM can label under-coordinated or strained atoms as OTHER. This is a structural indicator, not a stand-alone melting criterion.',
        'frame_count':int(nframes),'selected_frames':rows}
OUT.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
