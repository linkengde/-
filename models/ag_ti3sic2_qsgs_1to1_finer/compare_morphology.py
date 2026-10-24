#!/usr/bin/env python3
"""Compare QSGS phase chord lengths, without claiming mechanical neck sizes."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

ROOT=Path(__file__).resolve().parent
REPO=next((p for p in ROOT.parents if (p/'research/mace-v12-transfer').is_dir()),Path('/workspace/-'))
OLD=REPO/'models/ag_ti3sic2_qsgs_1to1_fine'
def metrics(path):
    data=np.load(path/'qsgs_morphology.npz')
    mask=data['skeleton'];size=data['size_A'];spacing=size/np.array(mask.shape)
    result={'skeleton_volume_fraction':float(mask.mean()),'mean_skeleton_voxel_chord_A':{}}
    for axis,name in enumerate('xyz'):
        if axis<2:
            transitions=np.count_nonzero(mask!=np.roll(mask,1,axis=axis))
        else:
            pad=np.pad(mask,((0,0),(0,0),(1,1)))
            transitions=np.count_nonzero(np.diff(pad.astype(np.int8),axis=2))
        result['mean_skeleton_voxel_chord_A'][name]=float(2*mask.sum()*spacing[axis]/transitions)
    result['mean_of_xyz_chords_A']=float(np.mean(list(result['mean_skeleton_voxel_chord_A'].values())))
    return mask,size,result
old,size,old_result=metrics(OLD)
new,_,new_result=metrics(ROOT)
result={'original':old_result,'finer':new_result,'change_in_mean_chord_percent':100*(new_result['mean_of_xyz_chords_A']/old_result['mean_of_xyz_chords_A']-1),'interpretation':'Directional mean line-intercept chords of the target voxel skeleton, not a measurement of minimum atomistic neck thickness.'}
(ROOT/'morphology_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
fig,axes=plt.subplots(1,2,figsize=(10,5.8))
for ax,mask,title in zip(axes,[old,new],['Previous fine skeleton','Further refined skeleton']):
    ax.imshow(mask[:,:,25].T,origin='lower',interpolation='nearest',extent=[0,size[0],0,size[1]],cmap=ListedColormap(['#bcbfc4','#205dcc']),vmin=0,vmax=1)
    ax.set(title=title,xlabel='x (A)',ylabel='y (A)')
fig.suptitle('Blue: Ti3SiC2 skeleton | Grey: Ag target domain\nSame box and Ag:Ti3SiC2 atomic mass ratio = 1:1')
fig.tight_layout();fig.savefig(ROOT/'morphology_comparison.png',dpi=180)
print(json.dumps(result,indent=2))
