from pathlib import Path
import json,numpy as np
from scipy.ndimage import distance_transform_edt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
ROOT=Path(__file__).resolve().parent
OLD=ROOT.parent/'qsgs_ag60_tsc40_300x300x150_20261021'/'model';NEW=ROOT/'model'
fig,axs=plt.subplots(1,2,figsize=(10,5),dpi=200);reports=[]
for ax,path,title in zip(axs,[OLD,NEW],['V1: coarse skeleton','V2: finer skeleton']):
 a=json.loads((path/'GEOMETRY_AUDIT.json').read_text());field=np.load(path/'qsgs_growth.npz')['arrival'];mask=field<=a['composition_trials'][-1]['threshold']
 ax.imshow(mask[:,:,37].T,origin='lower',extent=[0,300,0,300],cmap=ListedColormap(['#bfc1c3','#2158cf']),interpolation='nearest');ax.set_title(title);ax.set(xlabel='x / A',ylabel='y / A')
 metrics={}
 for phase,v in [('TSC',mask),('Ag',~mask)]:
  padded=np.pad(v,((30,30),(30,30),(0,0)),mode='wrap');padded=np.pad(padded,((0,0),(0,0),(30,30)),mode='edge')
  d=distance_transform_edt(padded,sampling=2.)[30:-30,30:-30,30:-30][v]
  metrics[phase]={'median_distance_to_interface_A':float(np.median(d)),'mean_distance_to_interface_A':float(d.mean()),'p90_distance_to_interface_A':float(np.quantile(d,.9))}
 reports.append({'version':title,'core_probability':a['parameters']['core_probability'],'atoms':a['atoms'],'Ag_mass_fraction':a['Ag_mass_fraction'],'mask_TSC_volume_fraction':float(mask.mean()),'distance_proxy':metrics,'connectivity':a['connectivity']})
fig.suptitle('Same XY plane at z=75 A | blue Ti3SiC2, gray Ag (QSGS masks)');fig.tight_layout();fig.savefig(NEW/'V1_V2_COMPARISON.png')
result={'versions':reports,'metric_note':'Voxel distance to nearest phase interface, not local strut diameter or relaxed atomic thickness. XY periodic padding and constant z extension used; 2 A grid resolution.'}
(NEW/'REFINEMENT_COMPARISON.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
