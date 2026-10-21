"""Seeded QSGS geometry only. No MD, no inherited trajectory or checkpoint."""
from pathlib import Path
import os,json,hashlib,time,csv
os.environ.setdefault('OMP_NUM_THREADS','1')
import numpy as np
from scipy.ndimage import convolve, label
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from ase import Atoms
from ase.io import read,write
from ase.build import make_supercell,bulk

ROOT=Path(__file__).resolve().parent; OUT=ROOT/'model';OUT.mkdir(exist_ok=True)
CIF=ROOT/'source/Ti3SiC2_COD_9009647.cif'
SEED=20261021;shape=(150,150,75);h=2.; rng=np.random.default_rng(SEED)
parameters={'seed':SEED,'voxel_A':h,'grid_shape':shape,'core_probability':0.0002,'growth_probability_face':0.18,'direction_weights':{'face':1.,'edge':2/3,'corner':4/9},'Ag_mass_fraction_target':0.6,'target_skeleton_volume_fraction':'calibrated to retained atomic masses after Ag interface exclusion','method':'QSGS seeded stochastic directional growth; no Gaussian random field','periodicity':[True,True,False]}
print('Growing QSGS',flush=True)
occupied=rng.random(shape)<parameters['core_probability'];arrival=np.full(shape,np.inf,np.float32)
arrival[occupied]=rng.random(occupied.sum()).astype(np.float32)
kernel=np.zeros((3,3,3),np.float32)
for ix in range(3):
 for iy in range(3):
  for iz in range(3):
   n=abs(ix-1)+abs(iy-1)+abs(iz-1)
   if n:kernel[ix,iy,iz]=-np.log1p(-0.18*{1:1.,2:2/3,3:4/9}[n])
for iteration in range(1,401):
 # wrap only x/y. Zero-pad z to prevent growth across free surfaces.
 padded=np.pad(occupied,((0,0),(0,0),(1,1)),mode='constant')
 exposure=convolve(padded.astype(np.float32),kernel,mode='wrap')[:,:,1:-1]
 chance=-np.expm1(-exposure)
 new=(~occupied)&(rng.random(shape)<chance)
 arrival[new]=iteration+rng.random(new.sum()).astype(np.float32);occupied|=new
 if iteration%10==0:print('growth',iteration,occupied.mean(),flush=True)
 if occupied.mean()>.97:break
assert occupied.mean()>.97
np.savez_compressed(OUT/'qsgs_growth.npz',arrival=arrival,parameters=json.dumps(parameters))

# c axis of Ti3SiC2 is normal to surface; basal orthorhombic lattice in xy.
t0=make_supercell(read(CIF),np.array([[1,0,0],[1,2,0],[0,0,1]]))
assert np.max(np.abs(t0.cell.array-np.diag(t0.cell.lengths())))<1e-8
t0.wrap();a,b,c=t0.cell.lengths();nt=(round(300/a),round(300/b),int(np.ceil(150/c)))
t=t0.repeat(nt);orig_t=np.array(t.cell.lengths());t.positions[:,0]*=300/orig_t[0];t.positions[:,1]*=300/orig_t[1]
t=t[t.positions[:,2]<150-1e-7];t.set_cell([300,300,250]);t.pbc=[1,1,0]
ag0=bulk('Ag','fcc',a=4.09,cubic=True);na=(round(300/4.09),round(300/4.09),int(np.ceil(150/4.09)))
ag=ag0.repeat(na);ag.positions[:,:2]*=300/(na[0]*4.09);ag=ag[ag.positions[:,2]<150-1e-7]
ag.set_cell([300,300,250]);ag.pbc=[1,1,0]
print('Parent sites',len(t),len(ag),flush=True)
# Keep the same morphology growth order for every composition trial.
def sample(pos):
 idx=np.minimum((pos/h).astype(int),np.array(shape)-1)
 return arrival[tuple(idx.T)]
ts=sample(t.positions);ags=sample(ag.positions)
# Large z period implements xy periodicity with no false z neighbors.
def tree_positions(pos):
 q=pos.copy();q[:,:2]%=300.;q[:,2]+=100.;return q
box=np.array([300.,300.,100000.])
tm=t.get_masses();am=107.8682;trials=[]
last=None
for step in range(17):
 if step==0:lo=float(np.quantile(arrival,.3));hi=float(np.quantile(arrival,.8))
 threshold=(lo+hi)/2
 keep_t=ts<=threshold;keep_a=ags>threshold
 tree=cKDTree(tree_positions(t.positions[keep_t]),boxsize=box)
 distance,_=tree.query(tree_positions(ag.positions[keep_a]),k=1,workers=2)
 selected=np.flatnonzero(keep_a)[distance>=2.2]
 mass_t=float(tm[keep_t].sum());mass_a=len(selected)*am;wt=mass_a/(mass_t+mass_a)
 trials.append({'threshold':threshold,'Ag_wt':wt,'TSC_atoms':int(keep_t.sum()),'Ag_atoms':len(selected),'Ag_interface_sites_removed':int(keep_a.sum()-len(selected))})
 print('composition',step,wt,keep_t.sum(),len(selected),flush=True)
 last=(keep_t,selected,threshold)
 if abs(wt-.6)<0.0001:break
 if wt>.6:lo=threshold
 else:hi=threshold
keep_t,selected,threshold=last
atoms=t[keep_t]+ag[selected];atoms.set_cell([300,300,250]);atoms.pbc=[1,1,0]
phase=np.r_[np.zeros(keep_t.sum(),np.int8),np.ones(len(selected),np.int8)]
atoms.arrays['phase_id']=phase+1
atoms.arrays['atom_id']=np.arange(1,len(atoms)+1)
z=atoms.positions[:,2]
# 1=fixed; 2=bath; 3=Newton. Material membership frozen at construction.
atoms.arrays['layer_id']=np.where(z<10,1,np.where(z<25,2,3)).astype(np.int8)
atoms.arrays['initial_top10A']=((z>=140)&(z<=150)).astype(np.int8)

print('Auditing atom connectivity',flush=True)
audits=[]
for ph,name,cutoff in [(0,'Ti3SiC2',2.75),(1,'Ag',3.2)]:
 ids=np.flatnonzero(phase==ph);pos=tree_positions(atoms.positions[ids]);tree=cKDTree(pos,boxsize=box)
 pairs=tree.query_pairs(cutoff,output_type='ndarray')
 rows=np.r_[pairs[:,0],pairs[:,1]];cols=np.r_[pairs[:,1],pairs[:,0]]
 graph=coo_matrix((np.ones(len(rows),np.int8),(rows,cols)),shape=(len(ids),len(ids))).tocsr()
 ncomp,labs=connected_components(graph,directed=False);sizes=np.bincount(labs);largest=int(sizes.argmax());main=pos[labs==largest]
 # Reach on each opposite physical face, plus periodic-crossing edges in xy.
 inmain=labs[pairs[:,0]]==largest
 dif=pos[pairs[inmain,1]]-pos[pairs[inmain,0]]
 crosses=np.any(np.abs(dif[:,:2])>150,axis=0)
 # Cross-boundary edges alone are not winding proof. Lift connected graph,
 # recording cycle residuals to prove periodic winding of the largest component.
 from scipy.sparse.csgraph import breadth_first_order
 order,pred=breadth_first_order(graph,int(np.flatnonzero(labs==largest)[0]),directed=False,return_predecessors=True)
 images=np.zeros((len(ids),2),np.int16)
 for node in order[1:]:
  parent=pred[node];shift=-np.rint((pos[node,:2]-pos[parent,:2])/300).astype(np.int16)
  images[node]=images[parent]+shift
 pp=pairs[inmain];shifts=-np.rint((pos[pp[:,1],:2]-pos[pp[:,0],:2])/300).astype(np.int16)
 winding=np.any(images[pp[:,1]]-images[pp[:,0]]!=shifts,axis=0)
 actual_z=main[:,2]-100
 info={'phase':name,'cutoff_A':cutoff,'atoms':len(ids),'components':int(ncomp),'largest_fraction':float(sizes[largest]/len(ids)),'largest_atoms':int(sizes[largest]),'minor_component_atoms':int(len(ids)-sizes[largest]),'periodic_winding_xy':winding.tolist(),'z_extent_A':[float(actual_z.min()),float(actual_z.max())],'spans_z':bool(actual_z.min()<5 and actual_z.max()>145)}
 audits.append(info);print(info,flush=True)
 del graph,pairs,tree,rows,cols
write(OUT/'AG_TSC_300x300x150.extxyz',atoms,format='extxyz')
write(OUT/'AG_TSC_300x300x150.data',atoms,format='lammps-data',atom_style='atomic',specorder=['Ti','Si','C','Ag'],masses=True)
with (OUT/'ATOM_GROUPS.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['atom_id','element','phase','layer','initial_z_A'])
 for i,(sym,ph,layer,zz) in enumerate(zip(atoms.get_chemical_symbols(),phase,atoms.arrays['layer_id'],z),1):w.writerow([i,sym,'Ag' if ph else 'Ti3SiC2',{1:'fixed',2:'thermostat',3:'newton'}[int(layer)],f'{zz:.8f}'])
counts={s:atoms.get_chemical_symbols().count(s) for s in ['Ti','Si','C','Ag']}
report={'status':'GEOMETRY_ONLY_NOT_RELAXED_OR_HEATING_READY','parameters':parameters,'atoms':len(atoms),'counts':counts,'Ag_mass_fraction':trials[-1]['Ag_wt'],'composition_trials':trials,'cell_A':[300,300,250],'solid_bounds_A':[[0,300],[0,300],[0,150]],'pbc':[True,True,False],'TSC_in_plane_strain':[300/orig_t[0]-1,300/orig_t[1]-1],'Ag_in_plane_strain':300/(na[0]*4.09)-1,'interface_rule':'remove Ag candidate sites within 2.20 A of any retained TSC site; TSC sites unchanged','minimum_cross_distance_required_A':2.2,'connectivity':audits,'layers':{'fixed':[0,10],'thermostat':[10,25],'newton':[25,150],'vacuum':[150,250]},'layer_counts':{str(i):int((atoms.arrays['layer_id']==i).sum()) for i in [1,2,3]},'potential':'calibration-v0 (not evaluated)','CIF_sha256':hashlib.sha256(CIF.read_bytes()).hexdigest(),'TTM':'not started; no material map generated','note':'QSGS parameters are newly chosen modeling assumptions, not recovered historical parameters. Minor components are reported, not silently removed.'}
(OUT/'GEOMETRY_AUDIT.json').write_text(json.dumps(report,indent=2)+'\n')
# voxel cross sections are a more legible morphology preview than a million-point projection.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
mask=arrival<=threshold
fig,axs=plt.subplots(1,3,figsize=(12,4),dpi=180)
for ax,im,title,extent in [(axs[0],mask[:,:,37].T,'XY, z=75 A',(0,300,0,300)),(axs[1],mask[:,75,:].T,'XZ, y=150 A',(0,300,0,150)),(axs[2],mask[75,:,:].T,'YZ, x=150 A',(0,300,0,150))]:
 ax.imshow(im,origin='lower',extent=extent,cmap=ListedColormap(['#bfc1c3','#2158cf']),interpolation='nearest');ax.set_title(title);ax.set_aspect('equal')
fig.suptitle('QSGS: blue Ti3SiC2 skeleton; gray Ag pore filling (geometry only)');fig.tight_layout();fig.savefig(OUT/'QSGS_SECTIONS.png');plt.close(fig)
print('DONE',len(atoms),flush=True)
