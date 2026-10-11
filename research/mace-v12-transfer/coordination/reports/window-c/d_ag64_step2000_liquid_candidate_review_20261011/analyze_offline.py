"""Read archived coordinates/logs only. No potential evaluation or MD."""
from pathlib import Path
import csv,io,json,hashlib,subprocess,tarfile,tempfile,collections
import numpy as np
from scipy.special import sph_harm_y
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ase import units
from ase.io import read
from ase.build import bulk
OUT=Path(__file__).resolve().parent
OLD_COMMIT='145ee8050cf95a3a1988e868c6b935e10538f723'
COMMIT='bc3ff0631db3b7d78dcc17556befa98c2f454480'
OLD_PATH='codex_warmup_archive/C-D-AGL-1250K-rho097-seed01/b62a818e58a2cf42/'
PATH='codex_segment_archives/D/C-D-AGL-1250K-rho097-seed01/step-000002000-c3a88a759e0bd258e74632aba0a68e4f1a765fcbe8324034d5bf127b5aa69965.tar.gz'
SHA=lambda b:hashlib.sha256(b).hexdigest()
def blob(commit,p):return subprocess.check_output(['git','show',commit+':'+p])
def members(b):
 with tarfile.open(fileobj=io.BytesIO(b)) as t:
  return {m.name:t.extractfile(m).read() for m in t.getmembers() if m.isfile()}
def csvrows(b):
 rows=list(csv.DictReader(io.StringIO(b.decode())))
 return [{k:(v if k=='phase' else float(v)) for k,v in row.items()} for row in rows]
def savecsv(name,rows):
 with (OUT/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for row in rows for k in row)));w.writeheader();w.writerows(rows)
def vectors(pos,cell):
 d=pos[None,:,:]-pos[:,None,:];s=d@np.linalg.inv(cell);s-=np.rint(s)
 return s@cell

def q6(pos,cell,cut):
 v=vectors(pos,cell);d=np.linalg.norm(v,axis=-1);qs=[];bondq=[]
 for i in range(len(pos)):
  vv=v[i,(d[i]>0)&(d[i]<cut)];rr=np.linalg.norm(vv,axis=1)
  theta=np.arccos(np.clip(vv[:,2]/rr,-1,1));phi=np.arctan2(vv[:,1],vv[:,0])
  y=np.array([sph_harm_y(6,m,theta,phi) for m in range(-6,7)])
  qq=y.mean(axis=1);qs.append(np.sqrt(4*np.pi/13*np.sum(abs(qq)**2)));bondq.append(y)
 globalq=np.concatenate(bondq,axis=1).mean(axis=1)
 return np.array(qs),float(np.sqrt(4*np.pi/13*np.sum(abs(globalq)**2)))

def longest_chain(adj):
 best=0
 def walk(start,u,seen,length):
  nonlocal best
  best=max(best,length)
  for v in np.flatnonzero(adj[u]):
   if v==start and length>=2:best=max(best,length+1)
   elif v not in seen:walk(start,v,seen|{v},length+1)
 for start in range(len(adj)):walk(start,start,{start},0)
 return best

def cna(pos,cell,cut):
 # Conventional fixed-cutoff CNA; local images are anchored to center i.
 v=vectors(pos,cell);d=np.linalg.norm(v,axis=-1);labels=[]
 for i in range(len(pos)):
  neigh=np.flatnonzero((d[i]>0)&(d[i]<cut));sig=[]
  for j in neigh:
   common=[k for k in neigh if k!=j and np.linalg.norm(v[i,k]-v[i,j])<cut]
   if not common:sig.append((0,0,0));continue
   xyz=v[i,common];dd=np.linalg.norm(xyz[:,None]-xyz[None,:],axis=-1)
   adj=(dd>0)&(dd<cut);sig.append((len(common),int(adj.sum()//2),longest_chain(adj)))
  count=collections.Counter(sig)
  if len(neigh)==12 and count[(4,2,1)]==12:label='FCC'
  elif len(neigh)==12 and count[(4,2,1)]==6 and count[(4,2,2)]==6:label='HCP'
  elif len(neigh)==12 and count[(5,5,5)]==12:label='ICOSAHEDRAL'
  else:label='OTHER_UNCLASSIFIED'
  labels.append(label)
 return labels

def stats(rows):
 t=np.array([r['time_fs'] for r in rows]);result={'n_logged_samples':len(rows),'start_step':int(rows[0]['step']),'end_step':int(rows[-1]['step'])}
 for key in ['T_K','sampler_energy_eV','kinetic_energy_eV','pressure_eV_A3']:
  x=np.array([r[key] for r in rows]);res={'mean':float(x.mean()),'std_ddof1':float(x.std(ddof=1)),'min':float(x.min()),'max':float(x.max()),'ols_slope_per_100fs':float(np.polyfit(t-t[0],x,1)[0]*100)}
  y=x-x.mean();corr=np.correlate(y,y,'full')[len(y)-1:]/np.dot(y,y)
  stop=next((i for i in range(1,len(corr)) if corr[i]<=0),len(corr))
  tau=.5+float(corr[1:stop].sum());res['descriptive_tau_int_fs']=tau*10;res['descriptive_N_eff']=min(len(x),len(x)/(2*tau))
  res['tau_int_status']='SHORT_WINDOW_UNCERTAINTY_NOT_ESTABLISHED';result[key]=res
 return result

oldb=blob(OLD_COMMIT,OLD_PATH+'warmup-b62a818e58a2cf42.tar.gz');oldm=json.loads(blob(OLD_COMMIT,OLD_PATH+'archive-manifest.json'))
assert SHA(oldb)==oldm['archive_sha256']=='b62a818e58a2cf423e4cb1b9ed0468ae73ad578a75c4b2cb9fb3f28baf75e7b1'
old=members(oldb)
for m in oldm['members']:assert SHA(old[m['path']])==m['sha256'] and len(old[m['path']])==m['size_bytes']
b=blob(COMMIT,PATH);assert SHA(b)=='c3a88a759e0bd258e74632aba0a68e4f1a765fcbe8324034d5bf127b5aa69965'
new=members(b);prefix=next(n.rsplit('/',1)[0]+'/' for n in new if n.endswith('/MANIFEST.json'))
manifest=json.loads(new[prefix+'MANIFEST.json'])
for line in new[prefix+'SHA256SUMS.txt'].decode().splitlines():
 h,n=line.split('  ',1);assert SHA(new[prefix+n])==h
for n,h in manifest['members_sha256'].items():assert SHA(new[prefix+n])==h
state=json.loads(new[prefix+'state.json']);assert SHA(new[prefix+'state.json'])=='5d03864f07a98f177c3c0cbdcb850158ed89a64fa2a6ab7f879e8c7ab003be2a'
assert state['step']==2000 and state['task_id']=='C-D-AGL-1250K-rho097-seed01'
assert state['execution_code_sha256']=='3033c3234c4074358f8ef54a6783c3c4c9495549c12c93641f4f777f3451ed76'
assert state['source_provenance']['potential_sha256']=='958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447'
assert state['source_provenance']['execution_release_config_sha256']=='70ca4190ecda29bb8212ba5ab144c36c6c99285b1ef9307393ef19dd409ed52c'
source_manifest=json.loads(old['source_manifest.json'])
oldrows=csvrows(old['thermal.csv']); newrows=csvrows(new[prefix+'thermal.csv'])
for row in oldrows:
 row['kinetic_energy_eV']=.5*189*units.kB*row['T_K'];row['kinetic_energy_source']='DERIVED_FROM_LOGGED_T_WITH_DOF189'
for row in newrows:row['kinetic_energy_source']='LOGGED'
rows=oldrows+newrows;assert len(rows)==101
assert [int(r['step']) for r in rows]==list(range(0,2001,20))
with tempfile.NamedTemporaryFile(suffix='.traj') as f:
 f.write(old['trajectory.traj']);f.flush();atoms=read(f.name,index=':')
assert len(atoms)==11
assert all(a.get_number_of_degrees_of_freedom()==189 for a in atoms)
z=np.load(io.BytesIO(new[prefix+'trajectory.npz']),allow_pickle=False)
assert np.array_equal(atoms[-1].positions,z['positions'][0]);assert np.array_equal(z['steps'],np.arange(1000,2001,100))
pos=np.concatenate([np.array([a.positions for a in atoms]),z['positions'][1:]])
steps=np.arange(0,2001,100);cell=z['cell'];L=cell[0,0];N=64;rho=N/abs(np.linalg.det(cell))
assert all(len(a)==64 and np.all(a.numbers==47) and np.array_equal(a.cell.array,cell) and np.all(a.pbc) for a in atoms)
assert np.array_equal(cell,np.eye(3)*L);assert np.array_equal(z['atom_ids'],np.arange(1,65));assert np.all(z['pbc'])
checkpoint=np.load(io.BytesIO(new[prefix+'state_arrays.npz']),allow_pickle=False)
assert np.array_equal(pos[-1],checkpoint['positions']);assert np.isfinite(pos).all()
savecsv('thermal_full.csv',rows)
last=[r for r in rows if 1500<r['step']<=2000]
edges=np.arange(0,L/2,.05);edges=np.append(edges,L/2);centers=(edges[:-1]+edges[1:])/2
shell=4*np.pi/3*np.diff(edges**3)
def rdf(frames):
 counts=np.zeros(len(centers))
 for p in frames:
  d=np.linalg.norm(vectors(p,cell),axis=-1);counts+=np.histogram(d[np.triu_indices(N,1)],edges)[0]*2
 return counts/(len(frames)*N*rho*shell)
late=pos[steps>1500];g=rdf(late);smooth=gaussian_filter1d(g,1)
peaks=find_peaks(smooth)[0];peak=max([i for i in peaks if 2.2<centers[i]<3.4],key=lambda i:smooth[i])
mins=find_peaks(-smooth)[0];valleys=[i for i in mins if centers[i]>centers[peak]+.35 and centers[i]<4.3]
valley=min(valleys,key=lambda i:smooth[i]);cut=float(centers[valley])
# Definition sanity checks on ideal FCC, no energy/force evaluation.
fcc=bulk('Ag','fcc',a=4.09,cubic=True).repeat((3,3,3));fq,_=q6(fcc.positions,fcc.cell.array,3.5)
assert np.allclose(fq,.5745242597,atol=1e-8);assert set(cna(fcc.positions,fcc.cell.array,3.5))=={'FCC'}
hcp=bulk('Ag','hcp',a=2.9,c=2.9*np.sqrt(8/3),orthorhombic=True).repeat((3,3,3));hq,_=q6(hcp.positions,hcp.cell.array,3.5)
assert np.allclose(hq,.4847616852,atol=1e-8);assert set(cna(hcp.positions,hcp.cell.array,3.5))=={'HCP'}
frame=[];atomrows=[]
for step,p in zip(steps,pos):
 d=np.linalg.norm(vectors(p,cell),axis=-1);np.fill_diagonal(d,np.inf);coord=np.sum(d<cut,axis=1)
 q,glob=q6(p,cell,cut);labels=cna(p,cell,cut);counts=collections.Counter(labels)
 frame.append({'step':int(step),'time_ps':step*.0005,'min_distance_A':float(d.min()),'coordination_mean':float(coord.mean()),'Q6_local_mean':float(q.mean()),'Q6_local_std':float(q.std(ddof=1)),'Q6_global':glob,'CNA_FCC':counts['FCC'],'CNA_HCP':counts['HCP'],'CNA_ICO':counts['ICOSAHEDRAL'],'CNA_OTHER':counts['OTHER_UNCLASSIFIED']})
 if step==2000:
  for i in range(N):atomrows.append({'atom_id':i+1,'Q6_local':float(q[i]),'coordination':int(coord[i]),'CNA_fixed_cutoff':labels[i],'nearest_neighbor_distance_A':float(d[i].min())})
savecsv('frame_metrics.csv',frame);savecsv('final_atom_structure_metrics.csv',atomrows)
savecsv('rdf.csv',[{'r_A':float(r),'g_r_late5frames':float(x),'smoothed_for_extrema_only':float(y),'g_r_final':float(v)} for r,x,y,v in zip(centers,g,smooth,rdf(pos[-1:]))])
# Minimum-image unwrapping of saved-frame increments (50fs spacing), not guaranteed continuous trajectories.
delta=np.diff(pos,axis=0);frac=delta@np.linalg.inv(cell);frac-=np.rint(frac);delta=frac@cell
unwrap=np.concatenate([pos[:1],pos[:1]+np.cumsum(delta,axis=0)])
disp=unwrap-unwrap[0];comdisp=disp-disp.mean(axis=1,keepdims=True)
msd=np.mean(np.sum(comdisp**2,axis=-1),axis=1);rawmsd=np.mean(np.sum(disp**2,axis=-1),axis=1)
# Late-origin MSD, averaged over all valid time origins starting at step1000.
msdrows=[]
for lag in range(1,11):
 dd=unwrap[10+lag:]-unwrap[10:-lag];dd-=dd.mean(axis=1,keepdims=True)
 msdrows.append({'lag_ps':lag*.05,'time_origins':len(dd),'MSD_multi_origin_A2':float(np.mean(np.sum(dd**2,axis=-1)))})
savecsv('msd_late_multi_origin.csv',msdrows)
savecsv('msd.csv',[{'step':int(s),'time_ps':s*.0005,'MSD_COM_removed_A2':float(m),'MSD_raw_A2':float(raw)} for s,m,raw in zip(steps,msd,rawmsd)])
savecsv('atom_migration.csv',[{'atom_id':i+1,'net_displacement_COM_removed_A':float(np.linalg.norm(comdisp[-1,i])),'saved_frame_path_length_A':float(np.linalg.norm(delta[:,i],axis=1).sum())} for i in range(N)])
metrics={'task_id':state['task_id'],'LIQUID_CANDIDATE_STATUS':'NO','thermalization_status':'THERMALIZATION_NOT_ESTABLISHED','target_temperature_K':1600,'required_final_sampling_temperature_K':1250,'duration_ps':1.0,'thermal_rows':len(rows),'trajectory_frames':len(pos),'late_window':'1500 < step <= 2000; 25 logs, 5 frames','thermal_full':stats(rows),'thermal_last500steps':stats(last),'thermal_previous500steps':stats([r for r in rows if 1000<r['step']<=1500]),'energy_definition':'sampler_energy_eV = Atoms.get_potential_energy(); kinetic_energy_eV logged for steps1020-2000; steps0-1000 derived from logged T with DOF189 and ASE kB (original log lacks KE); total is their sum. No force/energy recalculation.','RDF':{'dr_A':.05,'rmax_A':L/2,'normalization':'ordered pair counts / (frames*N*(N/V)*spherical shell volume); finite-N ideal gas baseline (N-1)/N','late_frames':5,'first_peak_A':float(centers[peak]),'first_minimum_A':cut,'late_mean_coordination':float(np.mean([f['coordination_mean'] for f in frame[-5:]]))},'final_frame':frame[-1],'CNA':{'method':'fixed-cutoff conventional CNA with MIC center-anchored local-image bond topology','cutoff_A':cut,'classifications':'FCC,HCP,ICOSAHEDRAL,OTHER_UNCLASSIFIED; BCC not assessed','PTM':'NOT_PERFORMED: ovito/PTM backend unavailable; no PTM labels invented','ideal_FCC_reference_test':'PASS Q6=0.5745242597 and 100% FCC; geometry only'},'dynamics':{'MSD_final_COM_removed_A2':float(msd[-1]),'MSD_raw_final_A2':float(rawmsd[-1]),'late_msd_slope_A2_per_ps':float(np.polyfit(steps[10:]*.0005,msd[10:],1)[0]),'max_saved_frame_displacement_A':float(np.linalg.norm(delta,axis=-1).max()),'median_net_atom_displacement_A':float(np.median(np.linalg.norm(comdisp[-1],axis=-1))),'diffusion_coefficient':None,'reason':'1ps nonstationary warmup; reconstructed MIC trajectory at 50fs intervals, insufficient diffusive regime and independent origins'},'safety':{'minimum_saved_frame_distance_A':min(f['min_distance_A'] for f in frame),'minimum_final_distance_A':frame[-1]['min_distance_A'],'runtime_gate_A':1.8,'abnormal_short_distance_in_saved_frames':any(f['min_distance_A']<1.8 for f in frame),'coverage':'21 saved frames only; not an every-integration-step distance audit'},'missing_conditions':['Not cooled/equilibrated at 1250K','Thermal and potential-energy stationarity with reliable autocorrelation not established','Insufficient long-time diffusive evidence','Only one 64-atom lineage; other five independent candidates and cross-lineage deduplication absent'],'source':{'step2000_commit':COMMIT,'archive_path':PATH,'archive_SHA256':SHA(b),'checkpoint_state_SHA256':SHA(new[prefix+'state.json']),'original_step1000_commit':OLD_COMMIT,'original_archive_SHA256':SHA(oldb),'archive_member_checks':'PASS','step1000_stitch_exact':'PASS','final_checkpoint_coordinates_exact':'PASS'},'MD_started':False,'DFT_started':False,'training_started':False,'candidate_written':False}
sensitivity=[]
for rc in [cut-.3,cut-.15,cut,cut+.15]:
 qq,gg=q6(pos[-1],cell,rc);cc=collections.Counter(cna(pos[-1],cell,rc));dd=np.linalg.norm(vectors(pos[-1],cell),axis=-1);np.fill_diagonal(dd,np.inf)
 sensitivity.append({'cutoff_A':rc,'coordination_mean':float((dd<rc).sum(axis=1).mean()),'Q6_local_mean':float(qq.mean()),'Q6_global':gg,'FCC':cc['FCC'],'HCP':cc['HCP'],'ICO':cc['ICOSAHEDRAL'],'OTHER':cc['OTHER_UNCLASSIFIED']})
metrics['CNA']['cutoff_sensitivity_final_frame']=sensitivity
metrics['CNA']['ideal_HCP_reference_test']='PASS Q6=0.4847616852 and 100% HCP; geometry only'
metrics['safety']['final_nearest_neighbor_quantiles_A']={str(q):float(np.quantile([r['nearest_neighbor_distance_A'] for r in atomrows],q)) for q in [0,.25,.5,.75,1]}
metrics['dynamics']['atoms_net_displacement_gt_1A']=int((np.linalg.norm(comdisp[-1],axis=-1)>1).sum())
metrics['dynamics']['atoms_net_displacement_gt_2A']=int((np.linalg.norm(comdisp[-1],axis=-1)>2).sum())
metrics['source']['parent_family']=source_manifest['trajectory_parent']
metrics['source']['initialization']=source_manifest['initialization']
metrics['source']['data_role']=source_manifest['fixed_identity']['row']['data_role']
metrics['source']['potential_SHA256']=source_manifest['fixed_identity']['potential']['sha256']
metrics['structure_conclusion']='Disordered short-range dense coordination observed; no predominant FCC/HCP detected by this fixed-cutoff CNA. OTHER is not a liquid label; equilibrium liquid not established.'
savecsv('cutoff_sensitivity.csv',sensitivity)
(OUT/'metrics.json').write_text(json.dumps(metrics,indent=2,sort_keys=True,allow_nan=False)+'\n')
# Figures.
t=np.array([r['step'] for r in rows]);fig,axs=plt.subplots(4,1,figsize=(9,10),sharex=True)
for ax,key,label in zip(axs,['T_K','sampler_energy_eV','kinetic_energy_eV','pressure_eV_A3'],['Temperature (K)','Potential energy (eV)','Kinetic energy (eV)','Pressure (eV/A^3)']):
 ax.plot(t,[r[key] for r in rows]);ax.set_ylabel(label);ax.axvspan(1500,2000,color='grey',alpha=.12);ax.grid(alpha=.2)
axs[0].axhline(1600,color='red',ls='--',label='Target 1600 K');axs[0].legend();axs[-1].set_xlabel('Global step (0.5 fs/step)');fig.tight_layout();fig.savefig(OUT/'thermal_full.png',dpi=160);plt.close(fig)
fig,axs=plt.subplots(2,2,figsize=(10,8));axs[0,0].plot(centers,g,label='late 5 frames');axs[0,0].plot(centers,rdf(pos[-1:]),alpha=.5,label='final frame');axs[0,0].axvline(cut,color='grey',ls='--');axs[0,0].set(xlabel='r (A)',ylabel='g(r)');axs[0,0].legend()
axs[0,1].plot(steps,[f['Q6_local_mean'] for f in frame],label='mean local Q6');axs[0,1].plot(steps,[f['Q6_global'] for f in frame],label='global Q6');axs[0,1].set(xlabel='Step',ylabel='Q6');axs[0,1].legend()
axs[1,0].plot(steps*.0005,msd,label='COM-removed, origin step0');axs[1,0].plot([r['lag_ps'] for r in msdrows],[r['MSD_multi_origin_A2'] for r in msdrows],label='late multi-origin');axs[1,0].set(xlabel='Time / lag (ps)',ylabel='MSD (A^2)');axs[1,0].legend()
axs[1,1].plot(steps,[f['min_distance_A'] for f in frame]);axs[1,1].axhline(1.8,color='red',ls='--');axs[1,1].set(xlabel='Step',ylabel='Minimum distance (A)')
fig.tight_layout();fig.savefig(OUT/'structure_dynamics.png',dpi=160);plt.close(fig)
print(json.dumps({k:metrics[k] for k in ['LIQUID_CANDIDATE_STATUS','thermal_last500steps','RDF','final_frame','dynamics','safety']},indent=2))
