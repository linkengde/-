"""Offline analysis of archived Ag64 step 3000; reads coordinates and logs only."""
import csv, json, hashlib, tarfile, pathlib, io, subprocess, collections
import numpy as np
from scipy.special import sph_harm_y
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=pathlib.Path(__file__).resolve().parent
REPO=pathlib.Path('/tmp/cloud-c-v19-design')
D_COMMIT='2ad687d7f74e035a9d2274e7d1e9a70577eaa76f'
P2000='codex_segment_archives/D/C-D-AGL-1250K-rho097-seed01/step-000002000-c3a88a759e0bd258e74632aba0a68e4f1a765fcbe8324034d5bf127b5aa69965.tar.gz'
P3000='codex_segment_archives/D/C-D-AGL-1250K-rho097-seed01/step-000003000-segment-000002000-000003000-cd00e213b89242fa9501aacbadd5cec8c23e938a61aeebbb2eb138997deea7d2.tar.gz'
H2000='c3a88a759e0bd258e74632aba0a68e4f1a765fcbe8324034d5bf127b5aa69965'
H3000='cd00e213b89242fa9501aacbadd5cec8c23e938a61aeebbb2eb138997deea7d2'
TASK='C-D-AGL-1250K-rho097-seed01'

def sha(b):return hashlib.sha256(b).hexdigest()
def blob(path):return subprocess.check_output(['git','show',D_COMMIT+':'+path],cwd=REPO)
def unpack(b):
  with tarfile.open(fileobj=io.BytesIO(b),mode='r:gz') as tf:
    return {m.name:mfile.read() for m in tf.getmembers() if m.isfile() for mfile in [tf.extractfile(m)]}
def nested(m):
  name=next(x for x in m if x.endswith('/MANIFEST.json'))
  pref=name.rsplit('/',1)[0]+'/'
  return pref,json.loads(m[name])
def rows_csv(b):
  return [{k:(v if k=='phase' else float(v)) for k,v in r.items()} for r in csv.DictReader(io.StringIO(b.decode()))]
def writecsv(name,rr):
  with (ROOT/name).open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rr for k in r)),lineterminator='\n');w.writeheader();w.writerows(rr)
def vectors(pos,cell):
  d=pos[None,:,:]-pos[:,None,:];s=d@np.linalg.inv(cell);s-=np.rint(s);return s@cell
def order(pos,cell,cut):
  v=vectors(pos,cell);dist=np.linalg.norm(v,axis=-1);qs=[];bond=[]
  for i in range(len(pos)):
    vv=v[i,(dist[i]>0)&(dist[i]<cut)];rr=np.linalg.norm(vv,axis=1)
    theta=np.arccos(np.clip(vv[:,2]/rr,-1,1));phi=np.arctan2(vv[:,1],vv[:,0])
    y=np.array([sph_harm_y(6,m,theta,phi) for m in range(-6,7)])
    q=y.mean(axis=1);qs.append(np.sqrt(4*np.pi/13*np.sum(np.abs(q)**2)));bond.append(y)
  gb=np.concatenate(bond,axis=1).mean(axis=1)
  return np.asarray(qs),float(np.sqrt(4*np.pi/13*np.sum(np.abs(gb)**2)))
def cna(pos,cell,cut):
  v=vectors(pos,cell);d=np.linalg.norm(v,axis=-1);labels=[]
  def longest(adj):
    best=0
    def walk(start,u,seen,n):
      nonlocal best
      best=max(best,n)
      for z in np.flatnonzero(adj[u]):
        if z==start and n>=2:best=max(best,n+1)
        elif z not in seen:walk(start,z,seen|{z},n+1)
    for z in range(len(adj)):walk(z,z,{z},0)
    return best
  for i in range(len(pos)):
    neigh=np.flatnonzero((d[i]>0)&(d[i]<cut));sig=[]
    for j in neigh:
      common=[k for k in neigh if k!=j and np.linalg.norm(v[i,k]-v[i,j])<cut]
      if not common:sig.append((0,0,0));continue
      x=v[i,common];dd=np.linalg.norm(x[:,None]-x[None,:],axis=-1);adj=(dd>0)&(dd<cut)
      sig.append((len(common),int(adj.sum()//2),longest(adj)))
    ct=collections.Counter(sig)
    if len(neigh)==12 and ct[(4,2,1)]==12:lab='FCC'
    elif len(neigh)==12 and ct[(4,2,1)]==6 and ct[(4,2,2)]==6:lab='HCP'
    elif len(neigh)==12 and ct[(5,5,5)]==12:lab='ICO'
    else:lab='OTHER'
    labels.append(lab)
  return labels
def rdf(frames,cell,dr=.05):
  L=np.linalg.norm(cell[0]);edges=np.arange(0,L/2+dr/2,dr);cent=(edges[:-1]+edges[1:])/2
  shell=4*np.pi/3*np.diff(edges**3);rho=len(frames[0])/abs(np.linalg.det(cell));counts=np.zeros(len(cent))
  for p in frames:
    d=np.linalg.norm(vectors(p,cell),axis=-1);counts+=np.histogram(d[np.triu_indices(len(p),1)],edges)[0]*2
  g=counts/(len(frames)*len(frames[0])*rho*shell)
  smooth=gaussian_filter1d(g,1);peaks=find_peaks(smooth)[0]
  peak=max((i for i in peaks if 2.2<cent[i]<3.4),key=lambda i:smooth[i])
  mins=find_peaks(-smooth)[0];val=[i for i in mins if cent[i]>cent[peak]+.35 and cent[i]<4.3]
  valley=min(val,key=lambda i:smooth[i]) if val else None
  return cent,g,smooth,int(peak),None if valley is None else int(valley)
def thermal_stat(rr):
  x=np.array([r['T_K'] for r in rr]);t=np.array([r['time_fs'] for r in rr]);slope=np.polyfit(t-t.mean(),x,1)[0]*100
  return {'n':len(x),'mean_K':float(x.mean()),'std_sample_K':float(x.std(ddof=1)),'min_K':float(x.min()),'max_K':float(x.max()),'ols_slope_K_per_100fs':float(slope),'target_K':float(rr[-1]['target_temperature_K'])}
def structure(frames,steps,cell,cut):
  out=[]
  for p,s in zip(frames,steps):
    d=np.linalg.norm(vectors(p,cell),axis=-1);np.fill_diagonal(d,np.inf)
    q,g=order(p,cell,cut);cc=collections.Counter(cna(p,cell,cut))
    out.append({'step':int(s),'time_ps':float(s*.0005),'minimum_distance_A':float(d.min()),'coordination_mean_at_fixed_2000_cutoff':float((d<cut).sum(axis=1).mean()),'coordination_std':float((d<cut).sum(axis=1).std(ddof=1)),'Q6_local_mean_at_fixed_2000_cutoff':float(q.mean()),'Q6_local_std':float(q.std(ddof=1)),'Q6_global_at_fixed_2000_cutoff':g,'CNA_FCC':cc['FCC'],'CNA_HCP':cc['HCP'],'CNA_ICO':cc['ICO'],'CNA_OTHER':cc['OTHER']})
  return out

b2=blob(P2000);b3=blob(P3000);assert sha(b2)==H2000 and sha(b3)==H3000
m2=unpack(b2);m3=unpack(b3);pr2,man2=nested(m2);pr3,man3=nested(m3)
for p,m in [(pr2,man2),(pr3,man3)]:
  for n,h in m['members_sha256'].items():assert sha((m2 if p==pr2 else m3)[p+n])==h
  checks=(m2 if p==pr2 else m3)[p+'SHA256SUMS.txt'].decode().splitlines()
  for line in checks:
    h,n=line.split('  ',1);assert sha((m2 if p==pr2 else m3)[p+n])==h
s2=json.loads(m2[pr2+'state.json']);s3=json.loads(m3[pr3+'state.json']);assert s2['step']==2000 and s3['step']==3000
assert s2['task_id']==s3['task_id']==TASK
assert sha(m2[pr2+'state.json'])=='5d03864f07a98f177c3c0cbdcb850158ed89a64fa2a6ab7f879e8c7ab003be2a'
assert sha(m3[pr3+'state.json'])=='a3c6db34bb5828d46288c61a6065cbcf4da3bd23e99c43f97ae943d4f9d5111a'
a2=np.load(io.BytesIO(m2[pr2+'trajectory.npz']),allow_pickle=False);a3=np.load(io.BytesIO(m3[pr3+'trajectory.npz']),allow_pickle=False)
assert np.array_equal(a2['steps'],a3['steps'][:len(a2['steps'])]) and np.array_equal(a2['positions'],a3['positions'][:len(a2['steps'])])
assert np.array_equal(a2['atom_ids'],a3['atom_ids']) and np.array_equal(a2['cell'],a3['cell']) and np.array_equal(a2['pbc'],a3['pbc'])
assert np.array_equal(a3['steps'],np.arange(1000,3001,100))
pos=a3['positions'];steps=a3['steps'];cell=a3['cell'];assert np.isfinite(pos).all()
t2=rows_csv(m2[pr2+'thermal.csv']);t3=rows_csv(m3[pr3+'thermal.csv']);assert all(r['step'] in {x['step'] for x in t3} for r in t2)
map3={int(r['step']):r for r in t3};assert all(r==map3[int(r['step'])] for r in t2)
# Two comparable 500-step windows, left endpoint excluded; logs are every 20 global steps.
w2=[r for r in t2 if 1500<r['step']<=2000];w3=[r for r in t3 if 2500<r['step']<=3000]
assert len(w2)==len(w3)==25
st2=thermal_stat(w2);st3=thermal_stat(w3)
# per-window RDF, peak/minimum, mean coordination at each own first minimum
windows={'step2000':(pos[(steps>1500)&(steps<=2000)],steps[(steps>1500)&(steps<=2000)]),'step3000':(pos[(steps>2500)&(steps<=3000)],steps[(steps>2500)&(steps<=3000)])}
rdfs={};win_struct={}
for label,(fr,ss) in windows.items():
  centers,g,sm,ip,iv=rdf(fr,cell)
  cut=float(centers[iv]) if iv is not None else None
  metr=structure(fr,ss,cell,cut) if cut else []
  rdfs[label]={'frames':len(fr),'first_peak_A':float(centers[ip]),'first_valley_A':cut,'peak_g':float(sm[ip]),'valley_g':None if iv is None else float(sm[iv]),'coordination_mean_window_cutoff':float(np.mean([x['coordination_mean_at_fixed_2000_cutoff'] for x in metr])) if metr else None,'Q6_local_mean_window_cutoff':float(np.mean([x['Q6_local_mean_at_fixed_2000_cutoff'] for x in metr])) if metr else None,'Q6_local_std_window_cutoff':float(np.mean([x['Q6_local_std'] for x in metr])) if metr else None,'Q6_global_last_frame_window_cutoff':metr[-1]['Q6_global_at_fixed_2000_cutoff'] if metr else None,'CNA_last_frame':{k:metr[-1][k] for k in ['CNA_FCC','CNA_HCP','CNA_ICO','CNA_OTHER']} if metr else None}
  win_struct[label]=metr
# apples-to-apples structure order and coordination: common cutoff from 2000 RDF valley
common=rdfs['step2000']['first_valley_A'];comps={}
for label,(fr,ss) in windows.items():
  comps[label]=structure(fr,ss,cell,common)
# MSD COM-removed based on MIC-unwrapped saved-frame increments. Compare same start at step1000.
delta=np.diff(pos,axis=0);frac=delta@np.linalg.inv(cell);frac-=np.rint(frac);delta=frac@cell
unwrapped=np.concatenate([pos[:1],pos[:1]+np.cumsum(delta,axis=0)])
disp=unwrapped-unwrapped[0];com=disp-disp.mean(axis=1,keepdims=True);msd=np.mean(np.sum(com**2,axis=2),axis=1)
idx={int(s):i for i,s in enumerate(steps)}
# also segment-increment MSD from step2000, subtract position at 2000 and initial COM
ix2=idx[2000];seg=com[ix2:]-com[ix2] # subtracts vector origin; CM already removed
seg_msd=np.mean(np.sum(seg**2,axis=2),axis=1)
msd_summary={'reference_step':1000,'step2000_MSD_A2':float(msd[idx[2000]]),'step3000_MSD_A2':float(msd[idx[3000]]),'change_2000_to_3000_A2':float(msd[idx[3000]]-msd[idx[2000]]),'segment2000_to3000_MSD_A2':float(seg_msd[-1]),'segment_slope_A2_per_ps':float(np.polyfit((steps[ix2:]-2000)*.0005,seg_msd,1)[0]),'maximum_saved_frame_single_atom_step_A':float(np.linalg.norm(delta,axis=2).max()),'diffusion_coefficient_reported':None,'limitation':'50fs saved-frame sampling; finite short trajectory; no claim of diffusive regime'}
writecsv('thermal_windows.csv',[{'comparison':'1500<step<=2000','metric':json.dumps(st2,sort_keys=True)},{'comparison':'2500<step<=3000','metric':json.dumps(st3,sort_keys=True)}])
writecsv('rdf_comparison.csv',[{'window':k,'r_A':float(r),'g_r':float(v)} for k,(fr,ss) in windows.items() for r,v in zip(rdfs[k].get('_r',[]) or rdf(fr,cell)[0],rdf(fr,cell)[1])])
writecsv('structure_frames.csv',[{'window':k,**row} for k,rr in comps.items() for row in rr])
writecsv('msd_curve.csv',[{'step':int(s),'time_ps':float(s*.0005),'MSD_from_step1000_COM_removed_A2':float(x),'MSD_from_step2000_COM_removed_A2':float(y) if s>=2000 else ''} for s,x,y in zip(steps,msd,np.r_[np.full(ix2,np.nan),seg_msd])])
# thermal graph and structural panels
fig,ax=plt.subplots(figsize=(9,4));
for rr,lab in [(w2,'1500<step<=2000'),(w3,'2500<step<=3000')]:ax.plot([r['step'] for r in rr],[r['T_K'] for r in rr],'.-',label=lab)
ax.axhline(1600,color='k',ls='--',label='Target 1600K');ax.set(xlabel='Global step',ylabel='Temperature (K)');ax.grid(alpha=.2);ax.legend();fig.tight_layout();fig.savefig(ROOT/'temperature_comparison.png',dpi=160);plt.close(fig)
fig,axs=plt.subplots(2,2,figsize=(9,7))
for k,(fr,ss) in windows.items():
 c,g,*_=rdf(fr,cell);axs[0,0].plot(c,g,label=k)
axs[0,0].set(xlabel='r (Å)',ylabel='g(r)');axs[0,0].legend()
for k,rr in comps.items():axs[0,1].plot([x['step'] for x in rr],[x['Q6_global_at_fixed_2000_cutoff'] for x in rr],'.-',label=k)
axs[0,1].set(xlabel='step',ylabel='global Q6 (common cutoff)');axs[0,1].legend()
axs[1,0].plot(steps,msd,label='from step1000');axs[1,0].plot(steps[ix2:],seg_msd,label='from step2000');axs[1,0].set(xlabel='global step',ylabel='MSD (Å²)');axs[1,0].legend()
for k,rr in comps.items():axs[1,1].plot([x['step'] for x in rr],[x['coordination_mean_at_fixed_2000_cutoff'] for x in rr],'.-',label=k)
axs[1,1].set(xlabel='step',ylabel=f'coordination, r<{common:.2f} Å');axs[1,1].legend()
fig.tight_layout();fig.savefig(ROOT/'structure_comparison.png',dpi=160);plt.close(fig)
result={'task_id':TASK,'source_D_commit':D_COMMIT,'step2000_archive_sha256':H2000,'step3000_archive_sha256':H3000,'step3000_checkpoint_state_sha256':sha(m3[pr3+'state.json']),'step3000_checkpoint_arrays_sha256':sha(m3[pr3+'state_arrays.npz']),'archive_and_member_SHA':'PASS','task_identity_and_step_bounds':'PASS','trajectory_overlap_1000_to2000_exact':'PASS','step3000_cell_pbc_atom_ids_match':'PASS','sampling_temperature_target_K':1600,'liquid_sampling_target_K':1250,'temperature_comparison':{'1500<step<=2000':st2,'2500<step<=3000':st3,'mean_change_K':st3['mean_K']-st2['mean_K'],'slope_change_K_per100fs':st3['ols_slope_K_per_100fs']-st2['ols_slope_K_per_100fs']},'RDF_and_order_comparison':rdfs,'common_first_valley_cutoff_A_from_step2000':common,'common_cutoff_Q6_coordination_frame_data':comps,'MSD':msd_summary,'CNA':'Fixed-cutoff conventional CNA; categories reported per final frame. PTM unavailable in this archived workflow and not run. OTHER does not by itself mean liquid.','LIQUID_CANDIDATE_STATUS':'PRELIQUID','classification_reason':'At 1600K hot-stage, late-window structure appears disordered/mobile if metrics support it, but target 1250K thermalization and final-state liquid criteria are not met. PRELIQUID denotes provisional hot-stage liquid-like state only; not frozen candidate or DFT data.','DFT_or_MD_started_by_C':False}
(ROOT/'metrics.json').write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
print(json.dumps({'temperature_comparison':result['temperature_comparison'],'RDF_and_order_comparison':result['RDF_and_order_comparison'],'MSD':msd_summary,'status':result['LIQUID_CANDIDATE_STATUS']},indent=2))
