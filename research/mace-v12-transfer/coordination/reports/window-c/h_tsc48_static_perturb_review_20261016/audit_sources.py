from pathlib import Path
import subprocess,hashlib,json,re,shlex,itertools,csv
import numpy as np
COMMIT='e923008453d71077421b651cf19c6705420d2356'
ROOT='research/mace-v12-transfer/coordination/reports/'
BASE='research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_pure_phase_controls/'
CAND=ROOT+'chatgpt-nep-ace-dft-repair-20261016/NEW_DFT_GEOMETRY_CANDIDATES/structures/TSC48_static_perturb_std0p025_seed2026101603_UNLABELED.extxyz'
OUT=Path(ROOT+'window-c/h_tsc48_static_perturb_review_20261016'); (OUT/'evidence').mkdir(parents=True,exist_ok=True)
records=[]
def get(path,save=None,commit=COMMIT):
 b=subprocess.check_output(['git','show',f'{commit}:{path}']); records.append({'commit':commit,'path':path,'sha256':hashlib.sha256(b).hexdigest(),'bytes':len(b)})
 if save:(OUT/'evidence'/save).write_bytes(b)
 return b
def parse(b):
 lines=b.decode().splitlines(); n=int(lines[0]); hdr=dict(x.split('=',1) for x in shlex.split(lines[1]) if '=' in x)
 props=hdr['Properties'].split(':'); offsets={}; offset=0
 for i in range(0,len(props),3): offsets[props[i]]=(offset,int(props[i+2]));offset+=int(props[i+2])
 rows=[l.split() for l in lines[2:2+n]]; assert len(rows)==n and all(len(r)==offset for r in rows)
 cell=np.array(list(map(float,hdr['Lattice'].split()))).reshape(3,3)
 pos=np.array([list(map(float,r[offsets['pos'][0]:offsets['pos'][0]+3])) for r in rows])
 return {'n':n,'header':hdr,'rows':rows,'offsets':offsets,'cell':cell,'positions':pos,'symbols':[r[0] for r in rows],'pbc':hdr['pbc'].split()}
def mic(d,cell):
 frac=d@np.linalg.inv(cell); central=np.rint(frac); shifts=np.array(list(itertools.product([-1,0,1],repeat=3)))
 vectors=d[...,None,:]-(central[...,None,:]+shifts)@cell
 norms=np.sum(vectors*vectors,axis=-1); ix=np.argmin(norms,axis=-1)
 return np.take_along_axis(vectors,ix[...,None,None],axis=-2)[...,0,:]
raw=get(CAND,'candidate_UNCHANGED.extxyz'); assert hashlib.sha256(raw).hexdigest()=='da3820c477023e0bb423f87550ed82ae0c2a571f53d294c055f78d7d052d5e8b'
a=parse(raw); p=parse(get(BASE+'inputs/Ti3SiC2_baseline_PROPOSAL.extxyz','parent_UNCHANGED.extxyz'))
assert a['n']==48 and a['symbols']==p['symbols'] and a['pbc']==['T','T','T'] and np.max(abs(a['cell']-p['cell']))<1e-9
ids=[int(r[a['offsets']['parent_row_id'][0]]) for r in a['rows']];assert ids==list(range(1,49))
d=mic(a['positions']-p['positions'],a['cell']); dn=np.linalg.norm(d,axis=1)
i,j=np.triu_indices(48,1); distances=np.linalg.norm(mic(a['positions'][j]-a['positions'][i],a['cell']),axis=1); k=int(np.argmin(distances))
labels=['Ti3SiC2_baseline_k10x10x4_v19']+[f'Ti3SiC2_{element}4f_id{atom}_z_{sgn}020A_k10x10x4_v19' for element,atom in [('Ti','0003'),('C','0009')] for sgn in ['m','p']]
comparisons=[]
for label in labels:
 folder=BASE+'calculations/'+label+'/'
 lb=get(folder+label+'_PW_PBE.extxyz'); la=parse(lb)
 summ=json.loads(get(folder+'summary.json'));ver=json.loads(get(folder+'verification.json'))
 assert summ['scf_converged'] is True and ver['status']=='PASS'
 assert hashlib.sha256(lb).hexdigest()==ver['sha256'][label+'_PW_PBE.extxyz']
 assert la['symbols']==a['symbols'] and np.max(abs(la['cell']-a['cell']))<1e-9
 v=mic(a['positions']-la['positions'],a['cell']); vv=np.linalg.norm(v,axis=1)
 comparisons.append({'label':label,'source_sha256':hashlib.sha256(lb).hexdigest(),'scf_converged':True,'verification':'PASS','rms_displacement_A':float(np.sqrt(np.mean(vv**2))),'max_displacement_A':float(vv.max()),'exact_duplicate_1e_6A':bool(vv.max()<1e-6)})
assert not any(c['exact_duplicate_1e_6A'] for c in comparisons)
mesh=json.loads(get(ROOT+'window-a/v19_Ti3SiC2_bulk_mesh_comparison/k8x8x4_k10x10x4_comparison.json','parent_mesh_comparison.json'))
for side,tag in [('left','k8x8x4'),('right','k10x10x4')]:
 folder=BASE+f'calculations/Ti3SiC2_baseline_{tag}_v19/'
 for name,expected in mesh['archive_verification'][side+'_files'].items():
  b=get(folder+name, 'parent_k10_gpaw.log' if tag=='k10x10x4' and name=='gpaw.log' else None)
  assert hashlib.sha256(b).hexdigest()==expected
get(ROOT+'chatgpt-nep-ace-dft-repair-20261016/TSC48_PREFLIGHT_20261016/H_TSC48_STATIC_PERTURB_SINGLEPOINT_PREFLIGHT.md','H_ORIGINAL_PREFLIGHT.md')
get(ROOT+'chatgpt-nep-ace-dft-repair-20261016/NEW_DFT_GEOMETRY_CANDIDATES/README_聊天生成与计算状态.md','GENERATOR_SOURCE_README.md')
assert get(CAND,commit='cc65a3e6b75be28c06fda8608b4064662d714117')==raw
assert get(BASE+'inputs/Ti3SiC2_baseline_PROPOSAL.extxyz',commit='4ee3daab860a846fec4ef0c0ca4b031a36cfc87e')==get(BASE+'inputs/Ti3SiC2_baseline_PROPOSAL.extxyz')
for name in ['VERSIONS_AND_RUNTIME.txt','MPI_4RANK_TEST.txt','GPAW_MPI_LINKAGE.txt']:
 get(ROOT+'chatgpt-nep-ace-dft-repair-20261016/CLOUD_DFT_RESULTS_20261016/H_Ti_gap_32_k1x10x8/preflight/'+name,name)
role=get(ROOT+'window-e/ti3sic2_backbone_e2_parameter_audit_20261014/role_freeze.json','E2_ROLE_FREEZE.json',commit='f3e09b7a3b595d21c1fd02cd72a6dda9cb74b8c1');assert hashlib.sha256(role).hexdigest()=='709d84b6f7205d3d437bae92dc1aaecd81278311f807b805e06c1b7d22acc27f'
metrics={'candidate_sha256':hashlib.sha256(raw).hexdigest(),'atom_count':48,'composition':{s:a['symbols'].count(s) for s in sorted(set(a['symbols']))},'cell_A':a['cell'].tolist(),'pbc':a['pbc'],'ids':'parent_row_id 1..48','parent_rms_displacement_A':float(np.sqrt(np.mean(dn**2))),'max_displacement_A':float(dn.max()),'minimum_periodic_distance_A':float(distances[k]),'minimum_pair_row_ids':[int(i[k]+1),int(j[k]+1)],'displacement_mean_xyz_A':d.mean(axis=0).tolist(),'species_displacement_rms_A':{s:float(np.sqrt(np.mean(dn[np.array(a['symbols'])==s]**2))) for s in sorted(set(a['symbols']))},'existing_five_comparison':comparisons,'parent_mesh_comparison':mesh,'candidate_kpoint_convergence':'UNVERIFIED','scientific_status':'NEEDS_REVISION','data_role':'HOLD','kgrid_candidate':'REVISE','resource_ready':False,'dft_execution_ready':False,'execution_authorized':False}
(OUT/'metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
with (OUT/'SOURCE_SHA256.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['commit','path','sha256','bytes'],lineterminator='\n');w.writeheader();w.writerows(records)
print(json.dumps({k:v for k,v in metrics.items() if k not in ['parent_mesh_comparison','existing_five_comparison']},indent=2))
print('SOURCE_RECORDS_VERIFIED',len(records))
