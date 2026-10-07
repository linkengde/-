"""Read geometry columns only; compare proposed V17 roles without label access."""
import importlib.util,json,csv,hashlib,itertools
from pathlib import Path
import numpy as np
from ase.geometry import find_mic
from scipy.optimize import linear_sum_assignment
ROOT=Path(__file__).resolve().parents[6];OUT=Path(__file__).resolve().parent
BASE=ROOT/'research/mace-v12-transfer';WB=BASE/'coordination/reports/window-b'
FINAL=WB/'v17_split_geometry_design/final_set'
spec=importlib.util.spec_from_file_location('geometry',FINAL/'generate_split_candidates.py');g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)
sha=g.sha256
for line in (FINAL/'SHA256SUMS.txt').read_text().splitlines():
 digest,name=line.split(None,1);name=name.strip().lstrip('*')
 if name in ['generate_split_candidates.py','role_manifest.json'] or name.startswith('geometries/'):
  assert sha(FINAL/name)==digest,(name,'frozen hash mismatch')
thresholds={'cell_tolerance_A':1e-5,'exact_max_displacement_A':1e-5,'near_Ag_RMS_A':0.15,'near_framework_RMS_A':0.15,'species_minimum_floor_rule':'historical V16 same-element minimum minus0.05 A; comparative triage only'}
inventory=[];structures=[]
def load(p,pin=None):
 h=sha(p);assert pin is None or h==pin,(p,pin,h)
 inventory.append({'path':str(p.relative_to(ROOT)),'sha256':h,'expected_sha256':pin})
 return g.extxyz_geometry_frames(p),h
train=BASE/'periodic_interface_v4/mace_periodic_v16_interface_energy/data/train.extxyz'
m=json.loads((train.parent/'dataset_manifest.json').read_text());aa,h=load(train,m['output_sha256']['train']);assert len(aa)==35
for i,a in enumerate(aa): structures.append({'id':f'V16_train_{i:02}','role':'train','origin':'V16','atoms':a,'file_sha256':h,'frame':i})
for directory in ['pbe_interface_v17_parallel_targeted_acquisition','pbe_interface_v17_main_targeted_acquisition']:
 mp=BASE/'periodic_interface_v4'/directory/'input_manifest.json';inventory.append({'path':str(mp.relative_to(ROOT)),'sha256':sha(mp)})
 for r in json.loads(mp.read_text())['records']:
  atoms,h=load(mp.parent/r['input'],r['input_sha256']);assert len(atoms)==1
  structures.append({'id':r['label'],'role':'train','origin':'V17','atoms':atoms[0],'file_sha256':h,'parent':r['source_parent']})
roles=json.loads((FINAL/'role_manifest.json').read_text());inventory.append({'path':str((FINAL/'role_manifest.json').relative_to(ROOT)),'sha256':sha(FINAL/'role_manifest.json')})
for r in roles['candidates']:
 atoms,h=load(FINAL/r['file'],r['file_sha256']);assert len(atoms)==1 and g.geometry_sha(atoms[0])==r['geometry_sha256']
 structures.append({'id':r['candidate_id'],'role':r['proposed_role'],'origin':'candidate','atoms':atoms[0],'file_sha256':h,'lineage':r['lineage']})
assert len(structures)==52
# Always supplement identity mapping with chemically constrained assignment.
# Framework-median translation alignment; no rotations/cell scaling. Assignment
# starts at the identity-based common translation and iterates squared MIC costs.
def species_match(a,b):
 if g.formula_key(a)!=g.formula_key(b) or not np.array_equal(a.pbc,b.pbc) or np.max(np.abs(a.cell.array-b.cell.array))>1e-5:return None
 sy=np.array(a.get_chemical_symbols());bt=np.array(b.get_chemical_symbols());base=g.compare_geometry(a,b)
 shift=np.zeros(3)
 amap,bmap=g.atom_id_map(a),g.atom_id_map(b)
 if amap and bmap and set(amap)==set(bmap):
  ia=np.array([amap[k] for k in sorted(amap)]);ib=np.array([bmap[k] for k in sorted(amap)])
  assert np.array_equal(sy[ia],bt[ib])
  delta,_=find_mic(a.positions[ia]-b.positions[ib],a.cell,pbc=a.pbc);fw=sy[ia]!='Ag';shift=np.median(delta[fw] if fw.any() else delta,axis=0)
 for _ in range(5):
  ia=[];ib=[]
  for s in sorted(set(sy)):
   x=np.where(sy==s)[0];y=np.where(bt==s)[0]
   _,dist=find_mic((a.positions[x,None]-b.positions[None,y]-shift).reshape(-1,3),a.cell,pbc=a.pbc)
   rr,cc=linear_sum_assignment(np.asarray(dist).reshape(len(x),len(y))**2);ia.extend(x[rr]);ib.extend(y[cc])
  ia=np.array(ia);ib=np.array(ib);delta,_=find_mic(a.positions[ia]-b.positions[ib],a.cell,pbc=a.pbc);fw=sy[ia]!='Ag';new=np.median(delta[fw] if fw.any() else delta,axis=0)
  if np.linalg.norm(new-shift)<1e-10:shift=new;break
  shift=new
 delta,_=find_mic(delta-shift,a.cell,pbc=a.pbc);mag=np.linalg.norm(delta,axis=1);ag=sy[ia]=='Ag'
 ar=float(np.sqrt(np.mean(mag[ag]**2))) if ag.any() else 0.;fr=float(np.sqrt(np.mean(mag[~ag]**2))) if (~ag).any() else 0.
 return {'Ag_RMS_A':ar,'framework_RMS_A':fr,'all_RMS_A':float(np.sqrt(np.mean(mag**2))),'max_A':float(max(mag)),'exact':bool(max(mag)<=1e-5),'near':bool(ar<=.15 and fr<=.15),'algorithm':'species Hungarian squared MIC; framework median translation; 5 iterations, no rotation; local assignment audit, not exhaustive symmetry optimization'}
pairs=[]
for x,y in itertools.combinations(structures,2):
 if x['role']==y['role']:continue
 a,b=x['atoms'],y['atoms'];identity=g.compare_geometry(a,b);species=species_match(a,b)
 rec={'left':x['id'],'right':y['id'],'left_role':x['role'],'right_role':y['role'],'left_origin':x['origin'],'exact_file_hash':x['file_sha256']==y['file_sha256'],'identity_comparison':identity,'species_comparison':species,'exact_geometry':bool(identity and identity['exact_duplicate']) or bool(species and species['exact']),'near_geometry':bool(identity and identity['near_duplicate']) or bool(species and species['near'])}
 amap,bmap=g.atom_id_map(a),g.atom_id_map(b)
 if amap and bmap and set(amap)==set(bmap):
  ids=sorted(amap);delta,_=find_mic(b.positions[[bmap[k] for k in ids]]-a.positions[[amap[k] for k in ids]],a.cell,pbc=a.pbc)
  rec['same_ID_displacements_A']={str(k):d.tolist() for k,d in zip(ids,delta) if np.linalg.norm(d)>1e-5};rec['same_ID_parent_family']=True
 else:rec['same_ID_parent_family']=False
 pairs.append(rec)
assert len(pairs)==405
floors={}
for a in aa:
 for pair,r in g.pair_minima(a).items():floors[pair]=min(floors.get(pair,float('inf')),r['distance_A'])
sanity=[]
for s in structures:
 minima=g.pair_minima(s['atoms']);flags=[p for p,r in minima.items() if r['distance_A']<floors.get(p,0)-.05]
 sanity.append({'id':s['id'],'role':s['role'],'species_minima':minima,'below_historical_floor_minus_0p05':flags})
counts={'structures':52,'cross_role_pairs':len(pairs),'exact_file_hits':sum(r['exact_file_hash'] for r in pairs),'exact_geometry_hits':sum(r['exact_geometry'] for r in pairs),'near_geometry_hits':sum(r['near_geometry'] for r in pairs),'new_V17_training_near_hits':sum(r['near_geometry'] and 'V17' in [r['left_origin'],r['right_origin']] for r in pairs)}
(OUT/'audit.json').write_text(json.dumps({'thresholds':thresholds,'counts':counts,'source_inventory':inventory,'pairs':pairs,'species_sanity':sanity,'historical_species_minima_A':floors},indent=2)+'\n')
with (OUT/'pairwise.csv').open('w') as f:
 fields=['left','right','left_role','right_role','left_origin','right_origin','exact_file_hash','exact_geometry','near_geometry','same_ID_parent_family'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows({k:r.get(k) for k in fields} for r in pairs)
hits=[r for r in pairs if r['exact_geometry'] or r['near_geometry']]
(OUT/'report.md').write_text(f'''# V17 geometry-only split leakage audit

## Scope and result

35 V16 training geometries, eight V17 diagnostic inputs, three development inputs and six withheld inputs: 52 structures, 405 cross-role pairs. Exact file matches: {counts['exact_file_hits']}; normalized exact geometry matches: {counts['exact_geometry_hits']}; near matches: {counts['near_geometry_hits']}; near matches involving new V17 training: {counts['new_V17_training_near_hits']}.

## Algorithms fixed before interpretation

Cell tolerance 1e-5 A, exact maximum displacement 1e-5 A, near Ag and framework RMS each <=0.15 A. Reused split-design identity mapping is supplemented by species-preserving Hungarian assignment on squared periodic MIC distances, framework-median translation alignment and up to five iterations. No rotation or cell scaling. The assignment is a local geometric screen, not proof against every global symmetry. ID-matched displacement vectors are reported independently. Shared IDs show common lineage, not automatically leakage. File hash equality across multi-frame training files and single-frame candidates is a weak check; normalized coordinate checks are the substantive comparison.

## Interpretation and recommendation

{json.dumps(counts)}

{'Cross-role exact/near matches require A review before role freezing; do not silently accept role independence.' if hits else 'No exact/near cross-role hit was found under both stated algorithms. Newly added paired-training inputs do not introduce a hit under these thresholds. A may retain proposed roles subject to provenance/label review.'}

All candidates inherit existing small-cluster parent motifs. Absence of duplicate geometries establishes only local separation, not new morphology or thermal coverage. Same-ID maps in audit.json quantify correlated families separately from leakage. Species minima are compared against the lowest observed historical V16 distance for that element pair minus0.05 A; these flags are comparative triage, not chemical validity or a universal 1.75 A cutoff. Full species minima and flags are provided for every geometry.

Only geometry columns were parsed. Candidate hashes match the frozen role manifest. No development/withheld calculation output, label, summary, log or model score was opened; no calculation/inference/training was started. No datasets or role assignments were changed.

## Reproduce

Run audit.py with the existing GPAW Python environment; it imports the committed geometry parser without running its generator. SHA256SUMS.txt covers the report artifacts.
''')
(OUT/'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='SHA256SUMS.txt'))
print(counts)
