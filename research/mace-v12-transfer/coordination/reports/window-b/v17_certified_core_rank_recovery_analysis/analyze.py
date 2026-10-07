import json,csv,itertools,importlib.util
from pathlib import Path
from fractions import Fraction
import numpy as np
ROOT=Path(__file__).resolve().parents[6];OUT=Path(__file__).resolve().parent;WB=ROOT/'research/mace-v12-transfer/coordination/reports/window-b'
spec=importlib.util.spec_from_file_location('geom',WB/'v17_split_geometry_design/final_set/generate_split_candidates.py');g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)
def rank(rows):
 a=[[Fraction(int(x)) for x in r] for r in rows];i=0
 for j in range(4):
  pivot=next((k for k in range(i,len(a)) if a[k][j]),None)
  if pivot is None:continue
  a[i],a[pivot]=a[pivot],a[i];p=a[i][j];a[i]=[x/p for x in a[i]]
  for k in range(i+1,len(a)):
   p=a[k][j];a[k]=[x-p*y for x,y in zip(a[k],a[i])]
  i+=1
 return i
species=['Ag','C','Si','Ti']
def counts(a):return [a.get_chemical_symbols().count(s) for s in species]
poolpath=WB/'v17_training_pool_provenance_audit/ledger.json';pool=json.loads(poolpath.read_text())
train=ROOT/'research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v16_interface_energy/data/train.extxyz';atoms=g.extxyz_geometry_frames(train);assert len(atoms)==35
matrix=[counts(a) for a in atoms];direct=[r['frame'] for r in pool['rows'] if r['evidence_status']=='directly_archive_verified'];partial=[r['frame'] for r in pool['rows'] if r['evidence_status']=='partially_supported'];core=[matrix[i] for i in direct];sources={str(train.relative_to(ROOT)):g.sha256(train),str(poolpath.relative_to(ROOT)):g.sha256(poolpath)}
for directory in ['pbe_interface_v17_main_targeted_acquisition','pbe_interface_v17_parallel_targeted_acquisition']:
 mp=ROOT/'research/mace-v12-transfer/periodic_interface_v4'/directory/'input_manifest.json';sources[str(mp.relative_to(ROOT))]=g.sha256(mp)
 for r in json.loads(mp.read_text())['records']:
  p=mp.parent/r['input'];assert g.sha256(p)==r['input_sha256'];core.append(counts(g.extxyz_geometry_frames(p)[0]));sources[str(p.relative_to(ROOT))]=g.sha256(p)
assert len(core)==29 and len(partial)==14
solutions=[];tested=0
for k in range(1,15):
 for subset in itertools.combinations(partial,k):
  tested+=1;rows=core+[matrix[i] for i in subset]
  if rank(rows)==4:
   solutions.append({'frames':list(subset),'rank':4,'singular_values':np.linalg.svd(np.array(rows,dtype=float),compute_uv=False).tolist(),'rows':[{'frame':i,'config_type':pool['rows'][i]['config_type'],'formula':atoms[i].get_chemical_formula(),'provenance':'partially_supported','method':pool['rows'][i].get('method_evidence'),'unresolved':pool['rows'][i].get('original_run_evidence')} for i in subset]})
 if solutions:break
result={'species_order':species,'algorithm':'exact Fraction Gaussian elimination; singular values descriptive only','core_count':29,'core_rank':rank(core),'core_singular_values':np.linalg.svd(np.array(core,float),compute_uv=False).tolist(),'core_integer_matrix':core,'partial_integer_matrix':{str(i):matrix[i] for i in partial},'minimum_additional_rows':k,'all_minimum_solutions':solutions,'subsets_tested_until_minimum_found':tested,'source_hashes':sources}
(OUT/'rank_analysis.json').write_text(json.dumps(result,indent=2)+'\n')
with (OUT/'rank_sensitivity.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['frames','configs','formulas','count','exact_rank']);w.writerow(['core','','',29,rank(core)])
 for r in solutions:w.writerow([str(r['frames']),';'.join(x['config_type'] for x in r['rows']),';'.join(x['formula'] for x in r['rows']),29+len(r['frames']),4])
text=f"Core29 exact rank={rank(core)}; minimum additional rows={k}; all minimum solutions="+str([r['frames'] for r in solutions])
(OUT/'interim_note.md').write_text(text+'\n')
(OUT/'report.md').write_text('# Certified-core rank recovery\n\n'+text+'\n\nExact integer compositions were rederived with a geometry-only parser. All subsets below the first successful cardinality and every subset at that cardinality were enumerated. Rational Gaussian elimination avoids floating tolerance choices. Source hashes, matrices, singular values and all candidate frame identities are recorded in JSON. No energy/force values from partial frames were used.\n\nEach proposed addition remains partially supported until original-run records are recovered or A authorizes an exact fixed-geometry relabel. Keeping it provisionally restores algebraic rank but does not certify historical method/convention. A-approved PW-PBE500/Gamma/Fermi0.1 relabels must retain atom order/cell/PBC and document native/free energies separately; Gamma adequacy for bulk unary cells requires review. These are recommendations, not launched jobs. Full rank does not prove force quality, chemical coverage or transferability. Removal of the sole added row returns this core to rank3. No production integration before A review.\n')
(OUT/'SHA256SUMS.txt').write_text(''.join(f'{g.sha256(p)}  {p.name}\n' for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='SHA256SUMS.txt'))
print(text)
