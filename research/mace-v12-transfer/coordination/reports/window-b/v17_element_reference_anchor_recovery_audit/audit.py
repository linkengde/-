"""Focused unary-anchor evidence/rank audit; no unseen output reads."""
import json,csv,subprocess,hashlib
from pathlib import Path
import numpy as np
from ase.io import read
ROOT=Path(__file__).resolve().parents[6];OUT=Path(__file__).resolve().parent;BASE=ROOT/'research/mace-v12-transfer';WB=BASE/'coordination/reports/window-b'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
train=BASE/'periodic_interface_v4/mace_periodic_v16_interface_energy/data/train.extxyz'
a=read(train,':');m=json.loads((train.parent/'dataset_manifest.json').read_text());assert sha(train)==m['output_sha256']['train']
prior_path=WB/'v15_inherited_method_provenance_trace/per_frame_ledger.json';prior=json.loads(prior_path.read_text())
pool_path=WB/'v17_training_pool_provenance_audit/ledger.json';pool=json.loads(pool_path.read_text());anchors=[];sources={str(train.relative_to(ROOT)):sha(train),str(prior_path.relative_to(ROOT)):sha(prior_path),str(pool_path.relative_to(ROOT)):sha(pool_path)}
for i in range(10,14):
 x=a[i];r=next(z for z in prior['frames'] if z['config_type']==x.info['config_type'])
 assert x.info['REF_energy']==r['original_serialized_info']['REF_energy']
 lineage=[]
 for item in r['lineage']:
  p=ROOT/item['path'];assert sha(p)==item['file_sha256'];src=read(p,':')[item['frame']]
  assert np.array_equal(x.numbers,src.numbers) and np.array_equal(x.positions,src.positions) and np.array_equal(x.cell,src.cell) and np.array_equal(x.pbc,src.pbc)
  assert x.info['REF_energy']==src.info['REF_energy'] and np.array_equal(x.arrays['REF_forces'],src.arrays['REF_forces'])
  sources[item['path']]=sha(p);lineage.append(item)
 anchors.append({'V16_frame':i,'config_type':x.info['config_type'],'formula':x.get_chemical_formula(),'REF_energy':x.info['REF_energy'],'REF_forces':x.arrays['REF_forces'].tolist(),'cell_A':x.cell.tolist(),'positions_A':x.positions.tolist(),'pbc':x.pbc.tolist(),'persistent_ids':None,'status':'partial_recovery','method_declaration':r['dft_method_declaration'],'parameter_evidence':r['method_parameter_evidence'],'lineage':lineage,'original_calculation_claims':r.get('original_calculation_claims'),'missing_materials':r.get('missing_materials'),'original_import_commit':r['original_import_commit'],'native_vs_free_convention':'unknown; no original getter/log evidence'})
# Search path names in reachable history, without opening blind output content.
cmd=['git','log','--all','--format=','--name-only','--','research/mace-v12-transfer/periodic_interface_v4']
paths=sorted(set(subprocess.check_output(cmd,cwd=ROOT,text=True).splitlines()))
matches=[p for p in paths if any(k in p.lower() for k in ['element_reference','elemental_reference','element_ref','reference_ag','reference_ti','reference_si','reference_c','v10'])]
(OUT/'history_path_search.json').write_text(json.dumps({'command':cmd,'method':'path-name inventory only; reused prior original-content trace; no blind output contents','matching_paths':matches,'scope_limitation':'not a new exhaustive content/hash scan of every Git blob; original missing-artifact findings reused from prior trace'},indent=2)+'\n')
# Compose 43 rows from V16 and previously verified diagnostic formulas only.
diag=WB/'v17_training_pool_provenance_audit/v17_acquisition_summary.json';child=json.loads(diag.read_text());parent=json.loads((WB/'v17_energy_convention_and_parent_hash_audit/audit.json').read_text());bylabel={c['label']:p['formula'] for p in parent['parents'] for c in p['child_inputs']}
from ase.formula import Formula
species=['Ag','C','Si','Ti']
counts=lambda f:[Formula(f).count().get(s,0) for s in species]
N=np.array([counts(x.get_chemical_formula()) for x in a]+[counts(bylabel[x['label']]) for x in child],dtype=float);assert N.shape==(43,4)
sensitivity=[]
for name,removed in [('full',[])]+[(f'remove_{a[i].get_chemical_formula()}',[i]) for i in range(10,14)]+[('remove_all_four',list(range(10,14)))]:
 B=np.delete(N,removed,axis=0);_,sv,vh=np.linalg.svd(B,full_matrices=False);rank=int(np.linalg.matrix_rank(B));sensitivity.append({'scenario':name,'removed_rows':removed,'frames':len(B),'rank':rank,'singular_values':sv.tolist(),'nullspace_basis_species_order_Ag_C_Si_Ti':vh[rank:].tolist()})
(OUT/'evidence_ledger.json').write_text(json.dumps({'sources_sha256':sources,'anchors':anchors},indent=2)+'\n')
(OUT/'rank_sensitivity.json').write_text(json.dumps({'species_order':species,'composition_matrix':N.tolist(),'scenarios':sensitivity},indent=2)+'\n')
with (OUT/'rank_sensitivity.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=['scenario','frames','rank','removed_rows']);w.writeheader();w.writerows({k:r[k] for k in w.fieldnames} for r in sensitivity)
(OUT/'report.md').write_text('''# Four unary reference anchor recovery audit

## Evidence

Ag4, Ti2, Si8 and C8 (V16 rows10–13) retain exact ordered geometry, cell/PBC, REF_energy and REF_forces through the previously documented V11–V14 lineage. Four-row ledger records original import commits, source hashes, declarations and exact missing-artifact claims. No persistent historical atom IDs exist; none were invented. Original calculator inputs, SCF logs, getter convention, GPAW/PAW version and elemental k-point meshes remain unverified. PW-PBE500/Fermi0.1 is a stored declaration, not an original-run certificate. Prior source tracing found original records absent; this focused path-history search found no independently certified exact-geometry replacement. It does not claim exhaustive Git-blob equivalence search.

## Rank and energy gauge

Rank sensitivity for all43 proposed rows is provided in CSV/JSON, including removal of each unary row and all four. For additive per-species energy offsets c, energies change by N c while coordinate forces are unaffected. Any null vector v of N leaves all composition energy equations unchanged under c -> c + t v. A rank-deficient composition matrix cannot identify all four independent offsets. Full rank alone does not certify reference method consistency. Unary solid energies are not isolated-atom chemical potentials; energy/atom here is a label quotient, not a transferability guarantee.

## Minimum safe recovery plan

First request the exact original elemental run records identified in the inherited ledger: fixed-cell input/script, k-point mesh, smearing distribution, XC/cutoff, GPAW/PAW versions, SCF stopping evidence, native/free energies and forces, with input/output hashes. Until obtained, retain explicit provisional provenance flags. If unavailable, A may approve four exact-geometry fixed-cell relabelings; preserve atom order, counts, cell and PBC, generate no synthetic historical IDs, record native extrapolated and free energies separately, and hash all inputs/outputs/logs. PW500/Gamma/Fermi0.1 is a proposed matched target only, not recovered historical settings; small bulk elemental cells require A to review k-point adequacy before adopting Gamma. Do not substitute another lattice, relaxed cell or method archive merely because the element name matches. Four relabelings close these four gaps; they do not certify the ten other partially supported interface frames. One unary row may suffice to resolve a rank-one gauge deficiency, but cannot resolve all four provenance uncertainties.

No DFT, inference, training, dataset or role edits were performed. No development-validation or withheld output content was accessed. Reproduce with audit.py in the existing GPAW Python environment. The history search is deliberately path-only and the full previous source trace is reused rather than repeated.
'''+ '\nRank results: '+json.dumps([{k:r[k] for k in ['scenario','rank']} for r in sensitivity])+'\n')
(OUT/'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='SHA256SUMS.txt'))
print([(r['scenario'],r['rank']) for r in sensitivity]);print('path matches',len(matches))
