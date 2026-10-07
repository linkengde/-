import json,hashlib,csv
from pathlib import Path
OUT=Path(__file__).resolve().parent;WB=OUT.parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
paths=['v17_certified_core_rank_recovery_analysis/rank_analysis.json','v17_training_pool_provenance_audit/ledger.json','v17_split_leakage_audit/audit.json','v17_element_reference_anchor_recovery_audit/rank_sensitivity.json']
for directory in {p.split('/')[0] for p in paths}:
 for line in (WB/directory/'SHA256SUMS.txt').read_text().splitlines():
  h,name=line.split('  ');assert sha(WB/directory/name)==h
rank=json.loads((WB/paths[0]).read_text());assert rank['core_count']==29 and rank['core_rank']==3 and rank['minimum_additional_rows']==1
solutions=[r['frames'][0] for r in rank['all_minimum_solutions']]
options=[{'option':'A_minimum_certified_extension','train_frames':30,'composition_rank':4,'new_DFT_jobs':1,'rank_lifting_frame_choices':solutions,'partial_rows_after_successful_relabel':0,'certainty':'conditional: new exact-geometry archive must pass; existing29 direct','energy_target':'retain native extrapolated target and preserve free_energy separately','limitation':'minimal rank restoration; no guarantee of force coverage; unary Gamma adequacy needs A review'}, {'option':'B_minimum_provisional_extension','train_frames':30,'composition_rank':4,'new_DFT_jobs':0,'rank_lifting_frame_choices':solutions,'partial_rows_after_successful_relabel':1,'certainty':'29 direct plus one explicitly partial','energy_target':'historical getter and some settings unverified for added row','limitation':'quick bookkeeping path, cannot claim fully archive-certified common target'}, {'option':'B_full_provisional_inheritance','train_frames':43,'composition_rank':4,'new_DFT_jobs':0,'rank_lifting_frame_choices':list(range(14)),'partial_rows_after_successful_relabel':14,'certainty':'29 direct plus14 partial','energy_target':'declared methods differ; Gamma versus1x4x2, four unary getter unknown','limitation':'preserves more coverage, retains provenance and numerical-method uncertainty'}, {'option':'C_relabel_all_partial','train_frames':43,'composition_rank':4,'new_DFT_jobs':14,'rank_lifting_frame_choices':list(range(14)),'partial_rows_after_successful_relabel':0,'certainty':'conditional after14 exact matched archives pass','energy_target':'A-approved common PBE PW500 Gamma Fermi0.1 native target; bulk kpoint adequacy must be reviewed','limitation':'closes inherited provenance gaps but not force/model quality; different approved meshes must be recorded if needed'}]
(OUT/'options.json').write_text(json.dumps({'source_hashes':{p:sha(WB/p) for p in paths},'options':options},indent=2)+'\n')
with (OUT/'options.csv').open('w') as f:
 keys=['option','train_frames','composition_rank','new_DFT_jobs','partial_rows_after_successful_relabel','certainty','energy_target','limitation'];w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows({k:r[k] for k in keys} for r in options)
(OUT/'report.md').write_text('''# V17 dataset integration recipe for A review

## Fastest archive-supported rank restoration

The 29 directly supported training rows have exact rank3. One exact matched relabel can provide a 30-frame rank4 pool. All13 eligible V16 frame choices are 0–8 and10–13; frame9 does not lift rank. Full identities, formulas and declared settings are in the rank report. Choose based on scientific coverage and cost, not rank alone; no evidence here selects a best force-learning candidate. Recovering an original valid archive is preferable to a new calculation if available. No jobs are authorized or launched by this recipe.

## Alternatives

A: core29 plus one successful exact-geometry relabel ->30 rows/rank4/one job/zero partial rows. B: core29 plus one provisional partial ->30 rows/rank4/zero jobs/one unresolved row, or all inherited ->43 rows/rank4/14 partial. C: relabel all14 partial exact geometries ->43 rows/rank4/14 jobs; relabeling only the minimum subset has the same bookkeeping outcome as A. Intermediate added rows must each justify coverage or provenance benefit. All counts assume replacement of historical labels for selected geometries, never double counting old and new labels.

## Builder acceptance gates

1. Freeze A-approved training membership and hashes. Keep the29 directly supported row identities and exact copied REF fields; map any replacement to its original frame and new archive without pretending it recovers original-run settings.
2. Require each new archive input geometry/order/species/cell/PBC identity, converged SCF, finite labels, original log and settings evidence, energy-field identity and exact force copy. Preserve historical absent IDs; do not invent them.
3. Use a declared native extrapolated target; retain free energy separately. Common target does not imply forces are exact derivatives of extrapolated energy at finite smearing. Prior paired-response evidence motivates preserving both and reviewing derivative consistency; it does not prove archives are wrong. Unary bulk Gamma adequacy must be reviewed before adopting interface numerical settings. Distinguish XC agreement from kpoint/smearing/version agreement.
4. Recompute exact composition rank, row counts, all source hashes, replacement/no-duplicate ledger and role isolation. Keep development and withheld labels outside training and model selection; compare permitted inputs under frozen algorithms before integration. Do not call shared-parent tests morphology-independent.
5. Publish explicit provenance tiers and uncertainty, then require A review before production build/training. Algebraic rank and provenance certainty are acceptance checks, not model force gates or permission for MD/TTM.

## Scope

Only previously published training audit metadata was reused; no development/withheld candidates or outputs, scores, model inference, DFT, dataset or role edits. Source artifacts were SHA256 verified. Reproduce with analyze.py. Options JSON/CSV contains exact counts and prerequisites.
''')
(OUT/'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='SHA256SUMS.txt'))
print('Verified sources; four integration options produced')
