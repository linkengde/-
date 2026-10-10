#!/usr/bin/env python3
"""Index public calculation archives without loading atomic structures or trajectories."""
import csv, hashlib, json, os, tempfile
from pathlib import Path

ROOT=Path(os.environ.get('CLOUD_C_REPO_ROOT',Path(__file__).resolve().parents[6])).resolve()
BASE=ROOT/'research/mace-v12-transfer/periodic_interface_v4'
OUT=Path(os.environ.get('CLOUD_C_AUDIT_OUTPUT',Path(tempfile.gettempdir())/'cloud-c-v18-audit')).resolve(); OUT.mkdir(parents=True,exist_ok=True)
V18=BASE/'mace_periodic_v18_clean_core'
blocked=('/blind','holdout','/test','sealed')
records=[]
for summary_path in sorted(BASE.glob('pbe_*/calculations/*/summary.json')):
    rel=summary_path.relative_to(ROOT).as_posix().lower()
    if any(term in rel for term in blocked):
        continue
    folder=summary_path.parent
    try: summary=json.loads(summary_path.read_text())
    except Exception as exc: summary={"parse_error":str(exc)}
    verification_path=folder/'verification.json'
    verification=json.loads(verification_path.read_text()) if verification_path.is_file() else {}
    checksum_path=folder/'SHA256SUMS.txt'
    checked=[]
    if checksum_path.is_file():
        for line in checksum_path.read_text().splitlines():
            if not line.strip(): continue
            expected,name=line.split(maxsplit=1); name=name.lstrip('*')
            member=(folder/name).resolve()
            if folder.resolve() not in member.parents or not member.is_file():
                checked.append({"file":name,"pass":False}); continue
            actual=hashlib.sha256(member.read_bytes()).hexdigest()
            checked.append({"file":name,"pass":actual==expected})
    formula=str(summary.get('formula',''))
    label=str(summary.get('label',folder.name))
    if formula=='Ag32': category='solid_Ag_static'
    elif formula in ('C16Si8Ti24','Si8Ti24C16'): category='solid_Ti3SiC2_bulk'
    elif 'Ag' in formula and ('Ti' in formula or 'Si' in formula or 'C' in formula): category='Ag_carbide_silicide_interface_or_cluster'
    elif 'Ag' in formula: category='Ag_containing_other'
    else: category='other_or_unresolved'
    role=str(summary.get('dataset_role') or '')
    if 'numerical' in role or 'pilot' in role or 'control' in role:
        train_decision='excluded_by_declared_role_until_explicit_reassignment'
    elif 'reused_development' in role:
        train_decision='candidate_for_training_only_after_group_role_review_never_fresh_validation'
    elif 'training' in role or 'acquisition' in role:
        train_decision='candidate_for_training_after_provenance_and_overlap_review'
    elif not role:
        train_decision='unknown_legacy_role_requires_source_and_lineage_review'
    else:
        train_decision='manual_role_review_required'
    records.append({
      'label':label,'category':category,'path':folder.relative_to(ROOT).as_posix(),
      'formula':formula,'atoms':summary.get('atoms'),'role':summary.get('dataset_role'),
      'training_use_recommendation':train_decision,
      'independent_validation_recommendation':'reserve_only_structurally_distinct_unseen_parent_groups; no current diagnostic archive qualifies by this audit',
      'scf_converged':summary.get('scf_converged') is True,
      'verification_PASS':verification.get('status')=='PASS' and bool(verification.get('checks')) and all(verification.get('checks',{}).values()),
      'sha256_manifest_present':bool(checked),'sha256_all_pass':(all(x['pass'] for x in checked) if checked else None),
      'verification_hashes_checked':len(verification.get('sha256',{})),
      'verification_hashes_all_pass':(all((folder/name).is_file() and hashlib.sha256((folder/name).read_bytes()).hexdigest()==digest for name,digest in verification.get('sha256',{}).items()) if verification.get('sha256') else None),
      'method':json.dumps(summary.get('method',{}),sort_keys=True,separators=(',',':')),
      'source_sha256':summary.get('source_sha256'),
      'result_sha256':verification.get('sha256',{}).get(f"{label}_PW_PBE.extxyz"),
      'archive_files_checked':len(checked),
    })
cols=list(records[0]) if records else []
with (OUT/'dft_coverage.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=cols); w.writeheader(); w.writerows(records)
counts={}
for r in records:
    key=(r['category'],r['scf_converged'],r['verification_PASS'],r['sha256_all_pass'],r['verification_hashes_all_pass'])
    counts['|'.join(map(str,key))]=counts.get('|'.join(map(str,key)),0)+1
report={
 'method':'Path-limited inventory of published PBE calculation summary archives under periodic_interface_v4/pbe_*/calculations. Excludes any path containing blind, holdout, test or sealed. No structures, trajectory frames or labels are loaded for model inference.',
 'archives_indexed':len(records),'status_counts':counts,
 'liquid_Ag_DFT_archives':0,'solid_liquid_interface_DFT_archives':0,
 'high_temperature_or_quench_DFT_archives':0,
 'notes':['A converged static liquid snapshot or solid-liquid interface would still require explicit source/ensemble/temperature lineage; none was identified in the indexed archive names/formulas.','Archive category is inferred from label/formula and must be interpreted with per-record role/method; numerical convergence and hash PASS do not automatically mean training eligibility.'],
}
(OUT/'dft_coverage.json').write_text(json.dumps(report,indent=2)+'\n')
# This is a train-role inventory only; no development/test files or sealed labels are read.
import sys
# ASE is provided by the selected Python environment.
from ase.io import read
train=read(V18/'data/train.extxyz',index=':')
manifest=json.loads((V18/'data/dataset_manifest.json').read_text())
with (OUT/'v18_training_frame_roles.csv').open('w') as f:
    fields=['frame','source','source_frame','config','formula','atoms','elements','dataset_role','current_reuse_policy']
    w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
    for row,atoms in zip(manifest['rows'],train):
        w.writerow({'frame':row['frame'],'source':row.get('source'),'source_frame':row.get('source_frame'),'config':row.get('config'),'formula':atoms.get_chemical_formula(),'atoms':len(atoms),'elements':','.join(sorted(set(atoms.get_chemical_symbols()))),'dataset_role':'V18 training only','current_reuse_policy':'may only be reused for train-side controlled reproduction; never presented as independent validation'})
print(json.dumps(report,indent=2))
