"""Training-only provenance reconciliation. Never traverses unseen output roots."""
import csv,json,re,hashlib
from pathlib import Path
import numpy as np
from ase.io import read
from inventory import ROOT,OUT,BASE,DATA,sha,run
run()
m=json.loads((DATA/'dataset_manifest.json').read_text()); frames=read(DATA/'train.extxyz',':')
rbase=BASE/'coordination/reports/window-b'
prior=json.loads((rbase/'v17_energy_convention_and_parent_hash_audit/audit.json').read_text())
ledger=[]; extra=[]
for i,a in enumerate(frames):
 row={'frame':i,'config_type':a.info['config_type'],'atoms':len(a),'formula':a.get_chemical_formula()}
 prov=next((x for x in m['per_frame_provenance'] if x['split']=='train' and x['frame']==i),None)
 if prov:
  src=read(ROOT/prov['source_path'],':')[prov['source_frame']]
  assert np.array_equal(a.numbers,src.numbers) and np.array_equal(a.positions,src.positions)
  assert np.array_equal(a.cell,src.cell) and np.array_equal(a.pbc,src.pbc)
  assert a.info['REF_energy']==src.info.get('REF_energy',src.info.get('PW_PBE_energy_eV')) and np.array_equal(a.arrays['REF_forces'],src.arrays.get('REF_forces',src.arrays.get('PW_PBE_forces')))
  row['inherited_source']=prov;row['current_row_copy_check']='PASS'
 if i<14:
  row.update(evidence_status='partially_supported',method_evidence=prov['parameter_evidence'],original_run_evidence=prov['original_run_evidence'],energy_convention=prov['energy_convention_evidence'],original_archive_REF_identity='unknown')
 elif 23<=i<=25:
  cached=next(x for x in prior['parents'] if x['config_type']==a.info['config_type'])
  row.update(evidence_status='directly_archive_verified',reused_parent_evidence=cached,energy_convention='native extrapolated; original log verified',method=cached['method'])
 else:
  if 14<=i<=17:
   label=a.info['config_type'].split('_periodic_')[0];folder=BASE/'periodic_interface_v4/pbe_interface_energy_additions_v12/calculations'/label
  elif 18<=i<=20:
   label=a.info['config_type'].split('_periodic_')[0];folder=BASE/'periodic_interface_v4/pbe_interface_energy_additions_v12/calculations/validation_distance_scans'/label
  elif 21<=i<=22:
   label=a.info['config_type'].split('_periodic_')[0];folder=BASE/'periodic_interface_v4/pbe_interface_v13_holdouts/calculations'/label
  elif 26<=i<=28:
   folder=(ROOT/prov['source_path']).parent;label=folder.name
  else:
   addition=m['V16_additions'][i-29];folder=ROOT/addition['source'];label=addition['label']
  output=folder/(label+'_PW_PBE.extxyz');summary=folder/'summary.json';log=folder/'gpaw.log'
  d=json.loads(summary.read_text());src=read(output)
  assert np.array_equal(a.numbers,src.numbers) and np.array_equal(a.positions,src.positions)
  assert np.array_equal(a.cell,src.cell) and np.array_equal(a.pbc,src.pbc)
  assert a.info['REF_energy']==src.info['PW_PBE_energy_eV']==d['energy_eV_cell']
  assert np.array_equal(a.arrays['REF_forces'],src.arrays['PW_PBE_forces'])
  assert d['scf_converged']
  matches=re.findall(r'Extrapolated:\s*([-+0-9.eE]+)',log.read_text());assert matches
  assert abs(float(matches[-1])-d['energy_eV_cell'])<=0.00000051
  hashes={str(p.relative_to(ROOT)):sha(p) for p in [output,summary,log]}
  inp=d.get('input',d.get('source'))
  if isinstance(inp,str) and '/research/' in inp: inp='research/'+inp.split('/research/',1)[1]
  if isinstance(inp,str) and (ROOT/inp).is_file(): 
   hashes[inp]=sha(ROOT/inp)
   pin=d.get('source_sha256',d.get('input_sha256'));assert pin is None or hashes[inp]==pin
  for p,h in hashes.items():
   expected=m['V16_sources_sha256'].get(p)
   if expected: assert expected==h
  v=folder/'verification.json'
  if v.is_file():
   vd=json.loads(v.read_text());hashes[str(v.relative_to(ROOT))]=sha(v);row['existing_verification']=vd
  row.update(evidence_status='directly_archive_verified',label=label,method=d.get('method'),scf_converged=True,mpi_ranks=d.get('mpi_ranks',d.get('execution',{}).get('mpi_ranks')),energy_convention='native extrapolated; original log verified',REF_energy_exact_source=True,REF_forces_exact_source=True,source_hashes=hashes,free_energy_eV_cell=d.get('free_energy_eV_cell'),log_extrapolated_eV_rounded=float(matches[-1]))
  extra.extend({'path':p,'sha256':h} for p,h in hashes.items())
 ledger.append(row)
children=[x for p in prior['parents'] for x in p['child_inputs']];assert len(children)==8 and all(x['archive_verification']=='PASS' for x in children)
composition=lambda rows: int(np.linalg.matrix_rank([[a.get_chemical_symbols().count(s) for s in ['Ag','C','Si','Ti']] for a in rows]))
counts={'V16_direct':sum(x['evidence_status']=='directly_archive_verified' for x in ledger),'V16_partial':sum(x['evidence_status']=='partially_supported' for x in ledger),'V17_diagnostics_direct':len(children),'full_composition_rank':composition(frames),'direct_only_composition_rank':composition([a for a,r in zip(frames,ledger) if r['evidence_status']=='directly_archive_verified'])}
(OUT/'ledger.json').write_text(json.dumps({'counts':counts,'rows':ledger,'additional_evidence_inventory':extra},indent=2)+'\n')
(OUT/'v17_acquisition_summary.json').write_text(json.dumps(children,indent=2)+'\n')
with (OUT/'v17_acquisition_summary.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=['label','interface','owner','archive_verification','energy_convention','native_energy_eV_cell','free_energy_eV_cell','input_sha256']);w.writeheader();w.writerows({k:r.get(k) for k in w.fieldnames} for r in children)
with (OUT/'ledger.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=['frame','config_type','atoms','formula','evidence_status','energy_convention']);w.writeheader();w.writerows({k:(json.dumps(r[k]) if isinstance(r.get(k),dict) else r.get(k)) for k in w.fieldnames} for r in ledger)
(OUT/'report.md').write_text(f'''# V17 training pool provenance audit

## Result

{counts['V16_direct']} of 35 V16 rows are directly archive verified; {counts['V16_partial']} remain partially supported. All eight V17 signed training diagnostics were already archive verified and retain native and free energies separately. Therefore the proposed 43-row pool has {counts['V16_direct']+8} directly supported rows and 14 partially supported rows. It cannot be certified as a completely original-run-proven homogeneous target.

## Evidence and method

Current V16 rows were compared exactly against their inherited source rows. Newly reconciled archives (rows 14–22) and residual/V16 archives were compared for ordered species, coordinates, cell/PBC, exact REF energy and force copies, converged summaries and original-log extrapolated energy (six-decimal log tolerance 5.1e-7 eV). Three parent archives and eight signed diagnostics reuse the completed hash-pinned parent audit. No development-validation or withheld output was read. Source declarations are retained separately from original evidence in ledger.json.

Directly supported interface labels use native extrapolated GPAW PBE energies, PW500/Gamma/0.1 eV smearing where recorded. Rows 0–13 retain their original declaration evidence and precise missing-source paths/hashes. Rows 0–9 have serialized native aliases, which do not independently prove the original getter; four elemental references have unproven energy convention. Older gap declarations use 1x4x2 rather than Gamma: common XC alone does not establish identical numerical settings. No evidence establishes that a free-energy label was substituted into the directly verified native target.

## Integration recommendation

Do not describe all 43 labels as fully certified. Preserve the 14 uncertain provenance flags and method differences; seek original logs/inputs/outputs or authorize matched relabeling before claiming complete consistency. Do not silently replace native energies with free energies, and keep both fields for future diagnostics. Full composition rank is {counts['full_composition_rank']}; archive-supported-only rank is {counts['direct_only_composition_rank']}, so removing uncertain elemental anchors is not a neutral filter. A must review target convention and reference-anchor policy before construction. Provenance verification does not prove force accuracy, model quality, or thermodynamic transferability.

## Reproduction

Run inventory.py, then audit_pool.py with the existing GPAW Python environment. Both scripts operate only on the listed training evidence. SHA256SUMS.txt covers artifacts, not itself.
''')
paths=sorted(p for p in OUT.iterdir() if p.is_file() and p.name!='SHA256SUMS.txt')
(OUT/'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in paths))
print(counts)
