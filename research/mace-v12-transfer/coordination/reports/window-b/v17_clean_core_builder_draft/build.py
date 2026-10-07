import json,hashlib,importlib.util,io,re,csv
from pathlib import Path
import numpy as np
from ase.io import read,write
ROOT=Path(__file__).resolve().parents[6];OUT=Path(__file__).resolve().parent;WB=OUT.parent;BASE=ROOT/'research/mace-v12-transfer/periodic_interface_v4'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
spec=importlib.util.spec_from_file_location('geom',WB/'v17_split_geometry_design/final_set/generate_split_candidates.py');g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)
spec=importlib.util.spec_from_file_location('pf',BASE/'mace_periodic_v17_interface_energy/preflight_v17.py');pf=importlib.util.module_from_spec(spec);spec.loader.exec_module(pf)
ledger_path=WB/'v17_training_pool_provenance_audit/ledger.json';ledger=json.loads(ledger_path.read_text());src=BASE/'mace_periodic_v16_interface_energy/data/train.extxyz';train=read(src,':')
manifest=json.loads((src.parent/'dataset_manifest.json').read_text());assert sha(src)==manifest['output_sha256']['train']
frames=[];rows=[];pins={str(src.relative_to(ROOT)):sha(src),str(ledger_path.relative_to(ROOT)):sha(ledger_path)}
for r in ledger['rows']:
 if r['evidence_status']=='directly_archive_verified':
  a=train[r['frame']].copy();a.calc=None;a.info['proposal_role']='train';a.info['proposal_provenance']='directly_archive_verified';frames.append(a);rows.append({'frame':len(rows),'source':'V16','source_frame':r['frame'],'config':a.info['config_type'],'formula':a.get_chemical_formula(),'evidence':'directly_archive_verified','source_hash':sha(src)})
assert len(frames)==21
for directory in ['pbe_interface_v17_main_targeted_acquisition','pbe_interface_v17_parallel_targeted_acquisition']:
 mp=BASE/directory/'input_manifest.json'
 for r in json.loads(mp.read_text())['records']:
  folder=mp.parent/'calculations'/r['label'];p=folder/(r['label']+'_PW_PBE.extxyz');v=json.loads((folder/'verification.json').read_text());assert v['status']=='PASS' and all(v['checks'].values()) and sha(p)==v['sha256'][p.name]
  a=read(p);a.calc=None;a.info['REF_energy']=a.info['PW_PBE_energy_eV'];a.arrays['REF_forces']=a.arrays['PW_PBE_forces'].copy();a.info['proposal_role']='train';a.info['proposal_provenance']='directly_archive_verified';a.info['source_method']='GPAW26.7 PW-PBE500 Gamma Fermi0.1';frames.append(a);pins[str(p.relative_to(ROOT))]=sha(p);rows.append({'frame':len(rows),'source':str(p.relative_to(ROOT)),'config':r['label'],'formula':a.get_chemical_formula(),'parent':r['source_parent'],'signed_design':r['design'],'source_hash':sha(p),'evidence':'directly_archive_verified'})
pack=BASE/'pbe_interface_v17_rank_recovery_parallel';rec=json.loads((pack/'input_manifest.json').read_text())['records'][0];folder=pack/'calculations'/rec['label'];p=folder/(rec['label']+'_PW_PBE.extxyz');v=json.loads((folder/'verification.json').read_text());assert v['status']=='PASS' and all(v['checks'].values()) and sha(p)==v['sha256'][p.name];a=read(p);a.calc=None;original=train[0];assert np.array_equal(a.positions,original.positions) and np.array_equal(a.cell,original.cell) and np.array_equal(a.numbers,original.numbers) and np.array_equal(a.pbc,original.pbc);a.info['REF_energy']=a.info['PW_PBE_energy_eV'];a.arrays['REF_forces']=a.arrays['PW_PBE_forces'].copy();a.info['source_method']='GPAW26.7 PW-PBE500 Gamma Fermi0.1';a.info['proposal_role']='train';a.info['proposal_provenance']='directly_archive_verified';frames.append(a);pins[str(p.relative_to(ROOT))]=sha(p);rows.append({'frame':29,'source':str(p.relative_to(ROOT)),'config':rec['label'],'formula':a.get_chemical_formula(),'source_frame':0,'source_hash':sha(p),'evidence':'directly_archive_verified','caveat':rec['caveat']});assert len(frames)==30
# Serialize using ASE schema, then restore all coordinate/force columns at17 digits.
blocks=[]
for a in frames:
 buf=io.StringIO();write(buf,a,format='extxyz',write_results=False);lines=buf.getvalue().splitlines();lines[1]=re.sub(r'Lattice="[^"]+"','Lattice="'+' '.join(format(v,'.17g') for v in a.cell.array.ravel())+'"',lines[1]);schema=re.search(r'Properties=([^ ]+)',lines[1]).group(1).split(':');offset=0
 fields={}
 for j in range(0,len(schema),3):fields[schema[j]]=(offset,int(schema[j+2]));offset+=int(schema[j+2])
 for i in range(len(a)):
  values=lines[i+2].split()
  for key in ['pos','REF_forces','PW_PBE_forces']:
   if key in fields:
    q,n=fields[key];array=a.positions if key=='pos' else a.arrays[key];values[q:q+n]=[format(v,'.17g') for v in array[i]]
  lines[i+2]=' '.join(values)
 blocks.append('\n'.join(lines)+'\n')
output=OUT/'train_candidate.extxyz';output.write_text(''.join(blocks));loaded=read(output,':')
for a,b in zip(frames,loaded):
 assert np.array_equal(a.positions,b.positions) and np.array_equal(a.cell,b.cell) and np.array_equal(a.numbers,b.numbers) and np.array_equal(a.pbc,b.pbc) and np.array_equal(a.arrays['REF_forces'],b.arrays['REF_forces']) and a.info['REF_energy']==b.info['REF_energy']
 assert np.isfinite(b.info['REF_energy']) and np.isfinite(b.arrays['REF_forces']).all()
N=np.array([[a.get_chemical_symbols().count(s) for s in ['Ag','C','Si','Ti']] for a in loaded]);assert pf.exact_rank(N)==4
# Parse only development geometry columns, use withheld role digests only.
dev=g.extxyz_geometry_frames(BASE/'mace_periodic_v17_interface_energy/data/valid.extxyz');roles=json.loads((WB/'v17_split_geometry_design/final_set/role_manifest.json').read_text());withheld={r['geometry_sha256'] for r in roles['candidates'] if r['proposed_role']=='withheld_test'};hits=[]
for i,a in enumerate(loaded):
 assert g.geometry_sha(a) not in withheld
 for j,b in enumerate(dev):
  match=g.compare_geometry(a,b)
  if match and (match['exact_duplicate'] or match['near_duplicate']):hits.append({'train':i,'dev':j,'match':match})
assert not hits
result={'status':'PROPOSAL_ONLY_NOT_INTEGRATED','frames':30,'composition_rank_exact':4,'sources_sha256':pins,'rows':rows,'output_sha256':sha(output),'development_geometry_overlap':hits,'withheld_role_hash_overlap':False,'energy_convention':'native extrapolated; free energy retained in source archives','frame0_kpoint_convergence':'UNRESOLVED; old declaration1x4x2 versus newGamma','limits':'archive supported does not mean uniform kpoints, force accuracy or transferability'}
(OUT/'manifest.json').write_text(json.dumps(result,indent=2)+'\n')
(OUT/'report.md').write_text('# Isolated clean-core builder proposal\n\nPROPOSAL_ONLY_NOT_INTEGRATED.30 frames=21 directly verified V16 +8 signed V17 diagnostics +new exact frame0 relabel. Exact rational composition rank4; source/output hashes, finite labels, ordered symbols, exact coordinates/cell/PBC and serialized label identity pass. Intended signed pairs retain parent/displacement lineage. No inherited partial-source frame is included; relabel replaces historical frame0 only within this proposal. Eight pairs are training diagnostics, not independent validation.\n\nGamma frame0 differs from the original1x4x2 declaration; convergence along3.08/5.33 A periodic directions remains unresolved. Historical IDs absent for frame0 remain absent. Direct archive support is provenance evidence, not a common-method or force-quality guarantee. Development geometry-only exact/near comparison passes under the published thresholds; withheld geometry role digests do not match, with no withheld inputs/outputs opened. No development label values/scores were read by this builder.\n\nA must review before integration or training. Existing official V17 data is untouched. Reproduce with build.py in existing GPAW Python environment; the script writes only this isolated report directory. No DFT/model execution.\n')
(OUT/'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='SHA256SUMS.txt'));print('proposal30 exact rank4; serialization/source/geometry checks PASS')
