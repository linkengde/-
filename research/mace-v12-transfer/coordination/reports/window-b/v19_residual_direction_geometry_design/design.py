import json,csv,hashlib,importlib.util
from pathlib import Path
import numpy as np
from ase.io import read
ROOT=Path(__file__).resolve().parents[6];OUT=Path(__file__).resolve().parent;WB=OUT.parent;BASE=ROOT/'research/mace-v12-transfer/periodic_interface_v4';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
spec=importlib.util.spec_from_file_location('g',WB/'v17_split_geometry_design/final_set/generate_split_candidates.py');g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)
parents=g.extxyz_geometry_frames(BASE/'mace_periodic_v17_interface_energy/data/valid.extxyz');vecpath=WB/'v17_frozen_dev_vector_diagnosis/per_atom.csv';vectors=list(csv.DictReader(vecpath.open()));roles=json.loads((WB/'v17_split_geometry_design/final_set/role_manifest.json').read_text());sealed={r['geometry_sha256'] for r in roles['candidates'] if r['proposed_role']=='withheld_test'}
inventory=[];sources={str(vecpath.relative_to(ROOT)):sha(vecpath)}
paths=[BASE/f'{name}/data/train.extxyz' for name in ['mace_periodic_v14_interface_energy','mace_periodic_v15_interface_energy','mace_periodic_v16_interface_energy','mace_periodic_v17_interface_energy','mace_periodic_v18_clean_core']]
paths+=[BASE/'mace_periodic_v17_interface_energy/data/valid.extxyz']
for folder in BASE.glob('pbe_interface_v*_acquisition'):
 if 'v17' not in folder.name or folder.name in ['pbe_interface_v17_main_targeted_acquisition','pbe_interface_v17_parallel_targeted_acquisition']:
  paths+=list((folder/'inputs').glob('*.extxyz'))
for name in ['v15_wide_span_sampling_design','v15_registry_candidates','v16_unseen_validation_geometry_design','targeted_transverse_framework_design']:
 paths+=list((WB/name).glob('**/geometries/*.extxyz'))
for p in sorted(set(paths)):
 if p.is_file():
  sources[str(p.relative_to(ROOT))]=sha(p)
  for i,a in enumerate(g.extxyz_geometry_frames(p)):inventory.append((str(p.relative_to(ROOT))+f'#{i}',a))
(OUT/'geometries').mkdir(exist_ok=True);records=[]
for interface,parent in zip(['AgC','AgSi','AgTi'],parents):
 ids=parent.arrays['lammps_id'];mapping={int(v):i for i,v in enumerate(ids)};central=np.flatnonzero(parent.arrays['central_pair']);ag=next(i for i in central if parent[i].symbol=='Ag');x=next(i for i in central if parent[i].symbol!='Ag');axis=parent.positions[x]-parent.positions[ag];axis/=np.linalg.norm(axis)
 if interface=='AgC':
  row=next(r for r in vectors if r['label'].startswith(interface+'_') and int(r['id'])==604600);e=np.array([float(row['error_'+d]) for d in 'xyz']);mode=e-np.dot(e,axis)*axis;mode/=np.linalg.norm(mode);moves={604600:mode}
 elif interface=='AgTi':moves={479027:axis,480805:-axis}
 else:moves={366934:axis,685172:-axis}
 parent_min=g.pair_minima(parent)
 for amplitude in [.03,.06,.10]:
  for sign in [-1,1]:
   a=parent.copy();displacements=[]
   for atom_id,direction in moves.items():
    delta=sign*amplitude*direction;a.positions[mapping[atom_id]]+=delta;displacements.append({'atom_id':atom_id,'species':a[mapping[atom_id]].symbol,'delta_xyz_A':delta.tolist()})
   label=f'{interface}_residual_mode_{"m" if sign<0 else "p"}_{int(amplitude*100):02d}_v19_proposal';p=OUT/'geometries'/f'{label}.extxyz'
   header='Lattice="'+' '.join(format(v,'.17g') for v in a.cell.array.ravel())+'" Properties=species:S:1:pos:R:3:lammps_id:I:1:central_pair:I:1 pbc="T T T"'
   p.write_text(str(len(a))+'\n'+header+'\n'+''.join(s+' '+' '.join(format(v,'.17g') for v in xyz)+f' {int(atom_id)} {int(c)}\n' for s,xyz,atom_id,c in zip(a.get_chemical_symbols(),a.positions,ids,a.arrays['central_pair'])))
   b=g.extxyz_geometry_frames(p)[0];assert np.array_equal(a.positions,b.positions) and np.array_equal(a.cell,b.cell) and np.array_equal(a.arrays['lammps_id'],b.arrays['lammps_id']);assert g.geometry_sha(b) not in sealed
   minima=g.pair_minima(b);flags=[pair for pair,r in minima.items() if r['distance_A']<parent_min[pair]['distance_A']-.05]
   exact=[];near=[]
   for name,c in inventory:
    match=g.compare_geometry(b,c)
    if match and match['exact_duplicate']:exact.append(name)
    elif match and match['near_duplicate']:near.append(name)
   # all pair contacts: not just global minima, check moved contacts against parent.
   pd=parent.get_all_distances(mic=True);bd=b.get_all_distances(mic=True);contact_flags=[]
   for atom_id in moves:
    i=mapping[atom_id]
    for j in range(len(a)):
     if i!=j and bd[i,j]<pd[i,j]-.15 and bd[i,j]<parent_min['-'.join(sorted([a[i].symbol,a[j].symbol]))]['distance_A']-.05:contact_flags.append([atom_id,int(ids[j])])
   fracdelta=(b.positions-parent.positions)@np.linalg.inv(parent.cell.array);boundary=bool(np.max(np.abs(fracdelta))<.5)
   records.append({'candidate_id':label,'interface':interface,'amplitude_A':amplitude,'sign':sign,'path':str(p.relative_to(ROOT)),'sha256':sha(p),'geometry_sha256':g.geometry_sha(b),'moved_atoms':displacements,'normalization':'max moved atom norm = nominal amplitude; opposite two-atom modes relative change=2*amplitude','parent_role':'already-scored V17 development reused as acquisition family only','parent_geometry_sha256':g.geometry_sha(parent),'species_minima_A':minima,'parent_relative_contact_flags':flags,'contact_pair_flags':contact_flags,'small_periodic_displacement_check':boundary,'exact_hits':exact,'near_hits':near,'screening_pass':not flags and not contact_flags and boundary and not exact,'near_overlap_sealed_evaluated':False})
assert len(records)==18
recommend=[]
for interface in ['AgC','AgSi','AgTi']:
 for amplitude in [.06,.03,.10]:
  pair=[r for r in records if r['interface']==interface and r['amplitude_A']==amplitude]
  if all(r['screening_pass'] for r in pair):recommend.extend(r['candidate_id'] for r in pair);break
proposal_pairs=[]
for i,r in enumerate(records):
 for t in records[i+1:]:
  match=g.compare_geometry(g.extxyz_geometry_frames(ROOT/r['path'])[0],g.extxyz_geometry_frames(ROOT/t['path'])[0])
  if match and (match['exact_duplicate'] or match['near_duplicate']):proposal_pairs.append({'left':r['candidate_id'],'right':t['candidate_id'],'comparison':match})
assert not any(x['comparison']['exact_duplicate'] for x in proposal_pairs)
result={'proposal_pair_overlap':proposal_pairs,'task':'v19_residual_direction_geometry_design','candidate_count':18,'recommended':recommend,'sources_sha256':sources,'records':records,'sealed_check':'exact role digests only; near sealed geometry not inspected','near_policy':'record intentional local training-family correlations, never claim independent validation; exact duplicates rejected','contact_policy':'species minimum compared to parent minus0.05A; extra moved-contact guard; geometric triage not physical certification'}
(OUT/'manifest.json').write_text(json.dumps(result,indent=2)+'\n');(OUT/'report.md').write_text('# V19 label-free residual geometry proposals\n\n18 proposals: signed amplitudes0.03/0.06/0.10A; one normalized mode per interface. AgC Ti604600 follows its observed transverse residual direction. AgTi Ti479027 and frameworkC480805 move oppositely along marked pair axis to probe differential framework response. AgSi Si366934 and Ag685172 move oppositely along axis. Two-atom relative displacement is twice nominal amplitude; all vectors/IDs/source hashes in manifest. IDs/order/cell/PBC preserved; inputs contain geometry and pair marker only.\n\nRecommended signed pairs: '+str(recommend)+'\n\nPrefer0.06A if both contacts pass, else0.03A, then0.10A; larger displacement improves finite-response signal but perturbs more contacts, and is reserve not an automatic compute request. Species minima use parent-specific0.05A allowance rather than all-element1.75A threshold. All passing near hits are reported: these deliberately probe a reused development family and are highly correlated with existing geometries. They are training-acquisition proposals, not blind or morphology-independent tests. Screening compares permitted historical train/dev/acquisition/prior candidates only; sealed comparison uses exact digests, so near sealed overlap remains unknown.\n\nA must review contact/near-correlation tradeoffs and explicitly assign DFT owners before labeling. Do not label18 redundant amplitudes automatically. No inference/DFT/training/MD/TTM or official edits performed.\n');(OUT/'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.relative_to(OUT)}\n' for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='SHA256SUMS.txt'));print('proposals18; passing',sum(r['screening_pass'] for r in records),'recommended',recommend)
