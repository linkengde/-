"""Inventory permitted training/proposal geometry only; no sealed inputs or labels."""
import collections,hashlib,json,subprocess
from pathlib import Path
import numpy as np
from ase.io import read
from ase.geometry import find_mic
OUT=Path(__file__).resolve().parent;REPO=OUT.parents[5];BASE=OUT.parents[3]/'periodic_interface_v4'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
root=BASE/'mace_periodic_v19_residual_response/training_candidate'
m=json.loads((root/'manifest.json').read_text());assert sha(root/'train.extxyz')==m['train_sha256']
frames=read(root/'train.extxyz',':');assert len(frames)==36
ledger=[]
for i,a in enumerate(frames):
 name=a.info.get('config_type','');s=np.array(a.get_chemical_symbols());ag=np.flatnonzero(s=='Ag')
 vec=a.positions[:,None,:]-a.positions[None,:,:];_,dist=find_mic(vec.reshape(-1,3),a.cell,a.pbc);dist=dist.reshape(len(a),len(a));np.fill_diagonal(dist,np.inf)
 contacts={}
 for x in ['Ti','Si','C']:
  ids=np.flatnonzero(s==x);d=dist[np.ix_(ag,ids)]
  contacts[x]={'minimum_A':float(d.min()),'Ag_framework_pairs_within_4p5A':int((d<4.5).sum()),'Ag_coordination_counts_within_4p5A':(d<4.5).sum(axis=1).tolist()}
 interface=next((x for x in ['AgSi','AgTi','AgC'] if name.startswith(x)),'AgTi_rank_recovery')
 tags=[k for k in ['registry','strain','distance','transverse','framework','residual','d2p'] if k in name]
 parent={k:str(a.info[k]) for k in ['source_frame_geometry_sha256','source_geometry_sha256','source_structure_sha256','source_config_type','source_dataset_sha256'] if k in a.info}
 evidence=m['base_rows'][i] if i<30 else m['new_rows'][i-30]
 ledger.append({'frame':i,'config':name,'interface_family_declaration':interface,'formula':a.get_chemical_formula(),'atoms':len(a),'pbc':a.pbc.tolist(),'cell':a.cell.array.tolist(),'declared_perturbation_tags':tags,'declared_mode':a.info.get('proposal_mode'),'thermal_temperature_K':None,'thermal_limitation':'No temperature/ensemble established by these static frames; local_noise is not a thermal trajectory.','contacts_4p5A_geometric_not_bond_assignment':contacts,'parent_metadata':parent,'provenance':evidence,'finite_positions':bool(np.isfinite(a.positions).all())})
proot=BASE/'v19_cif_parent_structures';pm=json.loads((proot/'manifest.json').read_text());props=[]
for r in pm['records']:
 p=REPO/r['path'];assert sha(p)==r['sha256'];a=read(p)
 props.append({'label':r['label'],'interface':r['interface'],'input_sha256':sha(p),'atoms':len(a),'formula':a.get_chemical_formula(),'source_COD':'9009647','source_cif_sha256':pm['source_cif_sha256'],'role':'source-derived label-free proposal; shared parent, not three independent sources'})
paths=subprocess.check_output(['git','ls-files'],cwd=REPO,text=True).splitlines()
assets=[p for p in paths if any(k in p.lower() for k in ['stage69','stage_69','stage-69','bicontinuous'])]
external=[]
for folder in ['/workspace/library-files','/workspace/shared','/workspace/scratch']:
 for p in Path(folder).rglob('*'):
  if p.is_file():external.append(str(p))
counts=dict(collections.Counter(r['interface_family_declaration'] for r in ledger))
result={'status':'COVERAGE_AUDIT_COMPLETE_INDEPENDENT_FAMILY_ASSETS_MISSING','train_count':36,'train_sha256':sha(root/'train.extxyz'),'training_authorized':m['training_authorized'],'roles':m['roles'],'family_counts_by_config_prefix_not_independent_source_counts':counts,'formula_counts':dict(collections.Counter(r['formula'] for r in ledger)),'ledger':ledger,'source_proposals':props,'shared_parent_limit':'All COD slabs derive from COD9009647. Different terminations and new IDs do not supply independent crystallographic sources. Training perturbations inherit few parent families; source identities must be split as groups.','stage69_assets':{'named_tracked_matches':assets,'external_file_names':external,'status':'No identifiable actual Stage69 bicontinuous structure supplied in inspected nonsealed inventory','scope':'Tracked filenames and allowed repository documentation; external library/shared/scratch filenames. Sealed contents and geometries never opened. Arbitrarily named unregistered files are not assumed to be Stage69.'},'missing_assets':['Stage69 atomic positions with species, cell/PBC, units, provenance and generation parameters','Stage69 phase/interface masks and parent/source identifiers','independently sourced alternative parent geometry/provenance for held-out Ag-Ti/Si/C families','role-freeze manifest excluding holdout parents and their perturbations from training and tuning','common numerical target decision and reviewed label/owner registration'],'pure_phase_coverage':{'pure_Ag_rows':0,'Ag_free_Ti3SiC2_rows':0,'limitation':'All36 rows contain Ag,C,Si,Ti; composition rank4 is not pure-phase or collision coverage.'},'recommendations':{'minimal_training_acquisition_proposal':{'count':9,'per_interface':'Three points: normal approach, lateral registry/shear, local framework displacement on a reviewed source-supported TRAIN-only parent. Reuse existing labeled points where geometry and target match; avoid duplicates. Values frozen after contact screen and target decision.','dependency':'reviewed matched numerical method, exact label-free input hashes, parent group roles and owner registrations'},'independent_development_proposal':{'minimum_groups':3,'description':'One distinct nontraining parent group per Ag-Ti/Si/C contact, at least2 geometries/group, development only; freeze before fitting. COD common-parent terminations must be treated as one correlated source group, not three independent sources.'},'independent_test_proposal':{'minimum_groups':3,'description':'Additional untouched parent groups per contact, at least2 geometries/group, completely excluded from training/tuning; labels scored once after selection. Never reuse scored V17 development parents as fresh tests.'},'stability_proposals':['Pure Ag fcc lattice/volume and small displacement response, plus undercoordinated surface/defect check','Ag-free Ti3SiC2 sourced bulk small volume/strain/displacement grid; compare forces and restoring behavior','Finite-temperature validation requires separately authorized trajectories/ensemble evidence; static local noise is insufficient'], 'collision_proposals':['Ag-Ag and Ag-Ti/Si/C short-range pair/embedded-environment scans with fixed separation grid and converged reference','Check repulsive force direction, energy-force agreement, continuity and avoidance of holes; compare physical reference/ZBL where justified without assuming pair-only transfer','No impact or MD launch under this audit']},'numerical_pilots_integration':'Excluded from model training/development; numerical pass does not authorize dataset integration','inputs_opened':[str((root/'train.extxyz').relative_to(REPO))]+[r['path'] for r in pm['records']]}
(OUT/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('PASS training36 hashes/geometry inventory; counts',counts,'Stage69 named matches',assets)
