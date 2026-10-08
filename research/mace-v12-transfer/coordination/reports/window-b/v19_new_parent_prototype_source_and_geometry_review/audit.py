import json,hashlib,subprocess
from pathlib import Path
import numpy as np
from ase.io import read
from ase.geometry import find_mic
ROOT=Path(__file__).resolve().parents[6];OUT=Path(__file__).resolve().parent;BASE=ROOT/'research/mace-v12-transfer/periodic_interface_v4';SRC=BASE/'v19_new_parent_prototypes';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
m=json.loads((SRC/'manifest.json').read_text());assert m['status']=='UNAPPROVED_PROTOTYPES_DFT_NOT_AUTHORIZED';pins={}
for line in (SRC/'SHA256SUMS.txt').read_text().splitlines():h,n=line.split(None,1);assert sha(SRC/n)==h;pins[str((SRC/n).relative_to(ROOT))]=h
# Source lookup is deliberately recorded without claiming success; no external request produced usable evidence.
reference_search={'attempted':'Crossref REST query for Ti3SiC2 crystal structure / lattice constants','endpoint':'https://api.crossref.org/works?query.title=Ti3SiC2%20crystal%20structure%20lattice%20parameters&rows=5','result':'UNAVAILABLE: network tunnel returned HTTP proxy 403; no source document retrieved','repository_search':'No committed citation/table for these values found in research/mace-v12-transfer','consequence':'P63/mmc, a=3.07,c=17.68,z(Ti4f)=0.135,z(C4f)=0.567 remain unverified assumptions; DFT readiness cannot pass crystallographic-source gate'}
# training36 geometry only
trainpath=BASE/'mace_periodic_v19_residual_response/training_candidate/train.extxyz';assert sha(trainpath)==json.loads((BASE/'mace_periodic_v19_residual_response/training_candidate/manifest.json').read_text())['train_sha256']
from importlib.util import spec_from_file_location,module_from_spec
spec=spec_from_file_location('geom',ROOT/'research/mace-v12-transfer/coordination/reports/window-b/v17_split_geometry_design/final_set/generate_split_candidates.py');g=module_from_spec(spec);spec.loader.exec_module(g)
train=g.extxyz_geometry_frames(trainpath);results=[]
for r in m['records']:
 p=ROOT/r['path'];assert sha(p)==r['sha256'];a=g.extxyz_geometry_frames(p)[0];assert len(a)==52 and a.get_chemical_formula()==r['formula'] and np.all(a.pbc)
 z=a.positions[:,2];unique=sorted(set(np.round(z,8)));layer_report={str(sym):sorted(set(np.round(z[np.array(a.get_chemical_symbols())==sym],5).tolist())) for sym in sorted(set(a.get_chemical_symbols()))}
 pairs=a.get_all_distances(mic=True);sy=a.get_chemical_symbols();mins={}
 for i in range(len(a)):
  for j in range(i):
   key='-'.join(sorted((sy[i],sy[j])));mins[key]=min(mins.get(key,float('inf')),float(pairs[i,j]))
 # in-plane cross-boundary minimum and Ag coherence
 vec=a.cell.array[:2];agb=np.array([a[i].position for i in range(len(a)) if sy[i]=='Ag']);lateral=[]
 for i in range(len(agb)):
  for j in range(i):
   d,_=find_mic(agb[i]-agb[j],a.cell,pbc=[True,True,False]);lateral.append(float(np.linalg.norm(d)))
 top=max(z);bottom=min(z);substrate_top=float(z[:48].max()); substrate_bottom=float(z[:48].min())
 # Manifest gap is the vertical spacing from the selected substrate surface to its Ag plane.
 agz=float(z[np.array(sy)=='Ag'][0]);vacuum=float(a.cell[2,2]-(top-bottom));assert vacuum>7.9
 substrate_top_species=sorted(set(sy[i] for i in range(48) if abs(z[i]-substrate_top)<1e-7));substrate_bottom_species=sorted(set(sy[i] for i in range(48) if abs(z[i]-substrate_bottom)<1e-7));actual_top_species=sorted(set(sy[i] for i in range(len(a)) if abs(z[i]-top)<1e-7))
 ids=a.arrays['lammps_id'];assert len(set(ids))==52 and len(a.arrays['central_pair'])==52 and int(a.arrays['central_pair'].sum())==2
 assert not any(k.lower().startswith(('ref_','pw_pbe','force','energy')) for k in list(a.info)+list(a.arrays))
 near=[]
 for j,b in enumerate(train):
  if a.get_chemical_formula()==b.get_chemical_formula():
   c=g.compare_geometry(a,b)
   if c and (c['exact_duplicate'] or c['near_duplicate']):near.append({'frame':j,'match':c})
 results.append({'interface':r['interface'],'input':r['path'],'sha256':sha(p),'atoms':len(a),'formula':a.get_chemical_formula(),'cell_A':a.cell.array.tolist(),'pbc':a.pbc.tolist(),'z_extent_A':[bottom,top],'substrate_z_extent_A':[substrate_bottom,substrate_top],'substrate_bottom_surface_species':substrate_bottom_species,'substrate_upper_surface_species':substrate_top_species,'Ag_plane_z_A':agz,'actual_total_top_surface_species':actual_top_species,'vertical_substrate_Ag_gap_A':agz-substrate_top,'vacuum_in_c_A':vacuum,'unique_layer_z_A':unique,'element_layer_z_A':layer_report,'pair_minima_recomputed_A':mins,'manifest_pair_minima_agreement_max_A':max(abs(mins[k]-r['species_pair_minima_A'][k]) for k in mins),'Ag_in_plane_pair_min_A':min(lateral),'Ag_in_plane_pair_max_A':max(lateral),'ideal_substrate_a_A':float(np.linalg.norm(vec[0])),'Ag_coherent_inplane_strain_percent_vs_assumed_a':100*(np.linalg.norm(vec[0])-3.07)/3.07,'marked_ids':ids[np.flatnonzero(a.arrays['central_pair'])].tolist(),'labels_present':False,'exact_or_near_same_formula_train36_hits':near,'nonAg_surface_at_zmax':sorted(set(sy[i] for i in range(len(a)) if abs(z[i]-top)<1e-7)),'geometry_checks':'PASS; physical equilibrium not assessed'})
assert len(results)==3
report={'status':'GEOMETRY_INDEPENDENCE_SCREEN_PASS_SOURCE_AND_PHYSICAL_APPROVAL_BLOCKED','manifest_hash':sha(SRC/'manifest.json'),'builder_hash':sha(SRC/'build_prototypes.py'),'source_search':reference_search,'prototypes':results,'training36_formula_set':sorted(set(x.get_chemical_formula() for x in train)),'training36_frame_count':len(train),'same_composition_exact_or_near_match_count':sum(len(x['exact_or_near_same_formula_train36_hits']) for x in results),'sealed_inputs_labels_opened':False,'interpretation':'Different composition and periodic layered topology from mixed small-cluster training frames supports geometry novelty within the repository, not independent validation or physical reliability. All share the same assumed ideal Ti3SiC2 source motif; three terminations are not three independent bulk source families.'}
(OUT/'audit.json').write_text(json.dumps(report,indent=2)+'\n');print([(r['interface'],r['vacuum_in_c_A'],r['pair_minima_recomputed_A']['C-Ti'],r['Ag_in_plane_pair_min_A'],r['Ag_in_plane_pair_max_A']) for r in results])
