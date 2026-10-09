"""Reuse verified input bytes and derive reviewed 3D mesh proposals; no DFT."""
from pathlib import Path
import collections,hashlib,json,shutil
import numpy as np
from ase.io import read
from entry import ROOT,REPO,ENERGY,sha
DESIGN=ROOT.parent/'v19_pure_phase_and_short_range_entry_design'
m=json.loads((DESIGN/'input_manifest.json').read_text())
setups=Path('/workspace/.venvs/gpaw-mpi/lib/python3.12/site-packages/gpaw_data/setups')
method={'xc':'PBE','basis':'plane wave','cutoff_eV':500,'smearing_eV':.1,'parallel':{'sl_auto':True},'scf_thresholds':'energy/density/eigenstates 1e-5','mixer':'Pulay beta=0.05, nmaxold=8, weight=100','maxiter':160,'GPAW_version':'26.7.0','ASE_version':'3.29.0','PAW_data_version':'1.2.1','PAW_setup_sha256':{n:sha(setups/n) for n in ['Ag.PBE.gz','Ti.PBE.gz','Si.PBE.gz','C.PBE.gz']},'energy_keys':{'native':'energy (force_consistent=False)','force_consistent':'free_energy'}}
plans={};records=[]
for phase in ['Ag','Ti3SiC2']:
 baseline=next(r for r in m['records'] if r['phase']==phase and r['operation']['kind']=='baseline')
 a=read(DESIGN/baseline['input']);B=2*np.pi*np.linalg.inv(a.cell.array).T;norm=np.linalg.norm(B,axis=1)
 steps=[]
 for spacing in [.40,.25,.16,.12]:
  raw=np.ceil(norm/spacing).astype(int);mesh=(2*np.ceil(raw/2)).astype(int);actual=norm/mesh
  assert np.all(actual<=spacing+1e-12)
  if any(s['kpts']==mesh.tolist() for s in steps):continue
  steps.append({'target_max_spacing_invA':spacing,'kpts':mesh.tolist(),'actual_axis_spacing_invA':actual.tolist(),'raw_kpoints_upper_bound':int(np.prod(mesh)),'relative_raw_work':float(np.prod(mesh))})
 for s in steps:s['relative_raw_work']/=steps[0]['raw_kpoints_upper_bound']
 plans[phase]={'reciprocal_vectors_2pi_invA':B.tolist(),'reciprocal_vector_lengths_invA':norm.tolist(),'supercell_atoms':len(a),'steps':steps,'selected_proposed_control_mesh':steps[1]['kpts'],'cost_limit':'Raw kpoint count ratios only, not measured runtime; symmetry, bands, plane waves, SCF and checkpoint costs unknown.'}
for r in m['records']:
 source=DESIGN/r['input'];destination=ROOT/r['input'];assert sha(source)==r['sha256'];assert sha(destination)==r['sha256']
 a=read(source);mesh=plans[r['phase']]['selected_proposed_control_mesh']
 variants=[mesh] if r['operation']['kind']!='baseline' else [s['kpts'] for s in plans[r['phase']]['steps']]
 for kpts in variants:
  label=r['label'].replace('_PROPOSAL','')+'_k'+'x'.join(map(str,kpts))+'_DRAFT'
  records.append({'label':label,'phase':r['phase'],'input':r['input'],'source_input':str(source.relative_to(REPO)),'input_sha256':r['sha256'],'ordered_symbols':a.get_chemical_symbols(),'ordered_ids':a.arrays['local_new_id'].tolist(),'composition':dict(collections.Counter(a.get_chemical_symbols())),'kpts':kpts,'pbc':[True,True,True],'poissonsolver':{},'method':dict(method,kpts=kpts,pbc=[True,True,True],poissonsolver={}),'owner_task':'window-b','owner_instance':None,'launch_enabled':False,'role':'numerical_pure_phase_control_not_training_or_independent_validation','status':'blocked_pending_A_registration','operation':r['operation']})
manifest={'launch_enabled':False,'method':method,'energy_convention':ENERGY,'records':records,'source_design_manifest_sha256':sha(DESIGN/'input_manifest.json'),'source_CIF_sha256':m['source_cif_sha256'],'disk_budget':{'minimum_start_bytes':3221225472,'minimum_before_checkpoint_bytes':1610612736}}
(ROOT/'input_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
(ROOT/'mesh_plan.json').write_text(json.dumps({'phases':plans,'convention':'reciprocal vectors include2pi; mesh ceil(norm/target) rounded up to even, fixed Monkhorst-Pack parity','comparisons':['native and free energy deltas/meV per atom','all-atom vector RMS and max force differences/eV per Angstrom','Optional stress only if explicitly implemented: ASE eV/A3, tension-positive, Voigt xx yy zz yz xz xy; no stress computed/exported by current draft'],'acceptance':'Budgets proposed2meV/atom,0.01eV/A vector; pair gate inapplicable. Consecutive meshes plus displaced/strained controls and width/cutoff checks required before convergence claim.','equilibrium_caveat':'Symmetry can make bulk forces zero at every mesh; baseline force pass is not displaced/strain force accuracy.','roles':'Pure-phase perturbations share parents; numerical controls only, not independent interface validation.'},indent=2)+'\n')
print('Prepared',len(records),'disabled records; plans', {p:[s['kpts'] for s in v['steps']] for p,v in plans.items()})
