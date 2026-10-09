"""Deterministic label-free pure-phase controls; no calculator or collision inputs."""
import hashlib,json,inspect
from pathlib import Path
import numpy as np
import ase
from ase.build import bulk
from ase.io import read,write
from ase.data import atomic_numbers,reference_states,covalent_radii
from ase.neighborlist import neighbor_list
OUT=Path(__file__).resolve().parent;REPO=OUT.parents[5]
SOURCE=REPO/'research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures/source/Ti3SiC2_COD_9009647.cif'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(SOURCE)=='6abbf0a40208c5985cacc3b9ac3564a944702203ad57e961010a674f30d598c3'
ag=bulk('Ag','fcc',a=float(reference_states[atomic_numbers['Ag']]['a']),cubic=True).repeat((2,2,2))
ti=read(SOURCE).repeat((2,2,1));assert ti.get_chemical_formula()=='C16Si8Ti24'
inputs=OUT/'inputs';inputs.mkdir(exist_ok=True)
records=[]
for phase,parent in [('Ag',ag),('Ti3SiC2',ti)]:
 for case in ['baseline','volume_m02','volume_p02','displacement_m01','displacement_p01','shear_m005','shear_p005']:
  a=parent.copy();operation={'kind':case}
  if case.startswith('volume'):
   factor=.98 if 'm02' in case else 1.02;a.set_cell(a.cell.array*factor**(1/3),scale_atoms=True);operation['volume_ratio']=factor
  elif case.startswith('displacement'):
   amount=-.01 if 'm01' in case else .01;a.positions[0,0]+=amount;operation.update(atom_index=0,direction=[1,0,0],displacement_A=amount)
  elif case.startswith('shear'):
   gamma=-.005 if 'm005' in case else .005;F=np.eye(3);F[0,1]=gamma;a.set_cell(a.cell.array@F.T,scale_atoms=True);operation.update(deformation_gradient=F.tolist(),engineering_shear=gamma)
  a.pbc=True;a.arrays['local_new_id']=np.arange(1,len(a)+1,dtype=int);a.info={'proposal_role':'pure_phase_numerical_stability_control_NOT_VALIDATED','phase':phase,'case':case,'identity_note':'Local new IDs only; not historical lammps identities'}
  assert np.isfinite(a.positions).all() and np.linalg.det(a.cell)>0
  i,j,d=neighbor_list('ijd',a,cutoff=8.,self_interaction=False)
  symbols=a.get_chemical_symbols();pairs={};screen=True
  for ii,jj,dd in zip(i,j,d):
   pair='-'.join(sorted([symbols[ii],symbols[jj]]));pairs[pair]=min(pairs.get(pair,float('inf')),float(dd))
   threshold=.55*(covalent_radii[a.numbers[ii]]+covalent_radii[a.numbers[jj]])
   screen=screen and dd>threshold
  assert screen
  label=f'{phase}_{case}_PROPOSAL';p=inputs/(label+'.extxyz');write(p,a,format='extxyz');b=read(p)
  assert np.max(np.abs(b.positions-a.positions))<1e-8 and np.max(np.abs(b.cell.array-a.cell.array))<1e-8
  assert np.array_equal(b.numbers,a.numbers) and np.array_equal(b.arrays['local_new_id'],a.arrays['local_new_id'])
  assert not any(k in b.info for k in ['energy','REF_energy','free_energy']) and b.calc is None
  records.append({'label':label,'input':str(p.relative_to(OUT)),'sha256':sha(p),'phase':phase,'formula':a.get_chemical_formula(),'atoms':len(a),'operation':operation,'volume_A3':a.get_volume(),'species_pair_periodic_minima_A':pairs,'screen':'PASS_GEOMETRY_TRIAGE_ONLY','role':'proposed_pure_phase_control','launch_enabled':False,'owner_instance':None,'numerical_target':{'xc':'PBE_PROPOSED','cutoff_eV':None,'kpts':None,'smearing_eV':None,'energy_key':None},'blocked_reason':'A fixed target decision, phase-specific numerical convergence, exact job/method/owner registration required'})
manifest={'launch_enabled':False,'source_cif_sha256':sha(SOURCE),'source_cif_path':str(SOURCE.relative_to(REPO)),'Ag_source':{'source':'ASE3.29 reference_states Ag fcc a4.09A','a_A':4.09,'status':'explicit modeling starting value, neither independently verified experiment nor optimized lattice'},'Ti_source':{'COD':'9009647','status':'uploaded CIF with internal provenance; remote checksum not independently verified','stoichiometry':'Ti3SiC2','supercell':[2,2,1]},'Ag_supercell':[2,2,2],'local_id_policy':'1..N within each phase input; inherited phase IDs only for paired controls in this pack, no historical IDs invented','screen_policy':{'periodic_neighbors':'ASE neighbor_list includes periodic image neighbors','cutoff_A':8,'pair_floor_A':'0.55*(ASE covalent radius_i+radius_j)','meaning':'conservative overlap triage, no equilibrium/stability assertion'},'roundtrip_tolerance_A':1e-8,'records':records,'tool_versions':{'ASE':ase.__version__,'numpy':np.__version__},'source_tool_sha256':{'ase_reference_data':sha(Path(inspect.getfile(__import__('ase.data',fromlist=['x'])))),'ase_bulk_builder':sha(Path(inspect.getfile(bulk)))}}
(OUT/'input_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
grid=[.8,1.,1.25,1.5,1.75,2.,2.5,3.,4.,5.]
collision={'launch_enabled':False,'owner_instance':None,'pairs':['Ag-Ag','Ag-Ti','Ag-Si','Ag-C'],'distance_grid_A':grid,'finite_difference_step_A':.005,'no_collision_inputs_generated':True,'isolated_pair_boundary_proposal':'nonperiodic neutral two-atom cell with vacuum>=10A each side, spin/PAW/cutoff convergence reviewed; do not transfer metallic slab mesh blindly','embedded_environment_proposal':'TRAIN-only reviewed pure/surface parent, move one local Ag relative to target Ti/Si/C or Ag; fixed surrounding ions; all other species-pair distances screened; exclude development/test groups','metrics':{'repulsive_force':'rhat points i->j, separating_force=dot(Fj-Fi,rhat)/2; positive for repulsion','energy_force':'separating_force compared with -(E(r+h)-E(r-h))/(2h) for symmetric pair displacement; compare native/free conventions separately and use force-consistent target for finite-width derivative','continuity':'finite E/F and smooth adjacent values/slopes; central differences plus h/2 spot checks, no arbitrary smoothness pass threshold before physical review','far_separation':'consistent asymptote after interaction isolation/boundary convergence'},'references_and_limits':['MACE multi-element model and current mixed potential must each declare species coverage, hybrid pair mapping, cutoff/units and any ZBL splice before comparable curves; current mixed-potential files are not established by this pack','ZBL assumes screened nuclear Coulomb repulsion with atomic numbers and screening; is a high-energy short-range asymptotic comparison, not an equilibrium chemical label','Sub-Angstrom PAW calculations may violate frozen-core/setup validity. Short grid points are design candidates only, not approved DFT or production inputs; hard-core subset requires reference/setup review','Pair scans do not demonstrate many-body or impact stability; embedded environments and pure-phase controls address different gaps'],'grid_point_status':[{'r_A':r,'status':'REVIEW_REQUIRED_NO_INPUT_NO_LAUNCH'} for r in grid]}
(OUT/'collision_scan_plan.json').write_text(json.dumps(collision,indent=2)+'\n')
report={'status':'PREPARATION_VERIFIED_LAUNCH_DISABLED','input_count':len(records),'phase_control_counts':{'Ag':7,'Ti3SiC2':7},'geometry_screens_pass':True,'labels_present':False,'DFT_executed':False,'manifest':manifest,'collision_plan':collision,'pure_phase_verifier_requirements':['exact ordered species/positions/cell/PBC/local IDs/hash and expected phase composition; no all-four-species or marked-interface-pair requirement','finite energies and Nx3 forces; explicit SCF-converged flag and actual iterations/MPI ranks','exact PAW/code/cutoff/mesh/width/energy/Poisson/magnetism convention and source hashes','store native energy/free energy separately; preserve potential energy per atom, force vector conventions and isolated atomic reference choices','symmetric volume/displacement/shear identification; test restoring response and derivatives with target-specific numerical errors, not arbitrary apparent energy decrease','checksummed immutable complete archive, owner/no-existing-run/disk guards and genuine convergence results'], 'missing_evidence':['independently verified Ag experimental source or optimized lattice; ASE4.09 is only starting value','phase-specific convergence and consistent numerical target after A6/8 decision','current mixed-potential source files/configuration and species/units mapping','short-range PAW/ZBL validity and splice/reference policy','Stage69 morphology/atomic inputs remain missing'],'scientific_limits':['Static controls do not establish physical stability, thermostat behavior or finite-temperature lifetime','No collision/overlapping production inputs or MD authorization','Pure phases are proposed controls, not automatically independent interface-validation families']}
(OUT/'report.json').write_text(json.dumps(report,indent=2)+'\n')
print('PASS: fourteen label-free controls, periodic screens and roundtrip; collision plan only, no launch')
