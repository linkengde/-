"""CIF-derived hexagonal slab proposals; source-supported, surface/DFT review pending."""
import hashlib,json
from pathlib import Path
import numpy as np
from ase.spacegroup import crystal
from ase import Atom
from ase.io import write,read
ROOT=Path(__file__).resolve().parent
CIF=ROOT/'source/Ti3SiC2_COD_9009647.cif'
bulk=read(CIF)
a,c=bulk.cell.lengths()[0],bulk.cell.lengths()[2]
zTi,zC=.13550,.07220
assert bulk.get_chemical_formula()=='C4Si2Ti6'
rows=[];ROOT.joinpath('inputs').mkdir(exist_ok=True)
for k,(contact,gap) in enumerate([('C',2.4),('Si',2.6),('Ti',2.6)]):
 candidates=np.unique(np.round(bulk.get_scaled_positions()[np.array(bulk.get_chemical_symbols())==contact,2],10));top=float(candidates.max());cut=(top+.001)%1
 unit=bulk.copy();f=unit.get_scaled_positions();f[:,2]=(f[:,2]-cut)%1;unit.set_scaled_positions(f)
 slab=unit.repeat((2,2,1));slab.positions[:,2]+=6;topz=slab.positions[:,2].max();assert {slab[i].symbol for i in range(len(slab)) if abs(slab.positions[i,2]-topz)<1e-6}=={contact}
 substrate_count=len(slab);oldcell=slab.cell.copy();newcell=oldcell.copy();newcell[2]=[0,0,topz+gap+8];slab.set_cell(newcell,scale_atoms=False)
 # A coherent hexagonal Ag monolayer with imposed in-plane strain; no relaxation.
 for i in range(2):
  for j in range(2):slab.append(Atom('Ag',position=(i+.37)*unit.cell[0]+(j+.21)*unit.cell[1]+np.array([0,0,topz+gap])))
 slab.pbc=True;slab.wrap();ids=np.arange(19100000+k*1000,19100000+k*1000+len(slab));slab.new_array('lammps_id',ids);markers=np.zeros(len(slab),dtype=int)
 d=slab.get_all_distances(mic=True);np.fill_diagonal(d,np.inf);ag=np.arange(substrate_count,len(slab));xs=np.array([i for i in range(substrate_count) if slab[i].symbol==contact]);ii,jj=np.unravel_index(np.argmin(d[np.ix_(ag,xs)]),(len(ag),len(xs)));markers[ag[ii]]=markers[xs[jj]]=1;slab.new_array('central_pair',markers)
 label=f'Ag{contact}_cod9009647_slab_v19_dev_proposal';slab.info={'config_type':label,'dataset_role':'CIF_DERIVED_DEVELOPMENT_PROPOSAL','source_method':'COD9009647 uploaded CIF; Kisi et al1998 DOI10.1016/S0022-3697(98)00226-1','new_ids_not_historical':True}
 path=ROOT/'inputs'/f'{label}.extxyz';write(path,slab,format='extxyz');stored=read(path);minima={}
 for i in range(len(slab)):
  for j in range(i):
   pair='-'.join(sorted([slab[i].symbol,slab[j].symbol]));minima[pair]=min(minima.get(pair,float('inf')),float(d[i,j]))
 rows.append({'interface':'Ag'+contact,'label':label,'path':str(path.relative_to(ROOT.parents[3])),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'atoms':len(slab),'formula':slab.get_chemical_formula(),'top_surface_species':contact,'gap_plane_A':gap,'species_pair_minima_A':minima,'new_id_range':[int(ids.min()),int(ids.max())],'lineage':'Constructed from uploaded COD9009647 CIF; no legacy cluster coordinates or labels used','status':'PROPOSAL_ONLY_PENDING_CRYSTAL_AND_GEOMETRY_REVIEW','marked_pair_ids':[int(ids[ag[ii]]),int(ids[xs[jj]])],'periodic_topology':'2x2 hexagonal in-plane substrate; one conventional cell thickness; four Ag monolayer atoms'} )
manifest={'status':'CIF_SUPPORTED_PROPOSALS_PENDING_SURFACE_AND_NUMERICAL_REVIEW','source_parameters':{'spacegroup':194,'a_A':a,'c_A':c,'Ti_4f_z':zTi,'C_4f_z':zC,'basis':[{'species':'Ti','fractional':[0,0,0]},{'species':'Ti','fractional':[0.66667,0.33333,zTi]},{'species':'Si','fractional':[0,0,.25]},{'species':'C','fractional':[0.33333,0.66667,zC]}]},'source_cif_sha256':hashlib.sha256(CIF.read_bytes()).hexdigest(),'COD_record':'9009647','DOI':'10.1016/S0022-3697(98)00226-1','source_URL':'https://www.crystallography.net/cod/9009647.cif','journal':'Journal of Physics and Chemistry of Solids59(1998)1437-1443','occupancy_note':'CIF omits explicit occupancy column; ASE uses full occupancy default; no defect occupation claimed','records':rows,'limitations':['User-uploaded CIF internally verified; remote COD checksum not independently retrieved.','Unrelaxed cleaved surfaces and commensurate Ag monolayer are stress probes, not equilibrium structures.','Thin slab, vacuum, k-point/cutoff convergence and Ag coherent strain require review.','New ID/atom count alone does not prove independence; geometric/topological audit still required.','No production candidate mass fraction or morphology claim.','No sealed inputs/labels accessed.']};(ROOT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print([(r['interface'],r['atoms'],min(r['species_pair_minima_A'].values())) for r in rows])
