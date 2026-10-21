from __future__ import annotations
import hashlib, json
from collections import Counter
from pathlib import Path
import numpy as np
import ase
from ase.build import surface
from ase.dft.kpoints import monkhorst_pack
from ase.io import read, write

ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[5]
CIF=REPO/'research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures/source/Ti3SiC2_COD_9009647.cif'
CIF_SHA='6abbf0a40208c5985cacc3b9ac3564a944702203ad57e961010a674f30d598c3'
VAC=18.0
NLAYERS=10
DISP=0.02

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def geom_sha(a):
    h=hashlib.sha256(); h.update(' '.join(a.get_chemical_symbols()).encode())
    h.update(np.asarray(a.cell,dtype='<f8').tobytes());h.update(np.asarray(np.round(a.positions,8),dtype='<f8').tobytes())
    h.update(np.asarray(a.pbc,dtype=np.uint8).tobytes());h.update(np.asarray(a.arrays['source_atom_id'],dtype='<i8').tobytes());return h.hexdigest()
def inversion_map(a):
    z=a.positions[:,2]; lo=float(z.min());hi=float(z.max());frac=a.get_scaled_positions(wrap=True);sym=a.get_chemical_symbols();pairs=[];errs=[]
    for i,s in enumerate(sym):
        targetxy=(-frac[i,:2])%1;targetz=lo+hi-z[i];cand=[]
        for j,t in enumerate(sym):
            if s!=t:continue
            dxy=(frac[j,:2]-targetxy+.5)%1-.5;drxy=dxy@a.cell[:2];dz=a.positions[j,2]-targetz
            cand.append((float(np.sqrt(drxy@drxy+dz*dz)),j))
        e,j=min(cand);errs.append(e)
        if i<j:pairs.append((i+1,j+1,e))
    return max(errs),pairs
def min_pair(a):
    d=a.get_all_distances(mic=True);np.fill_diagonal(d,np.inf);i,j=np.unravel_index(np.argmin(d),d.shape);return float(d[i,j]),int(i+1),int(j+1)
def plane(a,z):
    ix=np.flatnonzero(np.abs(a.positions[:,2]-z)<1e-6)
    return {'z_A':float(z),'atom_ids':[int(i+1) for i in ix],'composition':dict(sorted(Counter(a[i].symbol for i in ix).items()))}

if ase.__version__!='3.29.0': raise RuntimeError(f'ASE version mismatch: {ase.__version__}')
if sha(CIF)!=CIF_SHA: raise RuntimeError('source CIF SHA mismatch')
bulk=read(CIF,format='cif')
if len(bulk)!=12 or bulk.get_chemical_formula()!='C4Si2Ti6': raise RuntimeError('unexpected source cell')
# (1,1,0) on the conventional hexagonal cell corresponds to the (11-20) prism.
shifted=bulk.copy();shifted.positions+=0.5*shifted.cell[0]
slab=surface(shifted,(1,1,0),NLAYERS,vacuum=None,periodic=False)
slab.center(vacuum=VAC/2,axis=2);slab.set_pbc((True,True,True));slab.new_array('source_atom_id',np.arange(1,len(slab)+1,dtype=np.int64))
counts=dict(sorted(Counter(slab.get_chemical_symbols()).items()))
if len(slab)!=120 or counts!={'C':40,'Si':20,'Ti':60}: raise RuntimeError(f'count/composition mismatch: {len(slab)} {counts}')
inv_err,pairs=inversion_map(slab)
if inv_err>1e-5: raise RuntimeError(f'inversion symmetry failed {inv_err}')
md,mi,mj=min_pair(slab)
if md<1.5: raise RuntimeError(f'unphysical minimum distance {md} at {mi}/{mj}')
zlo=float(slab.positions[:,2].min());zhi=float(slab.positions[:,2].max())
terminal_ti=[i for i,a in enumerate(slab) if a.symbol=='Ti' and (abs(a.position[2]-zlo)<1e-6 or abs(a.position[2]-zhi)<1e-6)]
terminal_pairs=[]
for i,j,e in pairs:
    if i-1 in terminal_ti and j-1 in terminal_ti and slab[i-1].symbol=='Ti' and slab[j-1].symbol=='Ti':terminal_pairs.append((i,j,e))
if not terminal_pairs: raise RuntimeError('no symmetry-paired terminal Ti atoms')
lo_id,hi_id,_=min(terminal_pairs,key=lambda x:(x[0],x[1]))
ref=slab.copy();ref.info.update({'structure_id':'TSC120_DEV_PRISM_11m20_COD9009647_10L_REF','source_family':'surface_prismatic_11m20_cod9009647_10layer','candidate_role':'DEV_GEOMETRY_FAMILY_REFERENCE_CANDIDATE','material':'Ti3SiC2','surface_orientation_ASE_Miller':[1,1,0],'surface_orientation_hexagonal_4index_equivalent':'(11-20)','surface_termination':'symmetric mixed stoichiometric prismatic termination','source_cif':'research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures/source/Ti3SiC2_COD_9009647.cif','source_cif_sha256':CIF_SHA,'construction':'ASE surface; conventional cell translated by 0.5*a before (1,1,0) cleavage; 10 layers; no relaxation/strain','shared_crystallographic_root_with_TRAIN':True,'independence_scope':'distinct facet geometry parent; not independent experimental-source blind test','fixed_geometry':True,'vacuum_gap_angstrom':VAC})
resp=ref.copy();resp.info=dict(ref.info);resp.info.update({'structure_id':'TSC120_DEV_PRISM_11m20_COD9009647_10L_RESP_TIZ020','candidate_role':'DEV_GEOMETRY_FAMILY_RESPONSE_CANDIDATE','displacement_definition':'inversion-paired terminal Ti atoms moved oppositely along surface-normal z','displacement_magnitude_angstrom':DISP,'displaced_source_atom_ids':[lo_id,hi_id]})
resp.positions[lo_id-1,2]-=DISP;resp.positions[hi_id-1,2]+=DISP
for name,a in [('REF',ref),('RESP_TIZ020',resp)]:write(ROOT/'inputs'/f'{a.info["structure_id"]}.extxyz',a,format='extxyz')
# K mesh preserves or improves the old 12x2x1 reciprocal spacing of the (10-10) candidate.
kpts=(7,2,1);kp=monkhorst_pack(kpts);kp_sha=hashlib.sha256(np.asarray(kp,dtype='<f8').tobytes()).hexdigest()
entries=[]
for a in (ref,resp):
    entries.append({'structure_id':a.info['structure_id'],'input_path':f'inputs/{a.info["structure_id"]}.extxyz','input_sha256':sha(ROOT/'inputs'/f'{a.info["structure_id"]}.extxyz'),'geometry_sha256':geom_sha(a),'atom_count':len(a),'composition':counts,'source_family':a.info['source_family'],'role_candidate':a.info['candidate_role'],'displaced_source_atom_ids':a.info.get('displaced_source_atom_ids',[])})
(ROOT/'TASKS.json').write_text(json.dumps({'schema_version':'1.0','source_cif': 'research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures/source/Ti3SiC2_COD_9009647.cif','source_cif_sha256':CIF_SHA,'source_git_commit_at_preparation':'8886461261f8f81f775efbab89fe57a22bbd6a91','role_policy':'geometry-family candidates only; existing TRAIN/DEV/TEST assignments are unchanged','tasks':entries},indent=2)+'\n')
config={'schema_version':'1.0','scope':'Two fixed-geometry PBE/PAW single points: 120-atom (11-20) prismatic reference and symmetric terminal-Ti displacement response','role':'DEV_GEOMETRY_FAMILY_CANDIDATE_ONLY; same COD source root, not independent crystallographic source','source_family':entries[0]['source_family'],'method':{'xc':'PBE','mode':'plane-wave','cutoff_eV':500,'spin_polarized':False,'kpoint_type':'standard Monkhorst-Pack using GPAW kpts tuple; no explicit gamma flag or additional shift','kpoints_size':list(kpts),'kpoint_fractional_coordinates_sha256':kp_sha,'kpoint_mesh_rationale':'7x2x1 keeps reciprocal spacing at least as dense as the prior 12x2x1 mesh for in-plane lengths 3.0575x17.6235 A; new in-plane lengths are recorded in STRUCTURE_AUDIT.json. No k-point convergence study is implied.','fermi_dirac_width_eV':0.10,'scf_convergence':{'energy':1e-5,'density':1e-5,'eigenstates':1e-5},'maxiter':160,'eigensolver':'GPAW Davidson niter=2','mixer':{'type':'MSR1','nmaxold':10,'beta':0.05,'reg':0.005,'gb_scale':1.0,'max_A':0.75,'trust_scale':1.0,'soft_lim':4.0,'hard_lim':6.0,'weight':20.0,'sigma':0.015},'sl_auto':True,'fixed_coordinates':True,'fixed_cell':True,'relaxation':False,'pbc':[True,True,True],'vacuum_angstrom':VAC},'software':{'python':'3.12.14','GPAW_version':'26.7.0','ASE_version':'3.29.0','gpaw_data_version':'1.2.1','MPI_backend':'cgpaw'},'paw_sha256':{'Ti.PBE.gz':'bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f','Si.PBE.gz':'2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0','C.PBE.gz':'fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983'},'parallel':{'MPI_ranks':4,'OMP_NUM_THREADS':1,'GPAW_MPI_BACKEND':'cgpaw'},'resources':{'CPU_quota_cores':4,'memory_limit_bytes':34359738368,'memory_soft_guard_bytes':27917287424,'workspace_available_bytes_at_build':None,'walltime_seconds_per_point':18000,'checkpoint_policy':'One rotating full-wavefunction single-precision checkpoint every 20 SCF iterations; preserve prior good checkpoint until the new one is hashed and atomically installed; resume only after SHA and geometry validation.','estimated_checkpoint_bytes':3200000000,'resource_estimate_basis':'120 atoms and 14 k-points compared with archived 96-atom 24-kpoint checkpoint 2,422,680,386 bytes; quadratic atom/k-point scaling estimates 2.65 GB, rounded up to 3.2 GB. Estimate only; actual checkpoint and memory measured at runtime.'},'authorization':{'authorized_by_current_user':True,'scope':'Build this (11-20) facet geometry candidate and run only the listed two fixed-geometry single points sequentially if preflight passes; no role changes, TEST access, relaxation, parameter changes, or unrelated jobs.','c_approval_required':False}}
(ROOT/'SCIENTIFIC_CONFIG.json').write_text(json.dumps(config,indent=2)+'\n')
cell=slab.cell.array.tolist();plane_lo=plane(slab,zlo);plane_hi=plane(slab,zhi)
audit={'status':'GEOMETRY_CANDIDATE_PASS_NO_ROLE_CHANGE','ase_version':ase.__version__,'source_cif_sha256':CIF_SHA,'miller_ase_3index':[1,1,0],'miller_bravais_equivalent':'(11-20)','layers':NLAYERS,'atom_count':len(slab),'composition':counts,'cell_angstrom':cell,'pbc':slab.pbc.tolist(),'surface_span_angstrom':zhi-zlo,'vacuum_gap_angstrom':float(slab.cell[2,2]-(zhi-zlo)),'surface_termination_planes':[plane_lo,plane_hi],'inversion_symmetry_max_error_angstrom':inv_err,'minimum_distance_angstrom':md,'minimum_distance_atom_ids':[mi,mj],'terminal_ti_pair_source_ids':[lo_id,hi_id],'response_displacement_angstrom':DISP,'response_max_atom_displacement_angstrom':max(np.linalg.norm(resp.positions-ref.positions,axis=1)),'response_minimum_distance_angstrom':min_pair(resp)[0],'in_plane_cell_lengths_angstrom':slab.cell.lengths()[:2].tolist(),'kpoint_mesh':[7,2,1],'kpoint_mesh_sha256':kp_sha,'warning':'Distinct (11-20) facet geometry parent from the existing (10-10) slab, but both share COD 9009647 crystallographic root; not an independent source blind test. 10-layer thickness is a candidate and has not been thickness-converged.'}
(ROOT/'STRUCTURE_AUDIT.json').write_text(json.dumps(audit,indent=2,sort_keys=True)+'\n')
(ROOT/'KPOINT_AUDIT.json').write_text(json.dumps({'mesh':list(kpts),'convention':'ASE/GPAW standard Monkhorst-Pack; kpts tuple only; no gamma=True or extra shift','fractional_coordinates':kp.tolist(),'fractional_coordinate_sha256':kp_sha,'old_reference_mesh':[12,2,1],'old_in_plane_lengths_angstrom':[3.0575,17.6235],'new_in_plane_lengths_angstrom':slab.cell.lengths()[:2].tolist(),'new_mesh_spacing_not_coarser_by_axis':[float(slab.cell.lengths()[0]*kpts[0])>=3.0575*12,float(slab.cell.lengths()[1]*kpts[1])>=17.6235*2]},indent=2)+'\n')
(ROOT/'README.md').write_text('''# Ti3SiC2 (11-20) prismatic DEV geometry-family candidate\n\nA symmetric 120-atom, 10-layer slab from COD 9009647, with 18 A vacuum and an inversion-paired terminal-Ti response geometry. This is a new facet geometry family compared with the existing (10-10) prismatic candidate, but shares the COD crystallographic source. It is not an independent crystallographic-source blind test.\n\nThe geometry is unrelaxed and unstrained. Candidate role only; no TRAIN/DEV/TEST assignments are changed. The plan is two sequential PBE/PAW fixed-geometry single points, 500 eV, 7x2x1 standard Monkhorst-Pack, Fermi-Dirac 0.10 eV. The 7x2x1 mesh is selected to preserve reciprocal spacing relative to the previous 12x2x1 mesh; this is not a convergence study.\n\nRun only after the generated manifests, PAW hashes, resources, disk headroom and absence of competing compute jobs pass the supervisor preflight.\n''')
print(json.dumps(audit,indent=2))
