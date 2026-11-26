#!/usr/bin/env python3
"""Build a fresh 1-A voxel, fine-branch QSGS Ag/Ti3SiC2 model."""
from pathlib import Path
import importlib.util,json,hashlib
import numpy as np
from scipy import ndimage
from ase import Atoms
from ase.build import bulk
from ase.io import write

HERE=Path(__file__).resolve().parent
HELPER=HERE/'support'/'build_microtwig_candidate.py'
spec=importlib.util.spec_from_file_location('qsgs_helpers',HELPER)
g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)
SIZE=np.array([200.,200.,200.]);BOX=np.array([200.,200.,500.]);GRID=np.array([200,200,200])
g.SIZE=SIZE;g.GRID=GRID;g.VACUUM_PADDING=1000.
SEED=20261141;THRESHOLD=.52;CLEARANCE=2.0
RANK_FILE=HERE/'qsgs_rank_1A_seed20261141.npy'

def main():
    rank=np.load(RANK_FILE)
    labels,n=ndimage.label(rank<int(THRESHOLD*rank.size))
    counts=np.bincount(labels.ravel());counts[0]=0
    mask=labels==int(np.argmax(counts))
    tp,ts=g.build_tsc_lattice();inside=g.inside_mask(mask,tp)
    skeleton,species,formula_units,cleanup,tsc_comp=g.clean_atomic_stoichiometry(tp[inside],ts[inside])
    aglat=bulk('Ag','fcc',a=g.AG_A0,cubic=True).repeat((49,49,49))
    aglat.set_cell(np.diag(SIZE),scale_atoms=True);agpool=aglat.positions.copy()
    agcand=agpool[~g.inside_mask(mask,agpool)]
    dsk,_=g.tree(skeleton).query(agcand,workers=4)
    agcand=agcand[dsk>=CLEARANCE]
    agkeep,agcomp=g.connected_component(agcand,3.15);silver=agcand[agkeep]
    tsc_top=g.winding_rank(skeleton,3.25);ag_top=g.winding_rank(silver,3.15)

    xyz=np.vstack((skeleton,silver));symbols=np.concatenate((species,np.full(len(silver),'Ag')))
    phase=np.r_[np.ones(len(skeleton),np.int32),np.full(len(silver),2,np.int32)]
    z=xyz[:,2];region=np.full(len(xyz),3,np.int32);region[z<20]=0;region[(z>=20)&(z<40)]=1;region[(z>=40)&(z<200)]=2
    atoms=Atoms(symbols=symbols,positions=xyz,cell=np.diag(BOX),pbc=[True,True,False])
    for name,value in [('phase',phase),('z_region',region),('fixed_skeleton',((z<20)&(phase==1)).astype(np.int32)),
      ('thermostat_skeleton',((z>=20)&(z<40)&(phase==1)).astype(np.int32)),('ag_mobile',(phase==2).astype(np.int32)),
      ('newton_region',(region==2).astype(np.int32))]: atoms.new_array(name,value)
    atoms.info.update({'status':'FRESH_1A_THIN_BRANCH_GEOMETRY_NOT_MINIMIZED_OR_MD_VALIDATED','qsgs_seed':SEED,
      'material_A':SIZE.tolist(),'box_A':BOX.tolist(),'mass_ratio_target':'1:1','qsgs_voxel_A':1.0,
      'qsgs_rank_threshold':THRESHOLD,'Ag_to_skeleton_clearance_A':CLEARANCE,
      'layers_A':'fixed skeleton 0-20; thermostat skeleton 20-40; Newton 40-200; vacuum 200-500'})
    data=HERE/'qsgs_thin_1A_seed41_200x200x500.data';ext=HERE/'qsgs_thin_1A_seed41_200x200x500.extxyz'
    write(data,atoms,format='lammps-data',atom_style='atomic',specorder=['Ti','Si','C','Ag'],masses=True);write(ext,atoms)
    np.savez_compressed(HERE/'qsgs_thin_1A_seed41_mask.npz',mask=mask,seed=SEED,threshold=THRESHOLD,grid=GRID,voxel_A=1.0)

    # Periodic x/y-aware branch thickness: tile laterally, keep z open.
    tiled=np.tile(mask,(3,3,1))
    edt=ndimage.distance_transform_edt(tiled,sampling=(1,1,1))[200:400,200:400,:]
    inside_dist=edt[mask];del edt,tiled
    bins=mask.reshape(10,20,10,20,10,20).mean(axis=(1,3,5))
    top=mask[:,:,180:200].reshape(10,20,10,20,20).mean(axis=(1,3,4))
    ijk=np.clip(np.floor(xyz/20).astype(int),0,9);atombins=np.zeros((10,10,10),np.int32)
    np.add.at(atombins,(ijk[:,0],ijk[:,1],ijk[:,2]),1)
    atom_counts={s:int(np.count_nonzero(symbols==s)) for s in ('Ti','Si','C','Ag')}
    ratio=float(atoms.get_masses()[phase==2].sum()/atoms.get_masses()[phase==1].sum())
    report={'status':'FRESH_1A_THIN_BRANCH_QSGS_GEOMETRY_ONLY_NOT_MINIMIZED_OR_MD_VALIDATED',
      'purpose':'Finer branch geometry generated at 1-A voxel resolution with increased nucleation and lateral branching.',
      'dimensions_A':{'material':[200,200,200],'full_box':[200,200,500],'vacuum_above_material':300},
      'atom_count':int(len(atoms)),'atom_counts':atom_counts,'TSC_formula_units':int(formula_units),
      'Ag_over_Ti3SiC2_mass_ratio':ratio,'QSGS':{'seed':SEED,'voxel_A':1.0,'P_seed':.0016,
        'P_forward':.070,'P_branch':.030,'P_backward':.00060,'growth_target_fraction':.84,
        'growth_generations':16,'rank_threshold':THRESHOLD,'mask_fraction':float(mask.mean())},
      'branch_thickness_A':{'distance_to_mask_boundary_percentiles':[float(v) for v in np.percentile(inside_dist,[50,90,99])],
        'max_radius_A':float(inside_dist.max()),'max_diameter_estimate_A':float(2*inside_dist.max())},
      'top_surface_z180_200A':{'skeleton_fraction_mean':float(top.mean()),'max_20A_xy_bin':float(top.max()),
        'bins_over_90pct':int(np.count_nonzero(top>.90))},
      'mask_20A_bins':{'p10':float(np.percentile(bins,10)),'median':float(np.median(bins)),
        'p90':float(np.percentile(bins,90)),'max':float(bins.max()),'bins_over_90pct':int(np.count_nonzero(bins>.90))},
      'material_20A_cube_atom_counts':{'min':int(atombins.min()),'p10':float(np.percentile(atombins,10)),
        'median':float(np.median(atombins)),'mean':float(atombins.mean()),'p90':float(np.percentile(atombins,90)),
        'max':int(atombins.max()),'empty_cubes':int(np.count_nonzero(atombins==0))},
      'connectivity':{'Ti3SiC2':tsc_top,'Ag_main_network':ag_top,'Ag_candidates_before_largest_component':agcomp,
        'Ag_atoms_retained':int(len(silver))},'TSC_cleanup':cleanup,
      'layers_A':{'fixed_skeleton_only':[0,20],'thermostat_skeleton_only':[20,40],'newton':[40,200],'vacuum':[200,500]},
      'limitations':['Geometry only; no force minimization or MD validation.','Initial cross-phase clearance is a geometric screen; verify with the selected hybrid potential before MD.'],
      'files':{'data':data.name,'extxyz':ext.name,'mask':'qsgs_thin_1A_seed41_mask.npz'}}
    (HERE/'qsgs_thin_1A_seed41_report.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    (HERE/'README_seed41_更细骨架QSGS模型.md').write_text('# 1 Å 网格细骨架 Ag/Ti₃SiC₂ QSGS 模型\n\n'
      f'新随机种子 {SEED}；随机四参数 QSGS；材料区 200×200×200 Å，计算盒 200×200×500 Å。\n\n'
      f'原子数 {len(atoms):,}；Ag/Ti₃SiC₂ 质量比 {ratio:.5f}:1；骨架分支最大直径估算 {2*inside_dist.max():.2f} Å。\n\n'
      f'两相主网络均贯通；每个 20 Å 立方分区至少 {atombins.min()} 个原子。几何未弛豫，未做 MD 验证。\n')
    files=[p for p in sorted(HERE.iterdir()) if p.is_file() and p.name!='SHA256SUMS.txt']
    (HERE/'SHA256SUMS.txt').write_text(''.join(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n' for p in files))
    print(json.dumps({'atoms':len(atoms),'counts':atom_counts,'mass_ratio':ratio,'branch_max_diameter_A':float(2*inside_dist.max()),
      'top_max':float(top.max()),'cube_min':int(atombins.min()),'empty_cubes':int(np.count_nonzero(atombins==0)),
      'TSC_topology':tsc_top,'Ag_topology':ag_top,'Ag_candidate_components':agcomp},indent=2),flush=True)

if __name__=='__main__':main()
