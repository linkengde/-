from pathlib import Path
import csv, json, time, hashlib, gc, os, resource
import numpy as np
import torch
from ase.io import read, write
from ase import units
from ase.constraints import FixAtoms
from ase.md.nose_hoover_chain import NoseHooverChainNVT
from mace.calculators import MACECalculator
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[5] if False else None
BASE=ROOT.parents[5]
# Script lives under reports/window-f/<task>; locate repo root directly from current cwd.
REPO=Path.cwd()
SOURCE=REPO/'research/mace-v12-transfer/coordination/reports/window-f/bicontinuous_ag_ti3sic2_ttm_20261015/exploratory_melt_cool_20261018/input/NEW_MODEL_ATOMS.data'
COMMON=REPO/'research/mace-v12-transfer/coordination/reports/window-f/bicontinuous_ag_ti3sic2_ttm_20261015/exploratory_melt_cool_20261018/timestep_diagnostic_20261018/common_initial_velocities.csv'
MODELS={
 'MACE_MP': REPO/'research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v12_interface_energy/foundation_models/mace-mp-0b3-medium.model',
 'MACE_V17': REPO/'research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v17_interface_energy/training_provisional_03/models/MACE_periodic_v17_provisional_screening.model',
}
TYPE_TO_SYMBOL={1:'Ti',2:'Si',3:'C',4:'Ag'}
ELEMENT_TO_Z={'Ti':22,'Si':14,'C':6,'Ag':47}
N_STEPS=int(os.environ.get('MACE_LOCAL_STEPS','500'))
OUTPUT_EVERY=5
DT=0.2*units.fs
T_TARGET=300.0
TDAMP=100.0*units.fs
R_ACTIVE=8.0
R_HALO=16.0


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def min_image_delta(a,b,cell):
    d=a-b
    for j in (1,2): d[j]-=cell[j,j]*np.rint(d[j]/cell[j,j])
    return d

def write_patch():
    src=read(SOURCE,format='lammps-data',atom_style='atomic',Z_of_type={1:22,2:14,3:6,4:47})
    src.pbc=[False,True,True]
    ids=np.asarray(src.arrays['id'],dtype=int)
    id_index={int(aid):i for i,aid in enumerate(ids)}
    center_i=id_index[19219]; center=src.positions[center_i]
    ds=src.positions-center
    L=src.cell.lengths()
    for j in (1,2): ds[:,j]-=L[j]*np.rint(ds[:,j]/L[j])
    radii=np.linalg.norm(ds,axis=1)
    mask=radii<=R_HALO
    patch=src[mask]
    patch_ids=np.asarray(patch.arrays['id'],dtype=int)
    dr=ds[mask]; rr=np.linalg.norm(dr,axis=1)
    mobile=rr<=R_ACTIVE
    # Existing common 300 K velocity sample is retained for IDs in the local patch;
    # frozen buffer velocities are set to zero by the positional constraint.
    vels={}
    with COMMON.open() as f:
        for r in csv.DictReader(f): vels[int(r['id'])]=np.array([float(r['vx_A_per_ps']),float(r['vy_A_per_ps']),float(r['vz_A_per_ps'])])
    v=np.array([vels[int(aid)] for aid in patch_ids],dtype=float)
    v[~mobile]=0.0
    if mobile.any():
        v[mobile]-=v[mobile].mean(axis=0)
        # Use original velocity direction/distribution, scaled to exactly 300 K over mobile DOFs.
        patch.set_velocities(v)
        patch.set_constraint(FixAtoms(indices=np.flatnonzero(~mobile)))
        dof=patch.get_number_of_degrees_of_freedom()
        T0=patch.get_temperature()
        if not np.isfinite(T0) or T0<=0: raise ValueError(f'bad local mobile temperature {T0}')
        v[mobile]*=(T_TARGET/T0)**0.5
        patch.set_velocities(v)
    patch.arrays['source_id']=patch_ids.copy()
    patch.arrays['mobile']=mobile.astype(np.int32)
    patch.arrays['radius_from_C19219_A']=rr.astype(float)
    write(ROOT/'LOCAL_PATCH_INITIAL.extxyz',patch,format='extxyz')
    with (ROOT/'LOCAL_PATCH_GROUPS.csv').open('w',newline='') as f:
        w=csv.writer(f); w.writerow(['source_id','element','role','distance_from_C19219_A','x_A','y_A','z_A'])
        for i,aid in enumerate(patch_ids):
            w.writerow([int(aid),patch[i].symbol,'mobile_r_le_8A' if mobile[i] else 'fixed_buffer_r_8_to_16A',f'{rr[i]:.10f}',*[f'{q:.10f}' for q in patch.positions[i]]])
    manifest={'source_data':str(SOURCE.relative_to(REPO)),'source_sha256':sha(SOURCE),'atom_count_full':len(src),'mobile_center_source_id':19219,
              'mobile_region':'sphere radius <= 8.0 A around C ID 19219, using X nonperiodic and Y/Z minimum image',
              'fixed_buffer':'remaining atoms within radius <= 16.0 A; MACE cutoff 6 A, gives >6 A static environment around the mobile boundary',
              'patch_atom_count':len(patch),'mobile_atom_count':int(mobile.sum()),'fixed_buffer_atom_count':int((~mobile).sum()),
              'boundary':'X nonperiodic, Y/Z periodic; original full cell retained',
              'starting_coordinates':'exactly from NEW_MODEL_ATOMS.data; no relaxation or deletion from source; atoms outside local 16 A patch omitted only for this explicitly local diagnostic',
              'initial_velocities':'reused ID-matched prior 300 K velocity sample for mobile atoms, zeroed for fixed buffer, mobile COM removed and scaled to 300 K',
              'mobile_initial_temperature_K':patch.get_temperature(),'mobile_degrees_of_freedom':dof,
              'note':'This is not a full-model MD result; frozen outer buffer limits heat/strain transport.'}
    (ROOT/'LOCAL_PATCH_AUDIT.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return patch,patch_ids,mobile,center

def pair_metrics(atoms, ids, mobile, center):
    pos=atoms.positions; types=np.array([ELEMENT_TO_Z.get(a.symbol,0) for a in atoms])
    cell=atoms.cell.array; L=atoms.cell.lengths()
    # Cutoff 3.4 A covers compressed Ti-C and Ag cross-contact audit.
    points=pos.copy()
    points[:,1]%=L[1];points[:,2]%=L[2]
    tree=cKDTree(points,boxsize=[1e6,L[1],L[2]])
    pp=tree.query_pairs(3.4,output_type='ndarray')
    if len(pp):
        dv=points[pp[:,0]]-points[pp[:,1]]
        for j in (1,2): dv[:,j]-=L[j]*np.rint(dv[:,j]/L[j])
        dd=np.linalg.norm(dv,axis=1)
    else: dd=np.empty(0)
    specs=[('Ti_C',22,6),('Ag_C',47,6),('Ag_Ti',47,22)]
    out={}
    for label,z1,z2 in specs:
        m=(((types[pp[:,0]]==z1)&(types[pp[:,1]]==z2))|((types[pp[:,0]]==z2)&(types[pp[:,1]]==z1))) if len(pp) else np.zeros(0,bool)
        m &= (mobile[pp[:,0]]|mobile[pp[:,1]]) if len(pp) else m
        if m.any(): out['min_'+label+'_A']=float(dd[m].min())
        else: out['min_'+label+'_A']=float('nan')
    id_map={int(aid):i for i,aid in enumerate(ids)}
    for label,aid,bid in [('C19219_Ti25627',19219,25627),('C19219_Ti25628',19219,25628),('C19219_Ag38116',19219,38116)]:
        if aid in id_map and bid in id_map:
            out[label+'_A']=float(np.linalg.norm(min_image_delta(pos[id_map[aid]],pos[id_map[bid]],cell)))
        else: out[label+'_A']=float('nan')
    return out,pp,dd

def run_model(model_key,model_path,patch_template,ids,mobile):
    case=ROOT/model_key; case.mkdir(exist_ok=True)
    atoms=patch_template.copy()
    atoms.set_constraint(FixAtoms(indices=np.flatnonzero(~mobile)))
    t0=time.time(); rss0=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024/1024
    calc=MACECalculator(model_paths=str(model_path),device='cpu',default_dtype='float64',enable_cueq=False)
    atoms.calc=calc
    e0=atoms.get_potential_energy(); f0=atoms.get_forces()
    init_forces={'id':ids[mobile].tolist(),'fnorm':np.linalg.norm(f0[mobile],axis=1)}
    print(model_key,'loaded; initial PE',e0,'initial active Fmax/P95/P99',float(init_forces['fnorm'].max()),float(np.percentile(init_forces['fnorm'],95)),float(np.percentile(init_forces['fnorm'],99)),'calc_seconds',time.time()-t0,flush=True)
    # Reinitialize same common velocities on the active group for each model, scale to 300 K.
    v=atoms.get_velocities(); v[~mobile]=0; v[mobile]-=v[mobile].mean(axis=0); atoms.set_velocities(v)
    T0=atoms.get_temperature(); v=atoms.get_velocities(); v[mobile]*=(300.0/T0)**0.5; atoms.set_velocities(v)
    dyn=NoseHooverChainNVT(atoms,timestep=DT,temperature_K=T_TARGET,tdamp=TDAMP)
    frames=case/f'{model_key}_{N_STEPS*0.2/1000.0:.4f}ps.extxyz'; metrics=[]
    source_id=ids.copy(); active_idx=np.flatnonzero(mobile)
    def sample(step):
        forces=atoms.get_forces(); fmag=np.linalg.norm(forces,axis=1)
        temp=float(atoms.get_temperature()); pe=float(atoms.get_potential_energy())
        pm,_,_=pair_metrics(atoms,source_id,mobile,None)
        row={'model':model_key,'step':step,'time_ps':step*0.2/1000.0,'active_N':len(active_idx),'T_mobile_K':temp,'PE_eV':pe,
             'Fmax_mobile_eV_A':float(fmag[mobile].max()),'P95_mobile_eV_A':float(np.percentile(fmag[mobile],95)),
             'P99_mobile_eV_A':float(np.percentile(fmag[mobile],99)),
             'Fmax_mobile_source_id':int(source_id[active_idx[np.argmax(fmag[mobile])]]),**pm}
        metrics.append(row)
        atoms.arrays['source_id']=source_id.copy(); atoms.arrays['mobile']=mobile.astype(np.int32)
        atoms.arrays['mace_force_eV_A']=forces.copy(); atoms.arrays['mace_force_norm_eV_A']=fmag.copy()
        atoms.info['model']=model_key; atoms.info['step']=step; atoms.info['time_ps']=row['time_ps']; atoms.info['T_mobile_K']=temp; atoms.info['PE_eV']=pe
        write(frames,atoms,format='extxyz',append=(step!=0))
        print(model_key,'step',step,'t_ps',row['time_ps'],'T',temp,'Fmax',row['Fmax_mobile_eV_A'],'minTiC',row['min_Ti_C_A'],flush=True)
        return row
    sample(0)
    stop_reason='completed_requested_steps'
    for step in range(1,N_STEPS+1):
        dyn.run(1)
        forces=atoms.get_forces(); fmag=np.linalg.norm(forces,axis=1)
        if np.any(~np.isfinite(fmag)):
            stop_reason=f'NaN_or_Inf_force_at_step_{step}'; break
        # Stop criteria for the local diagnostic: severe spike or Ti-C contact below 1.50 A.
        if fmag[mobile].max()>50:
            stop_reason=f'Fmax_mobile_gt_50_at_step_{step}'
            sample(step); break
        if step%OUTPUT_EVERY==0 or step==N_STEPS:
            row=sample(step)
            if row['min_Ti_C_A']<1.50:
                stop_reason=f'min_Ti_C_below_1.50A_at_step_{step}'; break
    with (case/'METRICS.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(metrics[0]));w.writeheader();w.writerows(metrics)
    # Save endpoint checkpoint in extxyz for local replay.
    atoms.arrays['source_id']=source_id.copy(); atoms.arrays['mobile']=mobile.astype(np.int32)
    atoms.arrays['mace_force_eV_A']=atoms.get_forces().copy()
    write(case/'FINAL_LOCAL_PATCH.extxyz',atoms,format='extxyz')
    report={'model':model_key,'model_file':str(model_path.relative_to(REPO)),'model_sha256':sha(model_path),'MACE_version':'0.3.16','torch_version':torch.__version__,'device':'CPU','default_dtype':'float64','n_steps_requested':N_STEPS,'n_steps_completed':int(metrics[-1]['step']),
            'step_fs':0.2,'elapsed_wall_seconds':time.time()-t0,'initial_active_temperature_K':300.0,'thermostat':'NoseHooverChainNVT; target 300 K; tdamp 100 fs',
            'mobile_atoms':int(mobile.sum()),'local_patch_atoms':len(atoms),'stop_reason':stop_reason,
            'peak_rss_gib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024/1024,'note':'Local fixed-buffer trajectory only; cannot be interpreted as full-model stability or melt-pool test.'}
    (case/'RUN_REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
    del atoms,calc,dyn; gc.collect(); return report

def main():
    torch.set_num_threads(4)
    patch,ids,mobile,center=write_patch()
    print('patch atoms',len(patch),'mobile',int(mobile.sum()),'buffer',int((~mobile).sum()),'source SHA',sha(SOURCE),flush=True)
    for key,path in MODELS.items():
        if not path.exists(): raise FileNotFoundError(path)
    keys=[x for x in os.environ.get('MACE_MODEL_KEYS','MACE_MP,MACE_V17').split(',') if x]
    reports=[]
    for key in keys:
        if key not in MODELS: raise ValueError(f'unknown model key {key}')
        reports.append(run_model(key,MODELS[key],patch,ids,mobile))
    (ROOT/'SUMMARY.json').write_text(json.dumps({'models':reports,'local_patch_audit':json.loads((ROOT/'LOCAL_PATCH_AUDIT.json').read_text())},indent=2)+'\n')

if __name__=='__main__': main()
