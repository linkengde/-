import csv, hashlib, importlib.util, importlib.metadata, json, time, subprocess, sys, os, tempfile
from pathlib import Path
import numpy as np
import torch
from ase.io import read
from ase.neighborlist import neighbor_list
from mace.calculators import MACECalculator

ROOT=Path(os.environ.get('CLOUD_C_REPO_ROOT',Path(__file__).resolve().parents[6])).resolve()
BASE=ROOT/'research/mace-v12-transfer'
OUT=Path(os.environ.get('CLOUD_C_AUDIT_OUTPUT',Path(tempfile.gettempdir())/'cloud-c-v18-audit')).resolve(); OUT.mkdir(parents=True,exist_ok=True)
torch.set_num_threads(4)
torch.set_num_interop_threads(1)
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ids(a):
    for key in ('lammps_id','local_new_id'):
        if key in a.arrays: return np.asarray(a.arrays[key], dtype=int)
    raise RuntimeError('Explicit atom identity missing')
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module
iface=load('iface',BASE/'coordination/reports/window-b/v19_existing_AgSi_candidate_force_error_screen_install_and_infer_window_b/infer_and_evaluate.py')
pure=load('pure',BASE/'coordination/reports/window-a/v19_V18_pure_phase_force_screen/infer_pure_phase.py')
lineage,_=iface.verify_model_lineage()
foundation=BASE/'periodic_interface_v4/mace_periodic_v12_interface_energy/foundation_models/mace-mp-0b3-medium.model'
text=(BASE/'file_manifest.json').read_text(); manifest,end=json.JSONDecoder().raw_decode(text)
assert text[end:].strip() in ('',r'\n')
assert sha(foundation)==manifest[str(foundation.relative_to(BASE))]['sha256']

labels=['Ag_baseline_k8x8x8_v19','Ti3SiC2_baseline_k10x10x2_v19','Ti3SiC2_baseline_k10x10x4_v19','Ti3SiC2_Ti4f_id0003_z_m020A_k10x10x4_v19','Ti3SiC2_Ti4f_id0003_z_p020A_k10x10x4_v19','Ti3SiC2_C4f_id0009_z_m020A_k10x10x4_v19','Ti3SiC2_C4f_id0009_z_p020A_k10x10x4_v19']
refs={}; missing=[]; archive_audit={}
for label in labels:
    folder=ROOT/pure.CONTROL_ROOT/label
    if not folder.is_dir(): missing.append(label); continue
    record=pure.verify_archive(label)
    atoms=record['atoms']; force=record['reference_forces']
    assert np.array_equal(ids(atoms),ids(read(ROOT/record['source_path'],index=0)))
    assert force.shape==(len(atoms),3) and np.isfinite(force).all()
    refs[label]=(atoms,force)
    archive_audit[label]={k:v for k,v in record.items() if k not in ('folder','atoms','reference_forces')}
for label,folder in iface.REFERENCES.items():
    record=iface.verify_archive(label,folder)
    atoms=record['_atoms']; force=np.asarray(atoms.arrays['PW_PBE_forces'],float)
    assert np.isfinite(force).all()
    refs[label]=(atoms,force)
    archive_audit[label]={k:v for k,v in record.items() if not k.startswith('_')}

models={'V18':ROOT/iface.MODEL_PATH,'foundation':foundation}
metrics=[]; atomrows=[]; speciesrows=[]; predictions={}; models_audit={}; consistency={}
def summarize(label,name,a,p,f):
    d=p-f; norms=np.linalg.norm(d,axis=1); sym=np.array(a.get_chemical_symbols())
    row={'label':label,'model':name,'n':len(a),'vector_RMSE_eV_A':float(np.sqrt(np.mean(norms**2))),'max_error_eV_A':float(norms.max()),'component_RMSE_eV_A':np.sqrt(np.mean(d*d,axis=0)).tolist(),'z_squared_error_fraction':float(np.sum(d[:,2]**2)/max(np.sum(d*d),1e-30)),'DFT_vector_RMS_eV_A':float(np.sqrt(np.mean(np.sum(f*f,axis=1))))}
    try: row['interface']=iface.force_metrics(p,f,a)
    except RuntimeError: pass
    metrics.append(row)
    for s in sorted(set(sym)):
        mask=sym==s
        speciesrows.append({'label':label,'model':name,'element':s,'n':int(mask.sum()),'vector_RMSE_eV_A':float(np.sqrt(np.mean(norms[mask]**2))),'x_RMSE_eV_A':float(np.sqrt(np.mean(d[mask,0]**2))),'y_RMSE_eV_A':float(np.sqrt(np.mean(d[mask,1]**2))),'z_RMSE_eV_A':float(np.sqrt(np.mean(d[mask,2]**2))),'signed_z_mean_error_eV_A':float(d[mask,2].mean()),'max_error_eV_A':float(norms[mask].max())})
    atom_ids=ids(a)
    for i in range(len(a)):
        ar={'label':label,'model':name,'index_zero_based':i,'atom_id':int(atom_ids[i]),'element':sym[i],'z_A':float(a.positions[i,2]),'error_norm_eV_A':float(norms[i])}
        for k,axis in enumerate('xyz'): ar.update({f'DFT_F{axis}_eV_A':float(f[i,k]),f'model_F{axis}_eV_A':float(p[i,k]),f'error_F{axis}_eV_A':float(d[i,k])})
        atomrows.append(ar)
for name,path in models.items():
    started=time.perf_counter(); calc=MACECalculator(model_paths=str(path),device='cpu',default_dtype='float64')
    net=calc.models[0]
    if name=='V18':
        checkpoint=torch.load(ROOT/iface.CHECKPOINT_PATH,map_location='cpu',weights_only=False)['model']
        state=net.state_dict()
        assert state.keys()==checkpoint.keys()
        assert all(torch.equal(value.cpu(),checkpoint[key].cpu()) for key,value in state.items())
        lineage['checks']['epoch79_state_tensors_rechecked_exactly']=True
        lineage['checks']['epoch79_state_tensor_count']=len(state)
    elements=[int(x) for x in net.atomic_numbers.cpu().tolist()]
    radius=float(net.r_max)
    assert {6,14,22,47}.issubset(elements)
    models_audit[name]={'path':str(path.relative_to(ROOT)),'sha256':sha(path),'atomic_numbers_order':elements,'requested_elements_present':True,'r_max_A':radius,'parameter_dtypes':sorted(set(str(x.dtype) for x in net.parameters())),'num_interactions':int(net.num_interactions)}
    models_audit[name]['scale_shift']={k:v.detach().cpu().tolist() for k,v in net.scale_shift.state_dict().items()}
    models_audit[name]['average_neighbors']=[float(x.avg_num_neighbors) for x in net.interactions]
    for label,(source,f) in refs.items():
        a=source.copy(); a.calc=calc; p=a.get_forces(); e=float(a.get_potential_energy())
        assert np.isfinite(p).all() and np.isfinite(e)
        predictions[name,label]=p.copy(); summarize(label,name,a,p,f)
        print(name,label,metrics[-1]['vector_RMSE_eV_A'],flush=True)
    label='Ti3SiC2_baseline_k10x10x4_v19'; a=refs[label][0].copy(); a.calc=calc
    p=a.get_forces(); original=a.positions.copy(); h=1e-4; energies=[]
    for sign in [-1,1]:
        a.positions[:]=original; a.positions[3,2]+=sign*h; energies.append(float(a.get_potential_energy()))
    fd=-(energies[1]-energies[0])/(2*h)
    a.positions[:]=original
    order=np.arange(len(a))[::-1]; perm=a[order]; perm.calc=calc
    perm_error=float(np.max(np.abs(perm.get_forces()-p[order])))
    translated=a.copy(); translated.positions += [0.123,0.234,0.345]; translated.calc=calc
    translation_error=float(np.max(np.abs(translated.get_forces()-p)))
    consistency[name]={'finite_difference_Fz_eV_A':fd,'autograd_Fz_eV_A':float(p[3,2]),'finite_difference_abs_error_eV_A':abs(fd-p[3,2]),'permutation_max_component_error_eV_A':perm_error,'translation_max_component_error_eV_A':translation_error}
    assert abs(fd-p[3,2])<1e-4 and perm_error<1e-8 and translation_error<1e-8
    models_audit[name]['elapsed_s']=time.perf_counter()-started

# Local geometry descriptors, from the nonsealed training file only.
train=read(ROOT/iface.TRAIN_PATH,index=':'); training=[]; coordrows=[]
for i,a in enumerate(train):
    sym=a.get_chemical_symbols()
    training.append({'frame':i,'config_type':a.info.get('config_type'),'formula':a.get_chemical_formula(),'n':len(a),'pbc':a.pbc.tolist(),'contains_Ag':'Ag' in sym,'positions_and_forces_finite':bool(np.isfinite(a.positions).all() and np.isfinite(a.arrays['REF_forces']).all()) if 'REF_forces' in a.arrays else bool(np.isfinite(a.positions).all())})
for group,structures in [('V18_train',train),('bulk_control',[refs['Ti3SiC2_baseline_k10x10x4_v19'][0]])]:
    for index,a in enumerate(structures):
        ii,jj,dd=neighbor_list('ijd',a,models_audit['V18']['r_max_A'])
        sym=np.array(a.get_chemical_symbols())
        for element in ['Ti','C','Si']:
            for atom in np.flatnonzero(sym==element):
                adjacent=jj[ii==atom]; distances=dd[ii==atom]
                coordrows.append({'group':group,'frame':index,'atom_index':int(atom),'element':element,'neighbors_at_model_cutoff':len(adjacent),'Ag_neighbors':int(np.sum(sym[adjacent]=='Ag')),'Ti_neighbors':int(np.sum(sym[adjacent]=='Ti')),'C_neighbors':int(np.sum(sym[adjacent]=='C')),'Si_neighbors':int(np.sum(sym[adjacent]=='Si')),'nearest_distance_A':float(distances.min()) if len(distances) else None})
# Same-geometry k-z difference compared to model residual, independently recomputed.
b2=refs['Ti3SiC2_baseline_k10x10x2_v19']; b4=refs['Ti3SiC2_baseline_k10x10x4_v19']
assert np.array_equal(b2[0].positions,b4[0].positions) and np.array_equal(b2[0].numbers,b4[0].numbers)
mesh_difference=float(np.sqrt(np.mean(np.sum((b2[1]-b4[1])**2,axis=1))))
responses=[]
for site,index in [('Ti4f_id0003',3),('C4f_id0009',9)]:
    minus=f'Ti3SiC2_{site}_z_m020A_k10x10x4_v19'; plus=f'Ti3SiC2_{site}_z_p020A_k10x10x4_v19'
    if minus not in refs or plus not in refs:
        responses.append({'site':site,'status':'PAIR_INCOMPLETE','missing':[x for x in [minus,plus] if x not in refs]}); continue
    # IDs, not filename integers, resolve the displaced atom.
    parent=b4[0]; am=refs[minus][0]; ap=refs[plus][0]
    changed=np.flatnonzero(np.linalg.norm(ap.positions-am.positions,axis=1)>1e-8)
    assert len(changed)==1; i=int(changed[0]); span=ap.positions[i,2]-am.positions[i,2]; assert abs(span-0.04)<1e-8
    assert np.array_equal(am.numbers,ap.numbers) and np.array_equal(am.cell.array,ap.cell.array)
    assert np.array_equal(ids(am),ids(ap)) and np.array_equal(ids(am),ids(parent))
    row={'site':site,'atom_index':i,'atom_id':int(ids(parent)[i]),'span_A':float(span),'DFT_restoring_slope_eV_A2':float(-(refs[plus][1][i,2]-refs[minus][1][i,2])/span)}
    for name in models: row[name+'_restoring_slope_eV_A2']=float(-(predictions[name,plus][i,2]-predictions[name,minus][i,2])/span)
    responses.append(row)
for filename,rows in [('per_atom.csv',atomrows),('per_element.csv',speciesrows),('coordination.csv',coordrows)]:
    with (OUT/filename).open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
result={'commit':subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),'runtime':{'python':sys.version.split()[0],'packages':{n:importlib.metadata.version(n) for n in ['torch','mace-torch','ase','numpy','scipy']},'threads':torch.get_num_threads(),'dtype':'float64','units':{'positions':'Angstrom','forces':'eV/Angstrom'},'device':'cpu'},'model_lineage_checks':lineage['checks'],'models':models_audit,'archives':archive_audit,'missing_archives':missing,'metrics':metrics,'consistency':consistency,'training':training,'pure_Ti3SiC2_train_frames':sum('Ag' not in a.get_chemical_symbols() and set(a.get_chemical_symbols())=={'Ti','Si','C'} for a in train),'bulk_k10x10x2_to_k10x10x4_force_RMS_eV_A':mesh_difference,'signed_responses':responses,'sealed_files_read':False,'DFT_or_training_launched':False}
(OUT/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
print('AUDIT_COMPLETE',flush=True)
