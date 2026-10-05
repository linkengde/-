#!/usr/bin/env python3
"""Authorized post-evaluation CPU localization on exactly three already-scored probes."""
import argparse
import csv
import hashlib
import importlib.metadata
import json
from pathlib import Path
import numpy as np

PROJECT=Path('research/mace-v12-transfer/periodic_interface_v4')
REPORT=Path('research/mace-v12-transfer/coordination/reports/window-b/v15_force_error_localization')
CUTOFF=4.5
FORCE_TOL=1e-6
ENERGY_TOL_EV=1e-6

def require(ok,message):
    if not ok:raise ValueError(message)
def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def save(path,x):path.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')

def run(repo,out,inputs):
    from ase.io import read
    import torch
    from mace.calculators import MACECalculator
    project=repo/PROJECT;version=project/'mace_periodic_v15_interface_energy'
    def record(p):
        p=p.resolve();require(p.is_relative_to(repo),'Evidence path escapes repo');inputs[str(p.relative_to(repo))]=sha(p);return p
    def load(p):return json.loads(record(p).read_text())
    selection_path=version/'selected_model_record.json';selection=load(selection_path)
    evaluation_path=version/'evaluation_01/evaluation.json';evaluation=load(evaluation_path)
    require(evaluation['selection_record_sha256']==inputs[str(selection_path.relative_to(repo))],'Selection record differs from scored snapshot')
    require(selection['reviewed_by']=='window-a' and selection['completed_epochs']==80 and selection['training_exit_code']==0 and selection['blind_labels_used_for_selection'] is False,'Unreviewed/incomplete model selection')
    model=record(version/'training_01/checkpoints/MACE_periodic_v15_interface_energy_run-45.model')
    require(sha(model)==selection['model_sha256']==evaluation['model_sha256'],'Selected/scored model hash mismatch')
    expected={r['label']:r for r in evaluation['per_structure'] if r['role']=='fresh_local_registry'}
    labels={f'{x}_registry_holdout_v15_01' for x in ['AgC','AgSi','AgTi']};require(set(expected)==labels and len(expected)==3,'Expected exactly three already-scored probes')
    for line in record(version/'evaluation_01/SHA256SUMS.txt').read_text().splitlines():
        h,n=line.split('  ',1);require(Path(n).name==n,'Unsafe evaluation member');require(sha(version/'evaluation_01'/n)==h,'Published evaluation hash mismatch')
    previous=load(repo/'research/mace-v12-transfer/coordination/reports/window-b/v14_force_localization.json')
    record(repo/'research/mace-v12-transfer/coordination/reports/window-b/analyze_v14_force_localization.py')
    batch=load(repo/'research/mace-v12-transfer/coordination/reports/window-b/v15_acquisition_batch_priority/batch_recommendation.json')
    manifest=load(project/'pbe_interface_v15_registry_holdouts/input_manifest.json');meta={r['label']:r for r in manifest['records']}
    torch.set_default_dtype(torch.float64)
    loaded=torch.load(model,map_location='cpu',weights_only=False).double().eval()
    calc=MACECalculator(models=[loaded],device='cpu',default_dtype='float64')
    summaries=[];csv_rows=[]
    for label in sorted(labels):
        folder=project/'pbe_interface_v15_registry_holdouts/calculations'/label
        v=load(folder/'verification.json');summary=load(folder/'summary.json');progress=load(folder/'progress.json')
        require(v['status']=='PASS' and v['role']=='v15_blind_registry_holdout' and all(z is True for z in v['checks'].values()),'Unverified scored archive')
        require(v['label']==summary['label']==progress['label']==label and summary['scf_converged'] is True and progress['scf_converged'] is True,'Wrong/unconverged scored reference')
        require(v['input_sha256']==summary['source_sha256']==progress['source_sha256']==meta[label]['input_sha256'],'Input hash chain mismatch')
        for n,h in v['sha256'].items():
            require(Path(n).name==n,'Unsafe archive member');require(sha(record(folder/n))==h,'Reference archive member changed')
        path=record(folder/(label+'_PW_PBE.extxyz'));require(sha(path)==expected[label]['reference_sha256'],'Reference differs from original evaluation')
        inp=record(project/'pbe_interface_v15_registry_holdouts'/meta[label]['input']);require(sha(inp)==meta[label]['input_sha256'],'Input changed')
        a=read(path);original=read(inp)
        require(np.array_equal(a.numbers,original.numbers) and np.array_equal(a.pbc,original.pbc) and np.allclose(a.positions,original.positions,atol=1e-12,rtol=0) and np.allclose(a.cell,original.cell,atol=1e-12,rtol=0),'Reference geometry differs from input')
        for k in ['lammps_id','central_pair']:require(np.array_equal(a.arrays[k],original.arrays[k]),'ID/marked-pair mismatch')
        ids=a.arrays['lammps_id'].astype(int);require(len(set(ids))==len(a),'Duplicate IDs');symbols=np.array(a.get_chemical_symbols())
        fref=np.asarray(a.arrays['PW_PBE_forces'],dtype=float);eref=float(a.info['PW_PBE_energy_eV']);require(np.isfinite(fref).all() and np.isfinite(eref),'Nonfinite reference')
        predicted=a.copy();predicted.calc=calc;ep=float(predicted.get_potential_energy());fp=np.asarray(predicted.get_forces(),dtype=float)
        require(fp.shape==(len(a),3) and np.isfinite(fp).all() and np.isfinite(ep),'Invalid inference')
        delta=fp-fref;norm=np.linalg.norm(delta,axis=1);sse=norm**2
        dist=a.get_all_distances(mic=True);central=np.flatnonzero(a.arrays['central_pair']);require(len(central)==2,'Wrong marked pair')
        ia=next(int(i) for i in central if symbols[i]=='Ag');ix=next(int(i) for i in central if symbols[i]!='Ag')
        require(set(ids[central])==set(meta[label]['central_pair']['persistent_ids']),'Wrong marked IDs')
        vec=a.positions[ix]-a.positions[ia];u=vec/np.linalg.norm(vec)
        rref=float(np.dot(fref[ix]-fref[ia],u));rpred=float(np.dot(fp[ix]-fp[ia],u));signed=rpred-rref
        all_rmse=float(np.sqrt(np.mean(sse)));maxerr=float(norm.max());energy=abs(ep-eref)*1000/len(a)
        deltas={'force_vector_RMSE_eV_A':abs(all_rmse-expected[label]['force_vector_RMSE_eV_A']),'force_vector_max_error_eV_A':abs(maxerr-expected[label]['force_vector_max_error_eV_A']),'separating_force_abs_error_eV_A':abs(abs(signed)-expected[label]['separating_force_abs_error_eV_A']),'energy_prediction_eV_cell':abs(ep-expected[label]['energy_prediction_eV_cell'])}
        require(all(v <= (ENERGY_TOL_EV if k=='energy_prediction_eV_cell' else FORCE_TOL) for k,v in deltas.items()),'Published metric reproduction failed: '+json.dumps(deltas))
        local_distance=np.min(dist[:,central],axis=1);local=local_distance<=CUTOFF
        def region(mask):
            return {'atoms':int(mask.sum()),'force_vector_RMSE_eV_A':float(np.sqrt(np.mean(sse[mask]))) if mask.any() else None,'Cartesian_component_RMSE_xyz_eV_A':np.sqrt(np.mean(delta[mask]**2,axis=0)).tolist() if mask.any() else None,'force_error_SSE_share':float(sse[mask].sum()/sse.sum()) if sse.sum() else 0.}
        species={s:region(symbols==s) for s in sorted(set(symbols))}
        contributions=[]
        for i,sign in [(ia,-1),(ix,1)]:
            projected=float(np.dot(delta[i],u));transverse=delta[i]-projected*u
            contributions.append({'id':int(ids[i]),'element':str(symbols[i]),'force_error_vector_eV_A':delta[i].tolist(),'projected_error_along_Ag_to_X_eV_A':projected,'signed_contribution_to_separation_error_eV_A':sign*projected,'transverse_error_norm_eV_A':float(np.linalg.norm(transverse)),'total_vector_error_norm_eV_A':float(norm[i])})
        require(abs(sum(x['signed_contribution_to_separation_error_eV_A'] for x in contributions)-signed)<1e-10,'Projection decomposition mismatch')
        top=[]
        for i in np.argsort(norm)[::-1][:10]:
            neighbors=[int(j) for j in np.argsort(dist[i]) if j!=i][:4]
            top.append({'id':int(ids[i]),'index':int(i),'element':str(symbols[i]),'force_error_norm_eV_A':float(norm[i]),'force_error_vector_eV_A':delta[i].tolist(),'distance_to_marked_contact_A':float(local_distance[i]),'local':bool(local[i]),'neighbors':[{'id':int(ids[j]),'element':str(symbols[j]),'distance_A':float(dist[i,j])} for j in neighbors]})
        for i in range(len(a)):
            csv_rows.append({'label':label,'index':i,'id':int(ids[i]),'element':str(symbols[i]),'local_4p5A':bool(local[i]),'distance_to_marked_contact_A':float(local_distance[i]),**{f'DFT_F{axis}_eV_A':float(fref[i,j]) for j,axis in enumerate('xyz')},**{f'MACE_F{axis}_eV_A':float(fp[i,j]) for j,axis in enumerate('xyz')},**{f'error_F{axis}_eV_A':float(delta[i,j]) for j,axis in enumerate('xyz')},'force_error_norm_eV_A':float(norm[i])})
        old=next(x for x in previous['sources'] if x['label'].startswith(label.split('_')[0]+'_'))
        summaries.append({'label':label,'atoms':len(a),'reference_sha256':sha(path),'all':region(np.ones(len(a),dtype=bool)),'local':region(local),'outside':region(~local),'species':species,'max_force_error_eV_A':maxerr,'energy_error_meV_atom':energy,'marked_pair':{'Ag_id':int(ids[ia]),'X_id':int(ids[ix]),'unit_Ag_to_X_direct_cluster':u.tolist(),'DFT_separation_eV_A':rref,'MACE_separation_eV_A':rpred,'signed_error_eV_A':signed,'abs_error_eV_A':abs(signed),'atom_contributions':contributions},'top_10_atom_errors':top,'published_metric_abs_deltas':deltas,'V14_pattern_context':{'all_RMSE':old['force_vector_RMSE_all_eV_A'],'local_RMSE':old['force_vector_RMSE_within_4p5A_of_marked_pair_eV_A'],'outside_RMSE':old['force_vector_RMSE_outside_4p5A_eV_A'],'separation_abs_error':old['marked_pair']['absolute_error_eV_A'],'warning':'Different probes and trained models; not matched direct improvement/degradation evidence.'},'existing_candidate_ids':[x['candidate_id'] for x in batch['selected'] if x['interface']==label.split('_')[0]]})
    require(sha(model)==selection['model_sha256'] and sha(selection_path)==evaluation['selection_record_sha256'],'Frozen model/selection changed')
    report={'model_sha256':sha(model),'runtime_versions':{'torch':torch.__version__,'mace-torch':importlib.metadata.version('mace-torch'),'ase':importlib.metadata.version('ase')},'local_cutoff_A':CUTOFF,'local_geometry':'minimum-image distance to either marked atom; direct cluster vector for separation, consistent with original scorer','reproduction_tolerances':{'forces_eV_A':FORCE_TOL,'energy_eV_cell':ENERGY_TOL_EV,'justification':'CPU float64 and same selected model; small rounding/library-order tolerance, not gate relaxation'},'sources':summaries,'input_sha256':inputs,'scope':'Post-hoc localization of three already-scored V15 references; no new validation, model selection, training, DFT or MD/TTM. Scores remain FAIL; future reuse of scored probes is acquisition/history, not unseen validation.'}
    pooled={'force_vector_RMSE_eV_A':float(np.sqrt(sum(x['atoms']*x['all']['force_vector_RMSE_eV_A']**2 for x in summaries)/sum(x['atoms'] for x in summaries))),
            'energy_MAE_meV_atom':float(np.mean([x['energy_error_meV_atom'] for x in summaries])),
            'separating_force_max_abs_error_eV_A':max(x['marked_pair']['abs_error_eV_A'] for x in summaries),
            'force_vector_max_error_eV_A':max(x['max_force_error_eV_A'] for x in summaries)}
    target=evaluation['aggregate_by_role']['fresh_local_registry']
    differences={key:abs(value-target[key]) for key,value in pooled.items()}
    require(all(value<=1e-6 for value in differences.values()),'Pooled score reproduction failed')
    report['aggregate_metric_reproduction']={**pooled,'absolute_deltas':differences}
    save(out/'localization.json',report)
    with (out/'per_atom_errors.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(csv_rows[0]),lineterminator='\n');w.writeheader();w.writerows(csv_rows)
    print(json.dumps({'sources':[{k:x[k] for k in ['label','all','local','outside','species','marked_pair','top_10_atom_errors']} for x in summaries]},indent=2))

def main():
    p=argparse.ArgumentParser();p.add_argument('--repo-root',required=True,type=Path);p.add_argument('--output-dir',required=True,type=Path);a=p.parse_args();repo=a.repo_root.resolve();out=a.output_dir.resolve();require(out.is_relative_to(repo/REPORT) and out!=repo/REPORT,'Output must be a new report child');require(not out.exists() or not any(out.iterdir()),'Nonempty output refused');out.mkdir(parents=True,exist_ok=True);inputs={}
    try:run(repo,out,inputs)
    except Exception as e:
        save(out/'diagnostics_failure.json',{'error':type(e).__name__,'message':str(e),'input_sha256':inputs});raise
if __name__=='__main__':main()
