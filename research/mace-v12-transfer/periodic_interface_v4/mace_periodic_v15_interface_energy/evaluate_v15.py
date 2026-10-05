#!/usr/bin/env python3
"""A runtime only: evaluate a frozen selected model; never select from blind scores."""
import argparse
from hashlib import sha256
import csv
import json
from pathlib import Path
import numpy as np
from preflight import checks, labels, sha, require, DRAFT

GATES={'energy_MAE_meV_atom':10.,'force_vector_RMSE_eV_A':.05,'separating_force_max_abs_error_eV_A':.10}

def aggregate(rows):
    metrics={'structures':len(rows),'atoms':sum(r['atoms'] for r in rows),
            'energy_MAE_meV_atom':float(np.mean([r['energy_abs_error_meV_atom'] for r in rows])),
            'force_vector_RMSE_eV_A':float(np.sqrt(sum(r['atoms']*r['force_vector_RMSE_eV_A']**2 for r in rows)/sum(r['atoms'] for r in rows))),
            'force_vector_max_error_eV_A':max(r['force_vector_max_error_eV_A'] for r in rows),
            'separating_force_max_abs_error_eV_A':max((r['separating_force_abs_error_eV_A'] for r in rows if r['separating_force_abs_error_eV_A'] is not None),default=None)}

    require(all(value is None or np.isfinite(value) for value in metrics.values()),'Nonfinite aggregate metric')
    return metrics

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo-root',required=True,type=Path);ap.add_argument('--output-dir',required=True,type=Path)
    ap.add_argument('--model',required=True,type=Path);ap.add_argument('--selection-record',required=True,type=Path)
    args=ap.parse_args();repo=args.repo_root.resolve();out=args.output_dir.resolve();model=args.model.resolve()
    require(out.is_relative_to(repo/DRAFT) and out!=repo/DRAFT,'Output must be V15 run child')
    require(not out.exists() or (out.is_dir() and not any(out.iterdir())),'Nonempty evaluation output refused')
    require(model.is_file(),'Selected model missing')
    selection_bytes=args.selection_record.read_bytes()
    selection_digest=sha256(selection_bytes).hexdigest()
    selection=json.loads(selection_bytes)
    require(selection['model_sha256']==sha(model) and Path(selection['model_path']).resolve()==model,'Selected model/hash changed')
    require(selection['training_exit_code']==0 and selection['completed_epochs']==80 and 1<=selection['selected_epoch']<=80,'A must inspect 80-epoch completion and selection')
    require(selection['selection_uses_only_training_validation'] is True and selection['blind_labels_used_for_selection'] is False,'Blind selection forbidden')
    require(selection['reviewed_by']=='window-a' and selection['selection_reason'].strip(),'A reviewed selection record required')
    for field in ['training_stdout','epoch_completion_evidence']:
        require(sha(Path(selection[field]['path']))==selection[field]['sha256'],'Completion evidence changed')
    # This is the first blind-label access; it occurs only after frozen selection checks.
    _,splits,blind=checks(repo,archives=True,include_references=True)
    import torch
    from mace.calculators import MACECalculator
    torch.set_default_dtype(torch.float64)
    loaded=torch.load(model,map_location='cpu',weights_only=False).double().eval()
    calculator=MACECalculator(models=[loaded],device='cpu',default_dtype='float64')
    references=[(a.info.get('config_type',f'test_{i}'),'historical_regression',a,*labels(a),None) for i,a in enumerate(splits['test'])]
    references.extend((name,'fresh_local_registry',a,e,f,h) for name,a,e,f,h in blind)
    rows=[]
    for name,role,a,eref,fref,h in references:
        predicted=a.copy();predicted.calc=calculator
        ep=float(predicted.get_potential_energy());fp=np.asarray(predicted.get_forces(),dtype=float)
        require(np.isfinite(ep) and fp.shape==fref.shape and np.isfinite(fp).all(),'Invalid model prediction: '+name)
        delta=fp-fref;central=np.flatnonzero(np.asarray(a.arrays.get('central_pair',np.zeros(len(a))),dtype=bool))
        radial_ref=radial_pred=radial_error=None
        if role=='fresh_local_registry':require(len(central)==2,'Blind marked pair missing')
        if len(central)==2:
            ag=[int(i) for i in central if a[i].symbol=='Ag'];contact=[int(i) for i in central if a[i].symbol!='Ag']
            require(len(ag)==len(contact)==1,'Marked pair must be Ag-X')
            i,j=ag[0],contact[0]
            # Preserve historical evaluator convention: direct cluster coordinate vector.
            v=a.positions[j]-a.positions[i];d=float(np.linalg.norm(v));require(d>0,'Zero pair distance');u=v/d
            radial_ref=float(np.dot(fref[j]-fref[i],u));radial_pred=float(np.dot(fp[j]-fp[i],u));radial_error=abs(radial_pred-radial_ref)
        rows.append({'label':name,'role':role,'atoms':len(a),'formula':a.get_chemical_formula(),'reference_sha256':h,
                     'energy_DFT_eV_cell':eref,'energy_prediction_eV_cell':ep,'energy_abs_error_meV_atom':abs(ep-eref)*1000/len(a),
                     'force_vector_RMSE_eV_A':float(np.sqrt(np.mean(np.sum(delta**2,axis=1)))),
                     'force_vector_max_error_eV_A':float(np.linalg.norm(delta,axis=1).max()),
                     'separating_force_reference_eV_A':radial_ref,'separating_force_prediction_eV_A':radial_pred,
                     'separating_force_abs_error_eV_A':radial_error})
        derived=[rows[-1][key] for key in ['energy_abs_error_meV_atom','force_vector_RMSE_eV_A','force_vector_max_error_eV_A','separating_force_abs_error_eV_A'] if rows[-1][key] is not None]
        require(all(np.isfinite(value) for value in derived),'Nonfinite derived metric: '+name)
    groups={role:aggregate([r for r in rows if r['role']==role]) for role in ['historical_regression','fresh_local_registry']}
    interfaces={}
    for kind in ['AgC','AgSi','AgTi']:
        subset=[r for r in rows if r['role']=='fresh_local_registry' and r['label'].startswith(kind)]
        require(len(subset)==1,'Missing/duplicate interface probe')
        metrics=aggregate(subset);metrics['provisional_gates_met']=all(metrics[k] is not None and metrics[k]<=v for k,v in GATES.items());interfaces[kind]=metrics
    # Even meeting narrow gates does not establish reliable transfer or production PASS.
    status='FAIL' if any(not x['provisional_gates_met'] for x in interfaces.values()) else 'UNDETERMINED'
    require(sha(model)==selection['model_sha256'],'Model changed during evaluation')
    require(sha(args.selection_record)==selection_digest,'Selection record changed during evaluation')
    report={'model_path':str(model),'model_sha256':sha(model),'selection_record_sha256':selection_digest,
            'thresholds':GATES,'per_structure':rows,'aggregate_by_role':groups,'by_interface':interfaces,
            'overall_status':status,'all_provisional_gates_met':all(x['provisional_gates_met'] for x in interfaces.values()),
            'limitations':'14 partial and 12 declared inherited rows lack complete original-run revalidation; physical-method consistency UNKNOWN. Historical test correlated; fresh registry probes share small-cluster motifs, not independent morphology/thermal validation. No MD/TTM authorization.'}
    out.mkdir(parents=True,exist_ok=True)
    (out/'evaluation.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    with (out/'per_structure.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    (out/'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in sorted(out.iterdir()) if p.is_file()))
    print(json.dumps({'overall_status':status,'by_interface':interfaces},indent=2))

if __name__=='__main__':main()
