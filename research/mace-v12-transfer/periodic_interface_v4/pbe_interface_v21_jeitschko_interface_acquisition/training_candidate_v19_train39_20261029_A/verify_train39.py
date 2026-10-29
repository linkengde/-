#!/usr/bin/env python3
"""Verify all 80 epochs and evaluate the fixed epoch-79 model on frozen dev3 only."""
from __future__ import annotations
import argparse, hashlib, json, re
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import torch
from ase.io import read
from mace.calculators import MACECalculator

HERE = Path(__file__).resolve().parent
TASK = HERE.parent
TAG = 'MACE_periodic_v19_residual_response_run-45'
def digest(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()
def main() -> None:
    ap=argparse.ArgumentParser()
    ap.add_argument('--run-dir',type=Path,default=HERE/'run')
    args=ap.parse_args()
    run=args.run_dir.resolve()
    log=(run/f'logs/{TAG}.log').read_text()
    stdout=(run/'train.stdout').read_text()
    checkpoint_dir=run/'checkpoints'
    checkpoints={int(m.group(1)):p for p in checkpoint_dir.glob(f'{TAG}_epoch-*.pt') if (m:=re.fullmatch(rf'{re.escape(TAG)}_epoch-(\d+)\.pt',p.name))}
    epoch_lines=[int(x) for x in re.findall(r'INFO: Epoch (\d+):',log)]
    metric_pattern=re.compile(r'INFO: Epoch (\d+): head: Default, loss=\s*([0-9.eE+-]+), RMSE_E_per_atom=\s*([0-9.eE+-]+) meV, RMSE_F=\s*([0-9.eE+-]+) meV / A')
    metrics=[{'epoch':int(m.group(1)),'validation_loss':float(m.group(2)),'validation_energy_RMSE_meV_per_atom':float(m.group(3)),'validation_force_component_RMSE_meV_A':float(m.group(4))} for m in metric_pattern.finditer(log)]
    marker='Loaded Stage one model from epoch 79 for evaluation'
    receipt=json.loads((run/'launcher_completion.json').read_text())
    prior_report_path=run/'results/train39_verification.json'
    prior_report=json.loads(prior_report_path.read_text()) if prior_report_path.is_file() else None
    if epoch_lines!=list(range(80)): raise SystemExit(f'Expected epochs 0-79, got {epoch_lines}')
    if [x['epoch'] for x in metrics]!=list(range(80)): raise SystemExit('Validation metrics do not cover all 80 epochs')
    observed_checkpoint_epochs=sorted(checkpoints)
    full_checkpoint_set=observed_checkpoint_epochs==list(range(80))
    archived_pruned_set=(observed_checkpoint_epochs==[79] and prior_report is not None and
        prior_report.get('verification',{}).get('epoch_checkpoint_count')==80 and
        prior_report.get('verification',{}).get('epoch_checkpoint_coverage_pass') is True)
    if not (full_checkpoint_set or archived_pruned_set): raise SystemExit(f'Expected all checkpoints or verified archived epoch 79 only, got {observed_checkpoint_epochs}')
    if marker not in stdout: raise SystemExit('Export was not loaded from epoch 79')
    if receipt['training_command_exit_code']!=0 or not receipt['export_epoch79_marker_verified']: raise SystemExit('Training receipt failed')
    epoch79=checkpoints[79]
    if archived_pruned_set and digest(epoch79)!=prior_report['verification']['epoch79_checkpoint_sha256']:
        raise SystemExit('Retained epoch-79 checkpoint hash differs from the pre-pruning verification report')
    exported=run/'models/MACE_periodic_v19_residual_response.model'
    checkpoint_export=checkpoint_dir/f'{TAG}.model'
    if not (exported.is_file() and checkpoint_export.is_file()): raise SystemExit('Selected/final model files are missing')
    torch.set_default_dtype(torch.float64)
    ckpt=torch.load(epoch79,map_location='cpu',weights_only=False)
    model=torch.load(exported,map_location='cpu',weights_only=False).double().eval()
    state=model.state_dict()
    if set(ckpt['model'])!=set(state): raise SystemExit('Epoch-79 checkpoint and exported model tensor keys differ')
    differing=[k for k in state if not torch.equal(ckpt['model'][k].cpu(),state[k].cpu())]
    if differing: raise SystemExit(f'Exported model differs from epoch-79 checkpoint: {differing[:10]}')

    data=HERE/'dataset'
    dev=read(data/'valid.extxyz',index=':')
    baseline=json.loads((HERE/'baseline_v19_evaluation.json').read_text())
    baseline_rows={r['config_type']:r for r in baseline['per_interface']}
    calculator=MACECalculator(models=[model],device='cpu',default_dtype='float64')
    rows=[]; species_errors={}; aggregate_sq=0.0; aggregate_atoms=0; sep_errors=[]
    for atoms in dev:
        ref_e=float(atoms.info['REF_energy'])
        ref_f=np.asarray(atoms.arrays['REF_forces'],dtype=float)
        atoms.calc=calculator
        pred_e=float(atoms.get_potential_energy())
        pred_f=np.asarray(atoms.get_forces(),dtype=float)
        delta=pred_f-ref_f
        norms=np.linalg.norm(delta,axis=1)
        vector_rmse=float(np.sqrt(np.mean(norms**2)))
        component_rmse=float(np.sqrt(np.mean(delta**2)))
        aggregate_sq+=float(np.sum(delta**2));aggregate_atoms+=len(atoms)
        symbols=np.asarray(atoms.get_chemical_symbols())
        for element in sorted(set(symbols.tolist())):
            species_errors.setdefault(element,[]).extend(norms[symbols==element].tolist())
        central=np.flatnonzero(np.asarray(atoms.arrays.get('central_pair',np.zeros(len(atoms))),dtype=bool))
        if len(central)!=2: raise SystemExit(f'Missing central-pair marker: {atoms.info.get("config_type")}')
        ia,ib=map(int,central)
        vec=atoms.get_distance(ia,ib,mic=True,vector=True)
        dist=float(np.linalg.norm(vec))
        cfg=atoms.info['config_type']
        if cfg not in baseline_rows or not np.isclose(dist,baseline_rows[cfg]['central_pair_distance_A'],rtol=0,atol=1e-8):
            raise SystemExit(f'Development geometry/contact differs from frozen baseline: {cfg}')
        unit=vec/dist
        ref_sep=float(np.dot(ref_f[ib]-ref_f[ia],unit))
        pred_sep=float(np.dot(pred_f[ib]-pred_f[ia],unit))
        sep_error=abs(pred_sep-ref_sep);sep_errors.append(sep_error)
        row={
          'config_type':cfg,'interface':cfg.split('_v19_',1)[0],'atoms':len(atoms),
          'energy_error_eV_per_atom':(pred_e-ref_e)/len(atoms),
          'energy_abs_error_meV_per_atom':abs(pred_e-ref_e)*1000.0/len(atoms),
          'force_vector_RMSE_eV_A':vector_rmse,'force_component_RMSE_eV_A':component_rmse,
          'force_vector_max_error_eV_A':float(norms.max()),
          'central_pair_distance_A':dist,'separation_force_reference_eV_A':ref_sep,
          'separation_force_prediction_eV_A':pred_sep,'separation_force_abs_error_eV_A':sep_error,
        }
        old=baseline_rows[cfg]
        row['baseline_v19_force_vector_RMSE_eV_A']=old['force_vector_RMSE_eV_A']
        row['baseline_v19_separation_force_abs_error_eV_A']=old['separation_force_abs_error_eV_A']
        row['force_vector_RMSE_change_vs_v19_eV_A']=vector_rmse-old['force_vector_RMSE_eV_A']
        row['separation_force_error_change_vs_v19_eV_A']=sep_error-old['separation_force_abs_error_eV_A']
        rows.append(row)
    for element,values in species_errors.items(): species_errors[element]=float(np.sqrt(np.mean(np.square(values))))
    force_threshold=.05;sep_threshold=.10
    force_pass=all(r['force_vector_RMSE_eV_A']<=force_threshold for r in rows)
    sep_pass=all(x<=sep_threshold for x in sep_errors)
    dataset_manifest=json.loads((data/'manifest.json').read_text())
    result={
      'status':'TRAIN39_EXPERIMENT_DEV3_SCREEN_PASS' if force_pass and sep_pass else 'TRAIN39_EXPERIMENT_DEV3_SCREEN_FAIL',
      'verification':{
        'epoch_log_coverage':epoch_lines,'epoch_log_coverage_pass':True,
        'epoch_checkpoint_count':80,'epoch_checkpoint_coverage_pass':True,
        'retained_epoch_checkpoints':observed_checkpoint_epochs,
        'checkpoint_archive_pruned_after_full_verification':archived_pruned_set,
        'epoch79_export_tensor_identity_pass':not differing,
        'epoch79_checkpoint_sha256':digest(epoch79),'exported_model_sha256':digest(exported),
        'checkpoint_export_sha256':digest(checkpoint_export),'launcher_completion':receipt,
      },
      'inputs':{
        'training_candidate_manifest_sha256':digest(data/'manifest.json'),
        'train39_sha256':dataset_manifest['train_sha256'],
        'frozen_dev3_sha256':dataset_manifest['dev_sha256'],
        'foundation_model_sha256':dataset_manifest['foundation_model_sha256'],
        'test_or_sealed_data_opened':False,
      },
      'recipe':{'seed':45,'device':'cpu','batch_size':1,'epochs':80,'lr':1e-4,'weight_decay':5e-7,'energy_weight':100,'forces_weight':1000,'threads':4,'energy_key':'REF_energy','forces_key':'REF_forces'},
      'epoch_validation_metrics':metrics,
      'thresholds':{'per_interface_force_vector_RMSE_eV_A':.05,'per_interface_separation_force_abs_error_eV_A':.10},
      'aggregate_force_vector_RMSE_eV_A':float(np.sqrt(aggregate_sq/aggregate_atoms)),
      'per_species_force_vector_RMSE_eV_A':species_errors,
      'per_interface':rows,
      'gates':{'force_vector_pass':force_pass,'separation_force_pass':sep_pass,'overall_pass':force_pass and sep_pass},
      'limitations':[
        'This is an experiment-only train39 fit; official V19 train36/dev3 files were not changed.',
        'The three added training labels share one Jeitschko-Nwotny 1967 source family.',
        'The three frozen development labels share one COD 9009647 source family.',
        'Slab k-point, smearing, and vacuum convergence are not fully established.',
        'Passing these three development structures would not constitute blind test or production validation.'
      ]
    }
    receipt['status']=result['status']
    receipt['verification_completed_utc']=datetime.now(timezone.utc).isoformat().replace('+00:00','Z')
    receipt['development_gates']=result['gates']
    receipt['epoch79_export_tensor_identity_pass']=not differing
    (run/'launcher_completion.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    result['verification']['launcher_completion']=receipt
    out=run/'results/train39_verification.json'
    out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    curve=run/'results/train39_epoch_metrics.csv'
    curve.write_text('epoch,validation_loss,validation_energy_RMSE_meV_per_atom,validation_force_component_RMSE_meV_A\n'+''.join(f"{x['epoch']},{x['validation_loss']},{x['validation_energy_RMSE_meV_per_atom']},{x['validation_force_component_RMSE_meV_A']}\n" for x in metrics))
    print(json.dumps(result,indent=2))
if __name__=='__main__': main()
