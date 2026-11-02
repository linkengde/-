#!/usr/bin/env python3
"""Verify the fixed skeleton fit and score every held-out configuration."""
from __future__ import annotations
import hashlib, json, re
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import torch
from ase.io import read
from mace.calculators import MACECalculator

ROOT = Path(__file__).resolve().parent
RUN = ROOT / 'run'
DATA = ROOT / 'dataset'
TAG = 'MACE_TiSiC_skeleton_v25_run-45'
THRESHOLD = 0.07
STRICTER_REFERENCE_THRESHOLD = 0.05

def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main() -> None:
    log_path = RUN / f'logs/{TAG}.log'
    stdout_path = RUN / 'train.stdout'
    checkpoint_dir = RUN / 'checkpoints'
    selected_ckpt = checkpoint_dir / f'{TAG}.model'
    exported = RUN / 'models/MACE_TiSiC_skeleton_v25.model'
    receipt_path = RUN / 'launcher_completion.json'
    for path in (log_path, stdout_path, selected_ckpt, exported, receipt_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f'missing fit artifact: {path}')

    log = log_path.read_text()
    stdout = stdout_path.read_text()
    epochs = [int(x) for x in re.findall(r'INFO: Epoch (\d+):', log)]
    if epochs != list(range(80)):
        raise SystemExit(f'expected all epochs 0–79 exactly once, got {epochs}')
    marker = 'Loaded Stage one model from epoch 79 for evaluation'
    if marker not in stdout:
        raise SystemExit('export was not loaded from epoch 79')
    checkpoints = sorted(
        int(m.group(1)) for p in checkpoint_dir.glob('MACE_TiSiC_skeleton_v25_run-45_epoch-*.pt')
        if (m := re.fullmatch(r'MACE_TiSiC_skeleton_v25_run-45_epoch-(\d+)\.pt', p.name))
    )
    previous_path = RUN / 'results/skeleton_fit_verification.json'
    previous = json.loads(previous_path.read_text()) if previous_path.is_file() else {}
    epoch79_path = checkpoint_dir / f'{TAG}_epoch-79.pt'
    full_checkpoint_set = checkpoints == list(range(80))
    pruned_after_verified_full_set = (
        checkpoints == [79]
        and previous.get('verification', {}).get('epoch_checkpoint_coverage_pass') is True
        and previous.get('verification', {}).get('epoch_checkpoint_count') == 80
        and previous.get('verification', {}).get('epoch79_checkpoint_sha256') == digest(epoch79_path)
    )
    if not (full_checkpoint_set or pruned_after_verified_full_set):
        raise SystemExit(f'expected checkpoint coverage 0–79, got {checkpoints}')
    receipt = json.loads(receipt_path.read_text())
    if receipt.get('training_command_exit_code') != 0:
        raise SystemExit('training launcher receipt is not successful')

    torch.set_default_dtype(torch.float64)
    ckpt = torch.load(checkpoint_dir / f'{TAG}_epoch-79.pt', map_location='cpu', weights_only=False)
    selected = torch.load(selected_ckpt, map_location='cpu', weights_only=False).double().eval()
    model = torch.load(exported, map_location='cpu', weights_only=False).double().eval()
    state = model.state_dict()
    if set(ckpt['model']) != set(state):
        raise SystemExit('epoch-79 checkpoint and exported-model tensor keys differ')
    differing = [key for key in state if not torch.equal(ckpt['model'][key].cpu(), state[key].cpu())]
    if differing:
        raise SystemExit(f'exported model differs from epoch-79 checkpoint: {differing[:10]}')
    selected_state = selected.state_dict()
    if set(selected_state) != set(state):
        raise SystemExit('selected checkpoint and exported-model tensor keys differ')
    differing = [key for key in state if not torch.equal(selected_state[key].cpu(), state[key].cpu())]
    if differing:
        raise SystemExit(f'selected checkpoint differs from epoch-79 export: {differing[:10]}')

    dev = read(DATA / 'valid.extxyz', index=':')
    manifest = json.loads((DATA / 'manifest.json').read_text())
    if manifest.get('acceptance_standard', {}).get('per_configuration_force_vector_RMSE_max_eV_A') != THRESHOLD:
        raise SystemExit('dataset acceptance standard does not match the fixed 0.07 eV/A criterion')
    calculator = MACECalculator(models=[model], device='cpu', default_dtype='float64')
    rows = []
    for atoms in dev:
        ref_forces = np.asarray(atoms.arrays['REF_forces'], dtype=float)
        atoms.calc = calculator
        predicted_forces = np.asarray(atoms.get_forces(), dtype=float)
        error = predicted_forces - ref_forces
        norms = np.linalg.norm(error, axis=1)
        rows.append({
            'config_type': atoms.info['config_type'],
            'atom_count': len(atoms),
            'force_vector_RMSE_eV_A': float(np.sqrt(np.mean(norms ** 2))),
            'force_component_RMSE_eV_A': float(np.sqrt(np.mean(error ** 2))),
            'force_vector_max_error_eV_A': float(norms.max()),
            'pass_0p07_eV_A': bool(np.sqrt(np.mean(norms ** 2)) <= THRESHOLD),
            'pass_0p05_eV_A': bool(np.sqrt(np.mean(norms ** 2)) <= STRICTER_REFERENCE_THRESHOLD),
        })
    if len(rows) != 4:
        raise SystemExit(f'expected four fixed development configurations, got {len(rows)}')
    passed = all(row['pass_0p07_eV_A'] for row in rows)
    result = {
        'status': 'SKELETON_ACCEPTANCE_PASS' if passed else 'SKELETON_ACCEPTANCE_FAIL',
        'completed_utc': datetime.now(timezone.utc).isoformat().replace('+00:00', 'Z'),
        'role': 'non_official_tisic_skeleton_candidate',
        'training_recipe': {'seed': 45, 'epochs': 80, 'cpu_threads': 4, 'batch_size': 1,
                            'energy_weight': 100, 'force_weight': 1000},
        'acceptance_standard': {'per_configuration_force_vector_RMSE_eV_A': THRESHOLD},
        'additional_stricter_reference_screen': {
            'per_configuration_force_vector_RMSE_eV_A': STRICTER_REFERENCE_THRESHOLD,
            'status': 'PASS' if all(row['pass_0p05_eV_A'] for row in rows) else 'FAIL',
            'all_development_configurations_pass': all(row['pass_0p05_eV_A'] for row in rows),
        },
        'per_configuration': rows,
        'all_development_configurations_pass': passed,
        'aggregate_force_vector_RMSE_eV_A': float(np.sqrt(
            sum(row['force_vector_RMSE_eV_A'] ** 2 * row['atom_count'] for row in rows)
            / sum(row['atom_count'] for row in rows))),
        'verification': {
            'epoch_log_coverage_pass': True,
            'epoch_checkpoint_coverage_pass': True,
            'epoch_checkpoint_count': 80,
            'retained_epoch_checkpoints': checkpoints,
            'checkpoint_archive_pruned_after_full_verification': pruned_after_verified_full_set,
            'epoch79_export_marker_pass': True,
            'epoch79_checkpoint_export_tensor_identity_pass': True,
            'selected_checkpoint_sha256': digest(selected_ckpt),
            'epoch79_checkpoint_sha256': digest(checkpoint_dir / f'{TAG}_epoch-79.pt'),
            'exported_model_sha256': digest(exported),
            'dataset_manifest_sha256': digest(DATA / 'manifest.json'),
            'training_receipt': receipt,
        },
        'limitations': [
            'The standard is applied to four development configurations, which are not an independent blind test.',
            'Passing establishes only this skeleton-candidate acceptance; it does not establish transferability to Ag-containing interfaces or general production use.',
            'The 0.07 eV/A criterion is per-configuration force-vector RMSE, not a maximum atomwise error bound.',
        ],
    }
    (RUN / 'results/skeleton_fit_verification.json').write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    receipt.update({
        'status': result['status'],
        'fit_verification_status': result['status'],
        'verification_completed_utc': result['completed_utc'],
    })
    (RUN / 'launcher_completion.json').write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result, indent=2))

if __name__ == '__main__':
    main()
