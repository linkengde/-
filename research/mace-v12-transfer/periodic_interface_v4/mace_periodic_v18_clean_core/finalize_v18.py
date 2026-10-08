#!/usr/bin/env python3
"""Freeze epoch79 and score only the reused development set after runner success."""
import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--run-dir', type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo_root.resolve()
    run = args.run_dir.resolve()
    entry = Path(__file__).resolve().parent
    if not run.is_relative_to(entry) or run == entry:
        raise ValueError('Unexpected run directory')
    model = run / 'models/MACE_periodic_v18_clean_core.model'
    log = run / 'train.stdout'
    text = log.read_text()
    if [int(x) for x in re.findall(r'INFO: Epoch (\d+):', text)] != list(range(80)):
        raise ValueError('Incomplete epoch log')
    for marker in ['INFO: Training complete', 'INFO: Done', 'Loaded Stage one model from epoch 79 for evaluation']:
        if marker not in text:
            raise ValueError(f'Missing completion evidence: {marker}')
    record = {
        'model_path': str(model.relative_to(repo)),
        'model_sha256': sha(model),
        'training_stdout': {'path': str(log.relative_to(repo)), 'sha256': sha(log)},
        'training_exit_code': 0, 'completed_epochs': 80, 'selected_epoch': 79,
        'selection_uses_only_training_validation': True,
        'blind_labels_used_for_selection': False,
        'test_file_passed_to_training': False,
        'reviewed_by': 'window-a',
        'actual_cpu_threads': int(os.environ['OMP_NUM_THREADS']),
        'completion_evidence': 'Invoked by set-euo-pipefail training runner after successful train/export guards; no standalone launcher exit receipt.',
        'selection_reason': 'Predeclared fixed final epoch79 for clean30-versus-provisional43 comparison; runner verified the complete training process and exact exported epoch.',
        'dataset_manifest_sha256': sha(entry / 'data/dataset_manifest.json'),
        'training_input_sha256': sha(entry / 'data/train.extxyz'),
        'development_validation_input_sha256': sha(entry / 'data/valid.extxyz'),
        'frozen_at_utc': datetime.now(timezone.utc).isoformat(),
    }
    selection = run / 'selection_record.json'
    if selection.exists():
        raise ValueError('Existing selection record refused')
    selection.write_text(json.dumps(record, indent=2) + '\n')
    out = entry / 'evaluation_dev_01'
    result = subprocess.run([sys.executable, str(entry / 'evaluate_v18_dev.py'), '--repo-root', str(repo), '--output-dir', str(out), '--model', str(model), '--selection-record', str(selection)], cwd=repo, capture_output=True, text=True)
    (run / 'evaluation.stdout').write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f'Development evaluation failed; inspect {run / "evaluation.stdout"}')
    ev = json.loads((out / 'evaluation.json').read_text())
    summary = {'status': 'COMPLETED_PROVISIONAL_SCREENING', 'development_status': ev['status'], 'selected_epoch': 79, 'completed_epochs': 80, 'model_sha256': record['model_sha256'], 'by_interface': ev['by_interface'], 'withheld_labels_opened': False}
    (run / 'training_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    files = [log, model, selection, run / 'evaluation.stdout', run / 'training_summary.json']
    files += sorted((run / 'checkpoints').glob('*.pt'))
    files += sorted((run / 'logs').glob('*.log'))
    files += sorted((run / 'results').glob('*.txt'))
    (run / 'artifact_manifest.json').write_text(json.dumps({'files': [{'path': str(p.relative_to(repo)), 'bytes': p.stat().st_size, 'sha256': sha(p)} for p in files]}, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
