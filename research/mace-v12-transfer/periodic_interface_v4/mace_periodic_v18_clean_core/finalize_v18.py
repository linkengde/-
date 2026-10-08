#!/usr/bin/env python3
"""Freeze epoch79 and score only the reused development set after runner success."""
import argparse
import fcntl
import hashlib
import importlib.util
import json
import os
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def publish_completed(repo, entry, run, out, summary):
    """Publish only frozen final artifacts; never stage all epoch checkpoints."""
    coordination = repo / 'research/mace-v12-transfer/coordination'
    spec = importlib.util.spec_from_file_location('coordination_sync', coordination / 'sync_tasks.py')
    sync = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(sync)
    def git(*args):
        return subprocess.run(['git', *args], cwd=repo, text=True, capture_output=True, check=True)
    paths = [run / 'train.stdout', run / 'selection_record.json', run / 'training_summary.json',
             run / 'artifact_manifest.json', run / 'evaluation.stdout',
             run / 'models/MACE_periodic_v18_clean_core.model',
             run / 'checkpoints/MACE_periodic_v18_clean_core_run-45_epoch-79.pt']
    paths += sorted((run / 'logs').glob('*.log'))
    paths += sorted((run / 'results').glob('*.txt'))
    paths += sorted(p for p in out.iterdir() if p.is_file())
    for p in paths:
        if not p.is_relative_to(entry) or not p.is_file() or p.stat().st_size >= 100_000_000:
            raise ValueError(f'Unsafe/oversized publication artifact: {p}')
    with sync.LOCK.open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if git('diff', '--name-only').stdout.strip() or git('diff', '--cached', '--name-only').stdout.strip():
            raise RuntimeError('Preserve unrelated tracked/index changes before publication; frozen results remain local')
        task_path = coordination / 'tasks/window-a.json'
        task = json.loads(task_path.read_text())
        sync.ensure_owner('window-a', task)
        task['current_stage'] = f"V18 clean30 completed80 epochs and frozen epoch79; development screen {summary['development_status']}. Final model/checkpoint, full logs and per-atom metrics published. Withheld labels unopened; B V18 review is unblocked."
        task['current_jobs']['V18_clean_core_training'] = f"completed; epochs0-79, fixed epoch79 export and data/model hashes verified; actual{summary['actual_cpu_threads']} CPU threads."
        task['current_jobs']['V18_development_evaluation'] = 'completed; ' + summary['development_status'] + '; same3 reused development structures; no test/withheld labels.'
        task['next_action'] = 'Read B V18 matched/cumulative review and targeted proposals before freezing the next dataset/recipe. Preserve gates and independent withheld labels; no production MD/TTM from provisional screening.'
        task['updated_utc'] = datetime.now(timezone.utc).isoformat()
        task_path.write_text(json.dumps(task, indent=2, ensure_ascii=False) + '\n')
        log_path = repo / 'research/mace-v12-transfer/WORK_LOG_2026-10.md'
        note = '\n\n### V18 completed controlled screening (UTC ' + datetime.now(timezone.utc).isoformat() + ')\n\n'
        note += '- Clean30 fixed-recipe run completed80 epochs; selected/exported epoch79 verified. Frozen same3 development screen: ' + summary['development_status'] + '. No withheld labels opened.\n'
        for name, metrics in summary['by_interface'].items():
            note += f"- {name}: energy {metrics['energy_abs_error_meV_atom']:.6f} meV/atom; vector RMSE {metrics['force_vector_RMSE_eV_A']:.6f} eV/A; separation error {metrics['separating_force_abs_error_eV_A']:.6f} eV/A.\n"
        note += '- Final model/checkpoint, selection, full logs, artifact hash manifest and per-atom evaluation saved; only selected checkpoint is published. B V18 review dependency satisfied. This remains a reused-development/provisional comparison with coverage and frame0 mesh confounds.\n'
        log_path.write_text(log_path.read_text() + note)
        paths += [task_path, log_path]
        git('add', '--', *[str(p.relative_to(repo)) for p in paths])
        git('commit', '-m', 'Publish completed V18 fixed-epoch model and development vectors')
    # The coordination command owns its lock and retries an ordinary push/rebase.
    # Conflicts stop here with final results/commits preserved; never force/reset.
    subprocess.run([sys.executable, str(coordination / 'sync_tasks.py'), 'progress', 'window-a',
                    '--job', 'V18_clean_core_training', '--state', 'completed', '--iteration', '80',
                    '--note', 'V18 completed80 epochs; fixed epoch79/data/model checks PASS; development status ' + summary['development_status'] + '. Final model/checkpoint/logs/per-atom metrics published; B review unblocked; sealed labels unopened.'], cwd=repo, check=True)


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
    summary = {'status': 'COMPLETED_PROVISIONAL_SCREENING', 'development_status': ev['status'], 'selected_epoch': 79, 'completed_epochs': 80, 'model_sha256': record['model_sha256'], 'actual_cpu_threads': record['actual_cpu_threads'], 'by_interface': ev['by_interface'], 'withheld_labels_opened': False}
    (run / 'training_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    files = [log, model, selection, run / 'evaluation.stdout', run / 'training_summary.json']
    files += sorted((run / 'checkpoints').glob('*.pt'))
    files += sorted((run / 'logs').glob('*.log'))
    files += sorted((run / 'results').glob('*.txt'))
    files += sorted(p for p in out.iterdir() if p.is_file())
    files += [entry / n for n in ['run_v18_training.sh', 'build_v18_dataset.py', 'evaluate_v18_dev.py', 'finalize_v18.py', 'data/train.extxyz', 'data/valid.extxyz', 'data/dataset_manifest.json']]
    (run / 'artifact_manifest.json').write_text(json.dumps({'checkpoint_policy': 'Intermediate checkpoints remain local; only epoch79 is published.', 'files': [{'path': str(p.relative_to(repo)), 'bytes': p.stat().st_size, 'sha256': sha(p), 'published': p.parent != run / 'checkpoints' or p.name.endswith('_epoch-79.pt')} for p in files]}, indent=2) + '\n')
    publish_completed(repo, entry, run, out, summary)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
