#!/usr/bin/env python3
"""Claim and publish disjoint research tasks through ordinary Git pushes."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import uuid

HERE = Path(__file__).resolve().parent
REPO = Path(subprocess.check_output(['git', 'rev-parse', '--show-toplevel'], cwd=HERE, text=True).strip())
TASKS = HERE / 'tasks'
IDENTITY = Path('/workspace') / ('.mace-window-' + hashlib.sha256(str(REPO).encode()).hexdigest()[:12] + '.json')

def git(*args, check=True):
    return subprocess.run(['git', *args], cwd=REPO, text=True, capture_output=True, check=check)

def now():
    return datetime.now(timezone.utc).isoformat()

def identity():
    if not IDENTITY.exists():
        IDENTITY.write_text(json.dumps({'instance': str(uuid.uuid4()), 'created_utc': now()}, indent=2) + '\n')
    return json.loads(IDENTITY.read_text())['instance']

def task_path(name):
    if name not in ('window-a', 'window-b'):
        raise SystemExit('Choose window-a or window-b.')
    return TASKS / (name + '.json')

def load_task(name):
    return json.loads(task_path(name).read_text())

def save_task(name, task):
    task_path(name).write_text(json.dumps(task, ensure_ascii=False, indent=2) + '\n')

def ensure_owner(name, task):
    if task.get('owner_instance') != identity():
        raise SystemExit('Task belongs to another cloud instance. Do not start or repeat its calculations.')

def push_or_stop(name, retry):
    pushed = git('push', 'origin', 'HEAD:main', check=False)
    if pushed.returncode == 0:
        print(pushed.stdout + pushed.stderr)
        return
    if not retry:
        raise SystemExit('Claim was NOT published. Do not compute. Fetch main and inspect the winning owner.\n' + pushed.stderr)
    git('fetch', 'origin', 'main')
    rel = str(task_path(name).relative_to(REPO))
    remote = json.loads(git('show', 'origin/main:' + rel).stdout)
    ensure_owner(name, remote)
    rebased = git('rebase', 'origin/main', check=False)
    if rebased.returncode:
        raise SystemExit('Publication rebase has a conflict. Preserve results and resolve it before proceeding.\n' + rebased.stderr)
    pushed = git('push', 'origin', 'HEAD:main', check=False)
    if pushed.returncode:
        raise SystemExit('Publication is not confirmed. Preserve results and retry after fetching main.\n' + pushed.stderr)
    print(pushed.stdout + pushed.stderr)

def claim(name):
    # A failed ordinary push rejects simultaneous claims; no calculation starts first.
    if git('diff', '--name-only').stdout.strip() or git('diff', '--cached', '--name-only').stdout.strip():
        raise SystemExit('Publish or preserve existing tracked changes before claiming a task.')
    git('fetch', 'origin', 'main')
    merged = git('merge', '--ff-only', 'origin/main', check=False)
    if merged.returncode:
        raise SystemExit('Cannot fast-forward main. Reconcile local commits before claiming.\n' + merged.stderr)
    if git('rev-parse', 'HEAD').stdout != git('rev-parse', 'origin/main').stdout:
        raise SystemExit('Publish local commits before claiming; HEAD must equal origin/main.')
    task = load_task(name)
    if task['status'] == 'completed':
        raise SystemExit('Task is completed. Inspect archived results; do not repeat it.')
    if task.get('owner_instance'):
        ensure_owner(name, task)
        print('Task already belongs to this cloud instance:', name)
        return
    other = task.get('requires_separate_instance_from')
    if other and load_task(other).get('owner_instance') == identity():
        raise SystemExit('This task needs a separate cloud machine. Current instance already runs ' + other + '.')
    task.update(owner_instance=identity(), status='claimed', claimed_utc=now(), updated_utc=now())
    save_task(name, task)
    git('add', '--', str(task_path(name).relative_to(REPO)))
    git('commit', '-m', 'Claim ' + name + ' research task before computation')
    push_or_stop(name, retry=False)
    print('Claim published. Computation may start:', name)

def publish(name, label):
    task = load_task(name)
    ensure_owner(name, task)
    if label:
        if label not in task.get('labels', {}):
            raise SystemExit('Label is outside this task assignment: ' + label)
        folder = REPO / task['result_directory'] / 'calculations' / label
        summary = json.loads((folder / 'summary.json').read_text())
        if summary.get('scf_converged') is not True:
            raise SystemExit('Only converged, archive-verified labels may be marked completed.')
        if not (folder / 'verification.json').is_file():
            raise SystemExit('Run archive_parallel_pw.py before publishing this label.')
        verification = json.loads((folder / 'verification.json').read_text())
        if verification.get('status') != 'PASS':
            raise SystemExit('Archive verification did not pass.')
        task['labels'][label] = 'completed'
    labels_complete = task.get('labels') and all(v == 'completed' for v in task['labels'].values())
    task['status'] = 'completed' if labels_complete and task.get('auto_complete_after_labels', True) else 'running'
    task['updated_utc'] = now()
    prefixes = task['allowed_publish_paths']
    def allowed(path):
        return any(path == p.rstrip('/') or path.startswith(p.rstrip('/') + '/') for p in prefixes)
    changed = set(git('diff', '--name-only').stdout.splitlines()) | set(git('diff', '--cached', '--name-only').stdout.splitlines())
    forbidden = sorted(p for p in changed if not allowed(p))
    if forbidden:
        raise SystemExit('Preserve changes outside this task; do not publish them: ' + ', '.join(forbidden))
    save_task(name, task)
    git('add', '--', *prefixes)
    staged = git('diff', '--cached', '--name-only').stdout.splitlines()
    for name_in_git in staged:
        p = REPO / name_in_git
        if not allowed(name_in_git) or p.suffix == '.gpw' or (p.is_file() and p.stat().st_size >= 100_000_000):
            raise SystemExit('Unapproved path or oversized checkpoint was staged: ' + name_in_git)
    if not staged:
        print('No new task changes to publish.')
        return
    git('commit', '-m', 'Archive ' + (label or name) + ' and synchronize task status')
    push_or_stop(name, retry=True)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    sub.add_parser('identity')
    status = sub.add_parser('status'); status.add_argument('--remote', action='store_true')
    for action in ('claim', 'publish'):
        p = sub.add_parser(action); p.add_argument('task', choices=['window-a', 'window-b'])
        if action == 'publish': p.add_argument('--label')
    args = parser.parse_args()
    if args.command == 'identity': print(identity())
    elif args.command == 'status':
        if args.remote: git('fetch', 'origin', 'main')
        for name in ('window-a', 'window-b'):
            rel = str(task_path(name).relative_to(REPO))
            task = json.loads(git('show', 'origin/main:' + rel).stdout) if args.remote else load_task(name)
            print(json.dumps(task, ensure_ascii=False, indent=2))
    elif args.command == 'claim': claim(args.task)
    elif args.command == 'publish': publish(args.task, args.label)

if __name__ == '__main__': main()
