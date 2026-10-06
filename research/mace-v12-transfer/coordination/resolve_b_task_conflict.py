#!/usr/bin/env python3
"""Resolve only the B task-card rebase conflict, preserving remote assignments/local progress."""
import json
from pathlib import Path
import subprocess

TASK='research/mace-v12-transfer/coordination/tasks/window-b.json'
ASSIGNMENT_FIELDS={'current_stage','next_action','task_boundary','coordination_duties','result_directory','standing_work_instruction'}

def merge(base, remote, local, path=()):
    if remote==local:return remote
    if remote==base:return local
    if local==base:return remote
    if path==('live_progress',):
        return max((remote,local),key=lambda x:x.get('updated_utc',''))
    if len(path)==1 and path[0] in ASSIGNMENT_FIELDS:return remote
    if path==('updated_utc',):return max(remote,local)
    if isinstance(base,dict) and isinstance(remote,dict) and isinstance(local,dict):
        out={}
        for key in set(base)|set(remote)|set(local):
            if key not in base:
                if key in remote and key in local and remote[key]!=local[key]:
                    if isinstance(remote[key],dict) and isinstance(local[key],dict):out[key]=merge({},remote[key],local[key],path+(key,))
                    else:raise ValueError('Conflicting new key: '+'.'.join(path+(key,)))
                else:out[key]=remote.get(key,local.get(key))
            elif key not in remote or key not in local:
                raise ValueError('Deletion conflict: '+'.'.join(path+(key,)))
            else:out[key]=merge(base[key],remote[key],local[key],path+(key,))
        return out
    if len(path)==1 and path[0] in {'queued_work_order','allowed_publish_paths'}:
        return remote+[value for value in local if value not in remote]
    raise ValueError('Unresolved simultaneous change: '+'.'.join(path))

def main():
    repo=Path(subprocess.check_output(['git','rev-parse','--show-toplevel'],text=True).strip())
    unresolved=subprocess.check_output(['git','diff','--name-only','--diff-filter=U'],cwd=repo,text=True).splitlines()
    if unresolved!=[TASK]:raise SystemExit('Expected exactly one unresolved B task card; preserve current state and inspect')
    stages=[json.loads(subprocess.check_output(['git','show',f':{stage}:{TASK}'],cwd=repo)) for stage in (1,2,3)]
    base,remote,local=stages
    if remote['task_id']!='window-b' or local['task_id']!='window-b' or remote['owner_instance']!=local['owner_instance'] or remote['owner_instance']!='c5035b48-f43b-4da0-b8e4-e2862817f86a':raise SystemExit('B owner/task conflict; do not merge automatically')
    result=merge(base,remote,local)
    for snapshot in (remote,local):
        for job in snapshot.get('current_jobs',{}):
            if job not in result['current_jobs']:raise SystemExit('Lost assigned job')
        for label in snapshot.get('labels',{}):
            if label not in result['labels']:raise SystemExit('Lost DFT label')
    result['task_card_merge_note']='Preserved concurrent A assignment fields and latest B live progress; no owner change. Other scalar conflicts require manual review.'
    backup=Path('/tmp/window-b-task-conflict-before.json')
    if backup.exists():raise SystemExit('Backup already exists; use prior snapshot and inspect before retry')
    backup.write_text(json.dumps({'base':base,'remote':remote,'local':local},indent=2)+'\n')
    (repo/TASK).write_text(json.dumps(result,indent=2)+'\n')
    subprocess.run(['git','add',TASK],cwd=repo,check=True)
    print('B task card resolved and staged; original three-way snapshot in',backup)
    print('Review staged card, then GIT_EDITOR=true git rebase --continue. No push performed.')

if __name__=='__main__':main()
