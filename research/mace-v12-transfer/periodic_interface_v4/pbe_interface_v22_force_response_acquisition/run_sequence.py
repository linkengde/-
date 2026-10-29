#!/usr/bin/env python3
"""Run the six fixed-geometry labels sequentially, verifying each before continuing."""
from __future__ import annotations
import json,subprocess,sys
from datetime import datetime,timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parent
PYTHON='/workspace/.venvs/gpaw-mpi/bin/python'
LAUNCH=ROOT/'launch_preflight.py'
VERIFY=ROOT/'verify_preflight.py'
manifest=json.loads((ROOT/'input_manifest.json').read_text())
jobs=manifest['preflight_jobs']
priority={'AgTi':0,'AgSi':1,'AgC':2}
jobs=sorted(jobs,key=lambda x:(priority[x['label'].split('_',1)[0]], 0 if '_m0p20_' in x['label'] else 1))
status_path=ROOT/'sequence_status.json'
status={'status':'running','started_utc':datetime.now(timezone.utc).isoformat(),'completed_labels':[],'current_label':None,'automatic_retries':0}
def save(): status_path.write_text(json.dumps(status,indent=2)+'\n')
for job in jobs:
    label=job['label']; status['current_label']=label; save()
    print(f'=== {label}: resource preflight and DFT ===',flush=True)
    proc=subprocess.Popen([PYTHON,str(LAUNCH),label],cwd=ROOT.parents[3],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1)
    assert proc.stdout is not None
    for line in proc.stdout: print(line,end='',flush=True)
    code=proc.wait()
    if code!=0:
        status.update({'status':'failed_stop_no_retry','failed_label':label,'launcher_exit_code':code,'finished_utc':datetime.now(timezone.utc).isoformat()}); save()
        raise SystemExit(code)
    verified=subprocess.run([PYTHON,str(VERIFY),label],cwd=ROOT.parents[3],capture_output=True,text=True)
    if verified.returncode!=0:
        status.update({'status':'archive_verification_failed_stop_no_retry','failed_label':label,'verification_stderr':verified.stderr,'finished_utc':datetime.now(timezone.utc).isoformat()}); save()
        print(verified.stdout,flush=True); print(verified.stderr,file=sys.stderr,flush=True)
        raise SystemExit(verified.returncode)
    result=json.loads(verified.stdout)
    if result.get('status')!='PASS':
        status.update({'status':'archive_verification_failed_stop_no_retry','failed_label':label,'verification':result,'finished_utc':datetime.now(timezone.utc).isoformat()}); save()
        raise SystemExit('archive verifier did not pass')
    status['completed_labels'].append({'label':label,'archive_verification':'PASS','archive_sha256':result['archive_sha256']})
    status['current_label']=None; save()
status.update({'status':'all_six_archives_verified','finished_utc':datetime.now(timezone.utc).isoformat()}); save()
print(json.dumps(status,indent=2),flush=True)
