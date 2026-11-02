#!/usr/bin/env python3
"""Resume verified/pending v25 skeleton labels, reusing only matching checkpoints."""
from pathlib import Path
import json, subprocess, sys
ROOT=Path(__file__).resolve().parent
REPO=Path('/workspace/jeitschko-main-run')
PYTHON='/workspace/.venvs/gpaw-mpi/bin/python'
for p in Path('/proc').iterdir():
    if not p.name.isdigit() or int(p.name)==__import__('os').getpid(): continue
    try:
        cmd=(p/'cmdline').read_bytes().decode(errors='ignore')
        state=(p/'stat').read_text().rsplit(')',1)[1].strip().split()[0]
    except (OSError,IndexError): continue
    if state!='Z' and any(x in cmd for x in ('run_sequence.py','run_preflight_label.py')):
        raise SystemExit(f'refusing resume while active DFT process exists: PID {p.name}')
status=json.loads((ROOT/'sequence_status.json').read_text()) if (ROOT/'sequence_status.json').exists() else {}
if status.get('status') not in ('paused_by_user','failed_stop_no_retry','running','inputs_prepared'):
    raise SystemExit(f"unexpected resume status: {status.get('status')}")
p=subprocess.run([PYTHON,str(ROOT/'run_sequence.py'),'--resume'],cwd=REPO)
raise SystemExit(p.returncode)
