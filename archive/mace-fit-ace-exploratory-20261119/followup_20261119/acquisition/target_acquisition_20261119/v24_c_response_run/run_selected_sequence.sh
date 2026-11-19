#!/usr/bin/env bash
set -u
RUN_ROOT="/workspace/potential-comparison/mace_fit_data_20261118/target_acquisition_20261119/v24_c_response_run"
cd "$RUN_ROOT"
mkdir -p calculations
/workspace/.venvs/gpaw-diagnosis/bin/python checkpoint_keeper.py > calculations/checkpoint_keeper.log 2>&1 &
KEEPER_PID=$!
/workspace/.venvs/gpaw-diagnosis/bin/python recovery_monitor.py > calculations/recovery_monitor.log 2>&1 &
MONITOR_PID=$!
/workspace/.venvs/gpaw-diagnosis/bin/python run_sequence.py
SEQUENCE_RC=$?
wait "$MONITOR_PID"
MONITOR_RC=$?
kill "$KEEPER_PID" 2>/dev/null || true
wait "$KEEPER_PID" 2>/dev/null || true
/workspace/.venvs/gpaw-diagnosis/bin/python - <<'PY'
import json
from pathlib import Path
p=Path('sequence_status.json')
s=json.loads(p.read_text())
print(json.dumps({'sequence_status':s.get('status'),'current_label':s.get('current_label'),
                  'completed_labels':s.get('completed_labels',[]),'automatic_resume_attempts_by_label':s.get('automatic_resume_attempts_by_label',{}),
                  'sequence_return_code':int(__import__('os').environ.get('SEQUENCE_RC','0')),
                  'monitor_return_code':int(__import__('os').environ.get('MONITOR_RC','0'))},indent=2))
if s.get('status') != 'selected_three_archives_verified':
    raise SystemExit(1)
PY
