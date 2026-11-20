#!/workspace/venvs/dft-labels/bin/python
"""Resume the v12p checks after POD static; never launches heating."""
from __future__ import annotations
import json, subprocess, time, traceback
from pathlib import Path

ROOT=Path(__file__).resolve().parent
BASE=Path('/workspace/potential-fit/round7_pod_snap/models/v12_candidate_forceorder_fixed')
PY='/workspace/venvs/dft-labels/bin/python'
LOG=ROOT/'v12p_remaining_pipeline.log'
POD_COLD_PID=348231

def log(s):
    with LOG.open('a') as f:
        f.write(f'[{time.strftime("%Y-%m-%d %H:%M:%S UTC",time.gmtime())}] {s}\n'); f.flush()
def alive(pid):
    p=Path(f'/proc/{pid}/stat')
    if not p.exists(): return False
    try: state=p.read_text().rsplit(')',1)[1].strip().split()[0]
    except Exception: return False
    return state!='Z'
def run(pyfile,logfile):
    cmd=[PY,str(ROOT/pyfile)]
    log(f'RUN {cmd}')
    with logfile.open('w') as f: p=subprocess.run(cmd,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,text=True)
    log(f'EXIT {p.returncode} log={logfile}')
    return p.returncode

def audit(name,logfile):
    rc=run(name,logfile)
    return rc

LOG.write_text('')
log(f'Wait for active POD cold run PID {POD_COLD_PID}')
while alive(POD_COLD_PID): time.sleep(10)

# Audit POD cold if LAMMPS produced a complete dump.
pod_cold=ROOT/'cold_v12p_initial_458k_1fs'
if (pod_cold/'end_0p01ps.dump').exists():
    rc=audit('audit_cold_v12p_initial_458k.py',ROOT/'pod_v12p_cold_audit.log')
    if rc==0:
        try: pod_result=json.loads((pod_cold/'cold_screen_summary.json').read_text())
        except Exception: pod_result={'status':'AUDIT_FAILED'}
    else: pod_result={'status':'AUDIT_FAILED','rc':rc}
else:
    pod_result={'status':'COLD_RUN_FAILED_OR_INCOMPLETE','stdout':str(pod_cold/'stdout.log')}
log('POD_COLD '+json.dumps(pod_result,sort_keys=True))

# Run SNAP full-model static screen, and only run cold MD if static forces pass the emergency gate.
snap_result={'status':'STATIC_NOT_RUN'}
rc=run('run_static_snap_v12p_458k.py',ROOT/'snap_v12p_static.log')
static_path=ROOT/'static_snap_v12p_initial_458k'/'summary.json'
if rc==0 and static_path.exists():
    static=json.loads(static_path.read_text()); log('SNAP_STATIC '+json.dumps(static,sort_keys=True))
    if static.get('atoms')==458248 and static.get('nonfinite_force_rows')==0 and static.get('atoms_force_above_500',1)==0:
        rc=run('run_cold_snap_v12p_initial_458k.py',ROOT/'snap_v12p_cold.log')
        cold_dir=ROOT/'cold_snap_v12p_initial_458k_1fs'
        if rc==0 and (cold_dir/'end_0p01ps.dump').exists():
            arc=run('audit_cold_snap_v12p_initial_458k.py',ROOT/'snap_v12p_cold_audit.log')
            if arc==0 and (cold_dir/'cold_screen_summary.json').exists(): snap_result=json.loads((cold_dir/'cold_screen_summary.json').read_text())
            else: snap_result={'status':'COLD_AUDIT_FAILED','rc':arc}
        else: snap_result={'status':'COLD_RUN_FAILED_OR_INCOMPLETE','rc':rc}
    else: snap_result={'status':'STATIC_GATE_FAILED','static':static}
else: snap_result={'status':'STATIC_RUN_FAILED','rc':rc}
log('SNAP_COLD '+json.dumps(snap_result,sort_keys=True))
status={'status':'V12P_CHECKS_COMPLETED','heating_started':False,'POD_static':json.loads((ROOT/'static_v12p_initial_458k'/'summary.json').read_text()),'POD_holdout':json.loads((BASE/'v12p_10fs_ag_contacts_candidate'/'holdout_si01'/'summary.json').read_text()),'POD_cold':pod_result,'SNAP_holdout_and_cold':snap_result}
(ROOT/'v12p_validation_status.json').write_text(json.dumps(status,indent=2)+'\n')
log('DONE')
print(json.dumps(status,indent=2))
