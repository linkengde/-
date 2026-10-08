"""User-authorized post-completion owner transfer; no computation or forced Git writes."""
import fcntl,hashlib,json,subprocess,sys
from pathlib import Path
import sync_tasks as s
if len(sys.argv)!=2 or sys.argv[1] not in ('window-a','window-b'):raise SystemExit('usage: claim_restarted_window.py window-a|window-b')
name=sys.argv[1]
with s.LOCK.open('a') as lock:
 fcntl.flock(lock,fcntl.LOCK_EX)
 if s.git('diff','--name-only').stdout.strip() or s.git('diff','--cached','--name-only').stdout.strip():raise SystemExit('Preserve tracked/index changes before transfer')
 s.git('fetch','origin','main');s.git('merge','--ff-only','origin/main')
 h=json.loads((s.HERE/'RESTART_V19_HANDOVER.json').read_text());assert h['ready'] is True
 for r in h['critical_files']:
  p=s.REPO/r['path'];assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],'Handover evidence mismatch: '+r['path']
 current=s.identity();old=h['previous_owners'][name];task=s.load_task(name)
 if current in h['previous_owners'].values():raise SystemExit('This is an old registered instance. Use a genuinely new independent cloud environment; no owner change performed')
 if task['owner_instance']==current:print('Already transferred to this new instance');raise SystemExit(0)
 if task['owner_instance']!=old:raise SystemExit('Task already transferred to another instance; stop, no duplicate work')
 other=s.load_task('window-b' if name=='window-a' else 'window-a');assert other['owner_instance']!=current,'A/B require independent instances'
 for p in Path('/proc').iterdir():
  if not p.name.isdigit():continue
  try:args=(p/'cmdline').read_bytes().split(b'\0')
  except OSError:continue
  if any(Path(x.decode(errors='ignore')).name=='run_pw_reference.py' for x in args):raise SystemExit('Active GPAW runner detected; no transfer')
 progress=s.PROGRESS/(name+'.json');archive=s.PROGRESS/'previous_instances'/h['handover_id']/(name+'.json');archive.parent.mkdir(parents=True,exist_ok=True);assert not archive.exists(),'Old progress archive already exists'
 archive.write_bytes(progress.read_bytes())
 task.setdefault('owner_history',[]).append({'owner_instance':old,'handover_id':h['handover_id'],'transferred_utc':s.now(),'reason':'User requested two fresh environments after both owners finished/published mesh pilots'})
 task['owner_instance']=current;task['claimed_utc']=s.now();task['status']='claimed';task['current_stage']='New cloud instance claimed after completed/published four-mesh handover; no job launched';task['live_progress']=None;task['handover_consumed']=h['handover_id'];task['next_action']='Read coordination/RESTART_V19_NEXT_STEPS.md. Verify environment first. New A prepares owner-bound5x5 and6x6 entries after both new IDs registered; A runs5x5, B runs6x6. No Gamma/2x2/3x3/4x4 repetition; training remains blocked.'
 s.save_task(name,task);progress.write_text(json.dumps({'task_id':name,'owner_instance':current,'events':[],'updated_utc':s.now()},indent=2)+'\n')
 paths=[s.task_path(name),progress,archive];s.git('add','--',*[str(p.relative_to(s.REPO)) for p in paths]);s.git('commit','-m',f'Transfer {name} to fresh cloud instance after completed mesh handover')
 pushed=s.git('push','origin','HEAD:main',check=False)
 if pushed.returncode:
  s.git('fetch','origin','main');remote=json.loads(s.git('show','origin/main:'+str(s.task_path(name).relative_to(s.REPO))).stdout)
  if remote['owner_instance']!=old:raise SystemExit('Remote owner changed; preserve local commit, do not compute or force push')
  s.git('rebase','origin/main');s.git('push','origin','HEAD:main')
 print('Confirmed new registered owner:',name,current,'; no computation started')
