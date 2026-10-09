import json
from entry import ROOT,load
m=json.loads((ROOT/'input_manifest.json').read_text())
target='AgSi_rigid_Ag_z_parent_k6x6x1_sigma0p10_dipoleXY_v19'
assert m['launch_enabled'] is True and len(m['records'])==9
enabled=[r['label'] for r in m['records'] if r['launch_enabled']]
assert enabled==[target]
for r in m['records']:
 assert r['owner_task'] in ('window-a','window-b')
 assert r['launch_enabled'] is (r['label']==target)
 if r['label']==target:
  assert r['owner_task']=='window-a' and r['launch_status']=='enabled_A_parent_boundary_control_only'
 else:
  assert r['launch_status']=='blocked_pending_AgSi_parent_boundary_control_review'
 load(r['label'])
print('PASS: exactly one A-owned AgSi parent boundary control is enabled; other eight records remain disabled; hashes, transformations, geometry/source and method/ID guards pass')
