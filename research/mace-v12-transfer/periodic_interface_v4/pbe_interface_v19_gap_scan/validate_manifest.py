import json
from entry import ROOT,load
m=json.loads((ROOT/'input_manifest.json').read_text())
assert m['launch_enabled'] is False and len(m['records'])==9
for r in m['records']:
 assert r['launch_enabled'] is False and r['owner_task'] in ('window-a','window-b')
 load(r['label'])
print('PASS: nine disabled proposals, transformations, geometry/source hashes and method/ID guards')
