import json,csv,hashlib
from pathlib import Path
import numpy as np
from ase.io import read
ROOT=Path(__file__).resolve().parents[6];OUT=Path(__file__).resolve().parent;BASE=ROOT/'research/mace-v12-transfer/periodic_interface_v4';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
train=BASE/'mace_periodic_v19_residual_response/training_candidate/train.extxyz';a=read(train,':');manifest=json.loads((BASE/'mace_periodic_v19_residual_response/training_candidate/manifest.json').read_text());assert sha(train)==manifest['train_sha256'] if 'train_sha256' in manifest else True
entries=[]
for p in sorted((ROOT/'research/mace-v12-transfer').glob('**/*.extxyz')):
 rel=str(p.relative_to(ROOT))
 if any(x in rel.lower() for x in ['withheld','blind','holdout','dev_validation','sealed']):continue
 if 'v19_residual_direction_geometry_design/geometries' in rel or p==train:
  entries.append({'path':rel,'sha256':sha(p)})
# No candidate can be responsibly proposed: every authorized residual mode uses one
# of the three V17 scored development parent geometries, and new labels are local
# changes within those exact base motifs. Geometry-generation stops at this evidence.
parents={'AgC':'AgC_v17_dev_validation_01','AgSi':'AgSi_v17_dev_validation_01','AgTi':'AgTi_v17_dev_validation_01'}
roles=json.loads((ROOT/'research/mace-v12-transfer/coordination/reports/window-b/v17_split_geometry_design/final_set/role_manifest.json').read_text());
for r in roles['candidates']:
 if r['proposed_role']=='development_validation': assert r['candidate_id'] in parents.values()
rows=[]
for iface,pid in parents.items():
 labels=[x.get('config_type',x.get('config','')) for x in manifest['base_rows']+manifest['new_rows'] if iface in x.get('config_type',x.get('config',''))]
 candidates=[r for r in json.loads((ROOT/'research/mace-v12-transfer/coordination/reports/window-b/v19_residual_direction_geometry_design/manifest.json').read_text())['records'] if r['interface']==iface]
 rows.append({'interface':iface,'excluded_development_parent':pid,'training36_rows_same_interface':len(labels),'distinct_nonsealed_source_parent_identified':False,'v19_candidate_count':len(candidates),'independent_geometry_proposals':0,'blocker':'all mapped candidates are residual/finite-displacement descendants of reused development parent family; no public distinct source family with compatible geometry/lineage identified'})
result={'status':'BLOCKED_NO_INDEPENDENT_PARENT','recommendation':'Do not label a derivative of the reused V17 development parents as fresh development validation. Need a new public structural parent family, with owner-approved DFT labels after A role review.','per_interface':rows,'candidate_count':0,'source_inventory':entries,'sealed_geometry_contents_opened':False,'new_calculation':False}
(OUT/'manifest.json').write_text(json.dumps(result,indent=2)+'\n');(OUT/'report.md').write_text('# V19 independent development geometry design: blocked\n\nNo independent development candidate is proposed. The six V19 DFT labels and 18 proposals are explicitly generated from the same three structures already scored for V17/V18. They add local derivative coverage but do not create new geometry families. Selecting another small displacement of those parents would keep the development/ acquisition lineage coupled. The integrated V19 training36 candidate confirms exactly three interface composition families and contains the six descendants. Public lineage and hashes are in manifest. Sealed candidate geometry contents were not opened.\n\nMinimum remedy: A should source or construct one new interface parent family per interface with distinct registry/framework-neighbor topology, freeze roles before labeling, and have its assigned DFT owners label fixed structures under recorded settings. Then compare exact/periodic species-aware near geometry against the complete training36 and the six V19 structures; if a distinct topology cannot be made inside these small cells, report the limitation and expand the parent cell/source. Do not recycle a training derivative as an independent test. This design does not authorize DFT or role changes.\n');(OUT/'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='SHA256SUMS.txt'));print(result['status'],len(entries))
