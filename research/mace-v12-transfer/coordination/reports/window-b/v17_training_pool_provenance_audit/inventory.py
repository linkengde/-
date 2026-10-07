"""Inventory only explicitly permitted training evidence; no blind outputs."""
import csv, hashlib, json
from pathlib import Path
from ase.io import read
ROOT=Path(__file__).resolve().parents[6]
OUT=Path(__file__).resolve().parent
BASE=ROOT/'research/mace-v12-transfer'
DATA=BASE/'periodic_interface_v4/mace_periodic_v16_interface_energy/data'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def run():
 m=json.loads((DATA/'dataset_manifest.json').read_text())
 frames=read(DATA/'train.extxyz',':'); assert len(frames)==35
 assert sha(DATA/'train.extxyz')==m['output_sha256']['train']
 allowed={DATA/'train.extxyz':m['output_sha256']['train'],DATA/'dataset_manifest.json':None}
 reports=BASE/'coordination/reports/window-b'
 for name,file in [('v15_inherited_method_provenance_trace','per_frame_ledger.json'),('v16_dataset_integration_audit','audit.json'),('v17_energy_convention_and_parent_hash_audit','audit.json')]:
  allowed[reports/name/file]=m['sources_sha256'].get(str((reports/name/file).relative_to(ROOT)))
 for row in m['per_frame_provenance']:
  if row['split']=='train': allowed[ROOT/row['source_path']]=row['source_sha256']
 for path,digest in m['V16_sources_sha256'].items():
  if '/calculations/' in path and any(x['source']+'/' in path for x in m['V16_additions']): allowed[ROOT/path]=digest
 inventory=[]
 for path,pin in sorted(allowed.items()):
  assert path.is_file(),str(path)
  actual=sha(path); assert pin is None or pin==actual,(str(path),pin,actual)
  inventory.append({'path':str(path.relative_to(ROOT)),'sha256':actual,'expected_sha256':pin,'check':'PASS' if pin else 'recorded_current_bytes'})
 ledger=[{'frame':i,'config_type':a.info.get('config_type'),'atoms':len(a),'formula':a.get_chemical_formula(),'evidence_status':'pending_full_trace','training_file_sha256':sha(DATA/'train.extxyz')} for i,a in enumerate(frames)]
 (OUT/'source_inventory.json').write_text(json.dumps({'scope':'training evidence only','files':inventory},indent=2)+'\n')
 (OUT/'ledger_skeleton.json').write_text(json.dumps(ledger,indent=2)+'\n')
 with (OUT/'ledger_skeleton.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(ledger[0]));w.writeheader();w.writerows(ledger)
 paths=['inventory.py','source_inventory.json','ledger_skeleton.json','ledger_skeleton.csv']
 (OUT/'SHA256SUMS.txt').write_text(''.join(f'{sha(OUT/p)}  {p}\n' for p in paths))
 for line in (OUT/'SHA256SUMS.txt').read_text().splitlines():
  digest,name=line.split('  ');assert sha(OUT/name)==digest
 print(f'PASS: {len(inventory)} source files, {len(ledger)} training rows; artifact hashes verified')
if __name__=='__main__':run()
