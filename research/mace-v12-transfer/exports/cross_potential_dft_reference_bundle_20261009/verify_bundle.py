"""Verify the immutable export inventory, without importing simulation software."""
from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parent
names=set()
for line in (root/'SHA256SUMS.txt').read_text().splitlines():
 digest,name=line.split('  ',1)
 p=(root/name).resolve()
 assert p.is_relative_to(root) and p.is_file() and name not in names
 assert hashlib.sha256(p.read_bytes()).hexdigest()==digest,name
 names.add(name)
actual={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and p.name!='SHA256SUMS.txt' and '__pycache__' not in p.parts}
# Nested archive inventories themselves are included in the root list.
actual.update(str(p.relative_to(root)) for p in root.rglob('SHA256SUMS.txt') if p.parent!=root)
assert actual==names, 'Missing or extra exported files'
m=json.loads((root/'bundle_manifest.json').read_text())
for r in m['entries']:
 assert hashlib.sha256((root/r['input']).read_bytes()).hexdigest()==r['input_sha256']
 assert json.loads((root/r['summary']).read_text())['scf_converged'] is True
print(f"PASS: {len(names)} files; {len(m['entries'])} verified reference cases")
