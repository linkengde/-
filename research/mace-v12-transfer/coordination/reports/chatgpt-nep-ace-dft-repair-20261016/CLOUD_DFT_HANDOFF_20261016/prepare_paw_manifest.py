#!/usr/bin/env python3
"""Create SHA-pinned GPAW PAW setup manifest; execute on the *DFT node*."""
import argparse, hashlib, json, sys
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument("--output",type=Path,required=True)
p.add_argument("--check-only",action="store_true")
a=p.parse_args()
import gpaw_data
folder=Path(gpaw_data.datapath()).resolve()
items={}
for s in ("Ag","Ti","Si","C"):
    f=folder/(s+".PBE.gz")
    if not f.is_file():sys.exit("MISSING_PAW_SETUP: "+str(f))
    items[s]={"path":str(f),"sha256":hashlib.sha256(f.read_bytes()).hexdigest()}
print(json.dumps({"actual_paw_root":str(folder),"files":items},indent=2))
if not a.check_only:
    if a.output.exists():sys.exit("REFUSING_TO_OVERWRITE_PAW_MANIFEST")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(items,indent=2)+"\n",encoding="utf-8")
