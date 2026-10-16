#!/usr/bin/env python3
"""Verify GitHub-pinned frozen input structures. No ASE/GPAW required."""
import hashlib, json
from pathlib import Path
ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent
cfg=json.loads((ROOT/"INPUT_MANIFEST.json").read_text(encoding="utf-8"))
verified=[]
for item in cfg["records"]:
    file=SOURCE/item["relative_path"]
    data=file.read_bytes()
    sha=hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()
    if sha!=item["git_blob_sha1"]:
        raise SystemExit("INPUT_GIT_BLOB_SHA_MISMATCH: "+item["id"])
    lines=data.decode("utf-8").strip().splitlines()
    n=int(lines[0])
    if n!=item["atoms"] or len(lines)!=n+2: raise SystemExit("ATOM_COUNT_INVALID: "+item["id"])
    if 'pbc="T T T"' not in lines[1] or "Lattice=" not in lines[1]:
        raise SystemExit("PBC_CELL_MISSING: "+item["id"])
    if "REF_energy" in lines[1] or "REF_forces" in lines[1] or "PW_PBE_energy" in lines[1] or "PW_PBE_forces" in lines[1]:
        raise SystemExit("NONBLIND_LABEL_FOUND: "+item["id"])
    counts={}
    for ln in lines[2:]:
        fields=ln.split()
        if len(fields)<4:raise SystemExit("INVALID_ROW: "+item["id"])
        counts[fields[0]]=counts.get(fields[0],0)+1
        for x in fields[1:4]:
            v=float(x)
            if not (float("-inf")<v<float("inf")):raise SystemExit("NONFINITE_COORD")
    actual="".join(f"{s}{counts[s]}" for s in ("Ag","Ti","Si","C") if counts.get(s))
    if actual!=item["formula"]: raise SystemExit(f"FORMULA_MISMATCH: {item['id']} {actual} != {item['formula']}")
    verified.append(dict(id=item["id"],atoms=n,formula=actual,sha256=hashlib.sha256(data).hexdigest(),git_blob_sha1=sha))
print(json.dumps({"INPUT_CHECK":"PASS","count":len(verified),"DFT_STARTED":False,"verified":verified},indent=2,ensure_ascii=False))
