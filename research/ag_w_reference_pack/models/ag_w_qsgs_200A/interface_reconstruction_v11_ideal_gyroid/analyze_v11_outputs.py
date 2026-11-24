#!/usr/bin/env python3
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
L=200.;M_AG=107.8682;M_W=183.84
def sha(p):
 h=hashlib.sha256()
 with Path(p).open("rb") as f:
  for b in iter(lambda:f.read(1<<20),b""):h.update(b)
 return h.hexdigest()
def read_dump(path):
 lines=Path(path).read_text().splitlines();i=0
 while i<len(lines) and not lines[i].startswith("ITEM: ATOMS"):i+=1
 if i==len(lines):raise RuntimeError(f"No ATOMS section in {path}")
 cols=lines[i].split()[2:];arr=np.loadtxt(lines[i+1:],dtype=float)
 if arr.ndim==1:arr=arr[None,:]
 return {c:arr[:,j] for j,c in enumerate(cols)}
def cross_min(a,w):
 shifts=[[i*L,j*L,0.] for i in (-1,0,1) for j in (-1,0,1)]
 tw=cKDTree(np.concatenate([w+s for s in shifts]));ta=cKDTree(np.concatenate([a+s for s in shifts]))
 da,_=tw.query(a,k=1,workers=-1);dw,_=ta.query(w,k=1,workers=-1)
 return float(min(da.min(),dw.min())),float(np.quantile(da,.01)),float(np.quantile(dw,.01))
def parse_metrics(path):
 rows=[]
 for line in Path(path).read_text().splitlines():
  line=line.strip()
  if not line or line.startswith("#"):continue
  try:v=[float(x) for x in line.split()]
  except ValueError:continue
  if len(v)>=7:rows.append(v[:7])
 if not rows:raise RuntimeError("No thermo rows in cold metrics file")
 return np.asarray(rows)
def main():
 p=argparse.ArgumentParser()
 for x in ("out-dir","data","static-dump","cold-dump","cold-metrics","manifest"):p.add_argument("--"+x,type=Path,required=True)
 p.add_argument("--restart-step",type=int,default=8200);a=p.parse_args()
 out=a.out_dir;out.mkdir(parents=True,exist_ok=True)
 s=read_dump(a.static_dump);typ=s["type"].astype(int);z=s["z"];ids=s["id"].astype(np.int64)
 f=np.column_stack([s["fx"],s["fy"],s["fz"]]);fn=np.linalg.norm(f,axis=1);mobile=z>=20
 ag=typ==2;w=typ==1
 dmin,p01a,p01w=cross_min(np.column_stack([s["x"][ag],s["y"][ag],z[ag]]),np.column_stack([s["x"][w],s["y"][w],z[w]]))
 static={"atoms":int(len(ids)),"unique_ids":int(len(np.unique(ids))),"Ag_atoms":int(ag.sum()),"W_atoms":int(w.sum()),
  "mass_ratio_Ag_over_W":float(ag.sum()*M_AG/(w.sum()*M_W)),"force_max_mobile_eV_A":float(fn[mobile].max()),
  "force_rms_mobile_eV_A":float(np.sqrt(np.mean(fn[mobile]**2))),"force_p99_mobile_eV_A":float(np.quantile(fn[mobile],.99)),
  "min_cross_Ag_W_A":dmin,"cross_Ag_p01_A":p01a,"cross_W_p01_A":p01w,
  "atoms_outside_material_z_0_200_A":int(np.count_nonzero((z<0)|(z>200))),
  "static_exploratory_screen":"PASS" if fn[mobile].max()<5 and np.sqrt(np.mean(fn[mobile]**2))<.5 and dmin>2.3 else "FAIL"}
 (out/"v11_static_metrics.json").write_text(json.dumps(static,indent=2)+"\n")
 m=parse_metrics(a.cold_metrics);step=m[:,0];nvt=m[step<=8000];nve=m[step>8000];tail=nvt[nvt[:,0]>=7200]
 cold={"NVT_start_mobile_K":float(nvt[0,1]),"NVT_end_mobile_K":float(nvt[-1,1]),
  "NVT_end_Ag_K":float(nvt[-1,2]),"NVT_end_W_K":float(nvt[-1,3]),
  "NVT_end_Ag_W_delta_K":float(abs(nvt[-1,2]-nvt[-1,3])),
  "NVT_last_0p2ps_mobile_delta_K":float(tail[-1,1]-tail[0,1]),
  "NVT_last_0p2ps_PE_delta_eV":float(tail[-1,4]-tail[0,4]),
  "NVE_end_mobile_K":float(nve[-1,1]) if len(nve) else None,
  "NVE_end_Ag_K":float(nve[-1,2]) if len(nve) else None,
  "NVE_end_W_K":float(nve[-1,3]) if len(nve) else None,
  "NVE_end_Ag_W_delta_K":float(abs(nve[-1,2]-nve[-1,3])) if len(nve) else None,
  "NVE_energy_drift_eV":float(nve[-1,6]-nve[0,6]) if len(nve)>1 else None,
  "atom_count_preserved":int(len(read_dump(a.cold_dump)["id"]))==static["atoms"]}
 cold["cold_baseline_gate"]="PASS_EXPLORATORY" if (285<=cold["NVT_end_mobile_K"]<=315 and cold["NVT_end_Ag_W_delta_K"]<=15 and abs(cold["NVT_last_0p2ps_mobile_delta_K"])<=3 and cold["atom_count_preserved"]) else "FAIL_OR_CONTINUE_NVT"
 (out/"v11_cold_metrics.json").write_text(json.dumps(cold,indent=2)+"\n")
 manifest=json.loads(a.manifest.read_text())
 manifest["static_audit_status"]="STATIC_EFS_EXPLORATORY" if static["static_exploratory_screen"]=="PASS" else "STATIC_EFS_FAIL"
 manifest["md_start_gate"]="EXPLORATORY_ONLY";manifest["static_audit_metrics"]="v11_static_metrics.json"
 manifest["cold_baseline_metrics"]="v11_cold_metrics.json";manifest["restart_dump"]=a.cold_dump.name
 manifest["restart_dump_sha256"]=sha(a.cold_dump);manifest["restart_dump_timestep"]=a.restart_step
 a.manifest.write_text(json.dumps(manifest,indent=2)+"\n")
 print(json.dumps({"static":static,"cold":cold,"restart_sha256":manifest["restart_dump_sha256"]},indent=2))
 if static["static_exploratory_screen"]!="PASS":raise SystemExit("Static screen too severe for TTM smoke test")
if __name__=="__main__":main()
