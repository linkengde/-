#!/usr/bin/env python3
"""Read-only owner preflight; no installation, calculator, MD or DFT."""
import argparse, hashlib, importlib.metadata, json, os
from pathlib import Path
import socket, subprocess, sys

SETUP_SHA = {
 "Ag": "bff00ac79a7fa6e8c3a71c8b34d3a9c488e6bf1c00b3e9c4c91c231551145e97",
 "Ti": "bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f",
 "Si": "2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0",
 "C": "fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983"
}
def sha(p):
 h=hashlib.sha256()
 with Path(p).open("rb") as f:
  for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
 return h.hexdigest()

def machine_fingerprint(work_root):
 # Bind this owner report to the local machine and mounted work root.
 # This is identity evidence, not proof of physical resource independence.
 root=Path(work_root).resolve()
 if not root.is_dir(): return None
 try:
  identity={"hostname":socket.gethostname(),
   "boot_id":Path("/proc/sys/kernel/random/boot_id").read_text().strip(),
   "cgroup":Path("/proc/self/cgroup").read_text().strip(),
   "work_root":str(root),"work_root_device":os.stat(root).st_dev}
 except OSError: return None
 return hashlib.sha256(json.dumps(identity,sort_keys=True,separators=(",",":")).encode()).hexdigest()

def cg_values(name):
 # Include all visible ancestors; provider limits beyond the mount stay unknown.
 root=Path("/sys/fs/cgroup")
 relative=""
 try:
  for line in Path("/proc/self/cgroup").read_text().splitlines():
   if line.startswith("0::"): relative=line[3:].lstrip("/")
 except OSError: pass
 cur=root/relative
 if not cur.exists(): cur=root
 ans=[]
 while cur==root or root in cur.parents:
  p=cur/name
  try: ans.append(p.read_text().strip())
  except OSError: pass
  if cur==root: break
  cur=cur.parent
 return ans

def resources(work_root):
 affinity=len(os.sched_getaffinity(0)) if hasattr(os,"sched_getaffinity") else os.cpu_count()
 quotas=[]
 for s in cg_values("cpu.max"):
  a=s.split()
  if len(a)==2 and a[0]!="max": quotas.append(int(a[0])/int(a[1]))
 effective=min([affinity]+quotas) if quotas else None
 meminfo={}
 for line in Path("/proc/meminfo").read_text().splitlines():
  k,v=line.split(":",1)
  meminfo[k]=int(v.split()[0])*1024
 hostavail=meminfo.get("MemAvailable")
 lims=cg_values("memory.max"); used=cg_values("memory.current")
 heads=[]
 for lim,u in zip(lims,used):
  if lim!="max": heads.append(max(0,int(lim)-int(u)))
 avail=min(([hostavail] if hostavail is not None else [])+heads) if heads else None
 root=Path(work_root)
 path=root
 while not path.exists() and path!=path.parent: path=path.parent
 stat=os.statvfs(path)
 return {"cpu_visible":os.cpu_count(),"cpu_affinity":affinity,"cpu_max_visible_ancestors":cg_values("cpu.max"),
   "cpu_quota_effective":effective,"host_MemAvailable_bytes":hostavail,
   "memory_limit_visible":lims,"memory_current_visible":used,"memory_available_bytes":avail,
   "work_root_exists":root.exists(),"measurement_path":str(path),"disk_available_bytes":stat.f_bavail*stat.f_frsize,
   "disk_total_bytes":stat.f_blocks*stat.f_frsize,"inodes_available":stat.f_favail,
   "provider_limits":"UNVERIFIED","volume_quota_bytes":"UNVERIFIED"}

def main():
 a=argparse.ArgumentParser(description=__doc__)
 a.add_argument("--owner",choices=["C","D","E","F"],required=True)
 a.add_argument("--work-root",required=True)
 a.add_argument("--runtime",required=True)
 v=a.parse_args()
 root=Path(v.work_root)
 if not root.is_absolute(): a.error("work root must be absolute")
 runtime=Path(v.runtime)
 report={"schema_version":1,"owner":v.owner,"identity_grade":"OWNER_DECLARED; NOT_C_REMOTE_VERIFIED",
    "work_root":str(root.resolve()),"persistence_verified":False,"independent_resources_verified":False,
    "machine_fingerprint_sha256":machine_fingerprint(root),
    "persistence_evidence":"UNVERIFIED; attach provider/mount evidence for C review",
    "runtime_path":str(runtime),"runtime_sha256":sha(runtime) if runtime.is_file() else None,
    "runtime_syntax_pass":False,"resources":resources(root),"versions":{},
    "GPAW_MPI_BACKEND":os.environ.get("GPAW_MPI_BACKEND"),
    "GPAW_SETUP_PATH":os.environ.get("GPAW_SETUP_PATH"),"setup_hashes":{},
    "scientific_jobs_started":0}
 if runtime.is_file():
  check=subprocess.run(["bash","-n",str(runtime)],capture_output=True,text=True)
  report["runtime_syntax_pass"]=check.returncode==0
  report["runtime_syntax_error"]=check.stderr[:1000]
 for p in ["gpaw","ase","gpaw-data","numpy","scipy"]:
  try: report["versions"][p]=importlib.metadata.version(p)
  except importlib.metadata.PackageNotFoundError: report["versions"][p]="MISSING"
 paths=[Path(x) for x in os.environ.get("GPAW_SETUP_PATH","").split(os.pathsep) if x]
 for element,expected in SETUP_SHA.items():
  files=[p/(element+".PBE.gz") for p in paths if (p/(element+".PBE.gz")).is_file()]
  report["setup_hashes"][element]={"lookup_matches":len(files),"sha256":sha(files[0]) if files else None,
    "expected":expected,"matches_project":bool(files) and sha(files[0])==expected}
 report["scientific_release"]="HOLD; this helper never signs C release"
 print(json.dumps(report,ensure_ascii=False,indent=2))
if __name__=="__main__": main()
