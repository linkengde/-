#!/usr/bin/env python3
"""Full 458k initial-geometry run-0 force screen using POD v12k."""
from __future__ import annotations
import os, subprocess, time, json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE / "Ag_Ti3SiC2_bicontinuous_200x200x300.data"
POD = Path("/workspace/potential-fit/round7_pod_snap/models/v12_candidate_forceorder_fixed/v12k_ticpair_candidate/POD")
LMP = Path("/workspace/software/lammps-pod-snap-clean/build/lmp")
OUT = HERE / "static_v12k_initial_458k"
OUT.mkdir(exist_ok=True)
inp = OUT / "in.run0"
dump = OUT / "forces.dump"
inp.write_text(f"""clear
units metal
atom_style atomic
boundary p p m
processors 4 1 1
atom_modify sort 0 0
read_data {DATA}
# Source types Ti/Si/C/Ag -> model types Ag/C/Si/Ti.
group src_Ti type 1
group src_Si type 2
group src_C type 3
group src_Ag type 4
set group src_Ti type 4
set group src_Si type 3
set group src_C type 2
set group src_Ag type 1
mass 1 107.8682
mass 2 12.011
mass 3 28.085
mass 4 47.867
pair_style pod
pair_coeff * * {POD/'AgTiSiC_POD_v12k.pod'} {POD/'pod_v12k_coefficients.pod'} Ag C Si Ti
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
group skeleton type 2 3 4
region fixed_layer block INF INF INF INF 0.0 20.0 units box
group fixed_z region fixed_layer
group fixed intersect skeleton fixed_z
fix freeze fixed setforce 0.0 0.0 0.0
variable fmag atom sqrt(fx*fx+fy*fy+fz*fz)
compute fmax all reduce max v_fmag
thermo 1
thermo_style custom step atoms pe c_fmax
thermo_modify lost error flush yes format float %.16g
run 0
write_dump all custom {dump} id type x y z fx fy fz modify sort id
""")
env=os.environ.copy()
env["LD_LIBRARY_PATH"]="/workspace/software/openmpi-dev-5.0.7/usr/lib/x86_64-linux-gnu/openmpi/lib:"+str(LMP.parent)+":"+env.get("LD_LIBRARY_PATH","")
env["OMP_NUM_THREADS"]="1"
start=time.time()
res=subprocess.run(["/usr/bin/mpirun","-np","4",str(LMP),"-in",str(inp),"-log","log.lammps"],
                   cwd=OUT,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,
                   preexec_fn=lambda:os.nice(10))
(OUT/"stdout.log").write_text(res.stdout)
(OUT/"run_metadata.txt").write_text(f"POD v12k; N=458248; run 0; fixed Ti3SiC2 skeleton z=0-20 A; ranks=4; wall_s={time.time()-start:.3f}; rc={res.returncode}\\n")
print(res.stdout[-3000:])
if res.returncode: raise SystemExit(res.returncode)
lines=dump.read_text().splitlines(); h=next(i for i,s in enumerate(lines) if s.startswith("ITEM: ATOMS ")); c=lines[h].split()[2:]
cols={k:c.index(k) for k in ("id","type","fx","fy","fz")}
a=np.array([[float(v) for v in lines[j].split()] for j in range(h+1,h+1+458248)])
f=a[:,[cols[x] for x in ("fx","fy","fz")]]; typ=a[:,cols["type"]].astype(int); mag=np.linalg.norm(f,axis=1)
names={1:"Ag",2:"C",3:"Si",4:"Ti"}
summary={"atoms":int(len(a)),"unique_ids":int(len(np.unique(a[:,cols["id"]]))),"nonfinite_force_rows":int((~np.isfinite(f)).any(axis=1).sum()),
 "Fmax_eV_A":float(mag.max()),"F99_eV_A":float(np.quantile(mag,.99)),
 "Fmax_by_element_eV_A":{names[int(t)]:float(mag[typ==t].max()) for t in np.unique(typ)},
 "atoms_force_above_500":int((mag>500).sum()),"role":"static screen only; not MD or heating approval"}
(OUT/"summary.json").write_text(json.dumps(summary,indent=2)+"\\n")
print(json.dumps(summary,indent=2))
