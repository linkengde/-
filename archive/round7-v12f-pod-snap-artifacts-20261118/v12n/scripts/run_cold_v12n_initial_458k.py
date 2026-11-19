#!/usr/bin/env python3
"""No-laser 1 fs x 10 step cold stability screen on the full initial-height model."""
from __future__ import annotations
import os, subprocess, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "Ag_Ti3SiC2_bicontinuous_200x200x300.data"
POD = Path("/workspace/potential-fit/round7_pod_snap/models/v12_candidate_forceorder_fixed/v12n_ccpair_candidate/POD")
LMP = Path("/workspace/software/lammps-pod-snap-clean/build/lmp")
OUT = HERE / "cold_v12n_initial_458k_1fs"
OUT.mkdir(exist_ok=True)
inp = OUT / "in.cold_0p01ps_1fs"
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
pair_coeff * * {POD/'AgTiSiC_POD_v12n.pod'} {POD/'pod_v12n_coefficients.pod'} Ag C Si Ti
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
group skeleton type 2 3 4
region fixed_layer block INF INF INF INF 0.0 20.0 units box
region bath_layer block INF INF INF INF 20.0 40.0 units box
group fixed_z region fixed_layer
group fixed intersect skeleton fixed_z
group bath region bath_layer
group mobile subtract all fixed
velocity mobile create 300.0 882733 mom yes rot no dist gaussian
velocity fixed set 0.0 0.0 0.0
fix freeze fixed setforce 0.0 0.0 0.0
compute Tmobile mobile temp
variable fmag atom sqrt(fx*fx+fy*fy+fz*fz)
compute fmax all reduce max v_fmag
timestep 0.001
thermo 1
thermo_style custom step time atoms c_Tmobile pe ke etotal c_fmax
thermo_modify temp Tmobile lost error flush yes format float %.16g
fix integrate mobile nve
fix heat_bath bath langevin 300.0 300.0 0.1 639827 tally yes
run 0
run 10
write_dump all custom {OUT/'end_0p01ps.dump'} id type x y z vx vy vz fx fy fz modify sort id
write_data {OUT/'end_0p01ps.data'}
""")
env = os.environ.copy()
env["LD_LIBRARY_PATH"] = "/workspace/software/openmpi-dev-5.0.7/usr/lib/x86_64-linux-gnu/openmpi/lib:" + str(LMP.parent) + ":" + env.get("LD_LIBRARY_PATH", "")
env["OMP_NUM_THREADS"] = "1"
cmd = ["/usr/bin/mpirun", "-np", "4", str(LMP), "-in", str(inp), "-log", "log.lammps"]
start = time.time()
result = subprocess.run(cmd, cwd=OUT, env=env, stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT, text=True, preexec_fn=lambda: os.nice(10))
(OUT / "stdout.log").write_text(result.stdout)
(OUT / "run_metadata.txt").write_text(
    f"POD v12n; corrected-height full model N=458248; no laser; dt=0.001 ps (1 fs); "
    f"duration=0.01 ps (10 steps); fixed skeleton z=0-20 A; z=20-40 A Langevin 300 K; "
    f"remaining mobile atoms NVE; MPI ranks=4; wall={time.time()-start:.3f}s; rc={result.returncode}\\n")
print(result.stdout[-6000:])
raise SystemExit(result.returncode)
