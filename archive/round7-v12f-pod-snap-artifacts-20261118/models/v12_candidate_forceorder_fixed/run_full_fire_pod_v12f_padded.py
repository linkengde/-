#!/usr/bin/env python3
"""Launch an isolated, boundary-padded FIRE pre-relaxation for the full model."""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = Path(
    "/workspace/potential-fit/round7_pod_snap/models/target_200x200x300_qsgs/"
    "md_runs/attempt_008_POD_TTM_mpi_big/start_wrapped.data"
)
LMP = Path("/workspace/software/lammps-pod-snap-clean/build/lmp")
OUT = HERE / "full_fire_347k_padded" / "POD"
OUT.mkdir(parents=True, exist_ok=True)
IN = OUT / "in.fire_100_padded"

IN.write_text(
    f"""clear
units metal
atom_style atomic
boundary p p f
processors 4 1 1
atom_modify sort 0 0
read_data {DATA}
# Keep the physical layer coordinates and add numerical headroom below z=0.
change_box all z final -10.0 300.0 units box
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
pair_coeff * * {HERE / 'POD/AgTiSiC_POD_v12f.pod'} {HERE / 'POD/pod_v12f_coefficients.pod'} Ag C Si Ti
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
group skeleton type 2 3 4
region fixed_layer block INF INF INF INF 0.0 20.0 units box
group fixed_z region fixed_layer
group fixed_skeleton intersect fixed_z skeleton
fix freeze fixed_skeleton setforce 0.0 0.0 0.0
variable fmag atom sqrt(fx*fx+fy*fy+fz*fz)
compute fmax all reduce max v_fmag
thermo 10
thermo_style custom step atoms pe c_fmax
thermo_modify lost error flush yes format float %.16g
write_dump all custom {OUT / 'initial.dump'} id type x y z modify sort id
min_style fire
min_modify dmax 0.02
minimize 0.0 0.5 100 5000
write_dump all custom {OUT / 'forces_minimized.dump'} id type x y z fx fy fz modify sort id
write_data {OUT / 'minimized.data'}
"""
)

env = os.environ.copy()
env["LD_LIBRARY_PATH"] = (
    "/workspace/software/openmpi-dev-5.0.7/usr/lib/x86_64-linux-gnu/openmpi/lib:"
    + str(LMP.parent)
    + ":"
    + env.get("LD_LIBRARY_PATH", "")
)
env["OMP_NUM_THREADS"] = "1"
stdout_path = OUT / "stdout.log"
with stdout_path.open("w") as stdout:
    proc = subprocess.Popen(
        ["/usr/bin/mpirun", "-np", "4", str(LMP), "-in", str(IN), "-log", "log.lammps"],
        cwd=OUT,
        env=env,
        stdout=stdout,
        stderr=subprocess.STDOUT,
        start_new_session=True,
    )
(OUT / "run_metadata.txt").write_text(
    "POD v12f; full model N=347464; FIRE max=100 iterations, ftol=0.5 eV/A, "
    "dmax=0.02 A; computational z bounds=-10..300 A; fixed skeleton z=0..20 A; "
    f"MPI ranks=4; launcher_pid={proc.pid}; source unchanged={DATA}\n"
)
print(f"launched pid={proc.pid}")
print(f"input={IN}")
print(f"output={OUT}")
