# Basilisk build and reproduction

The P0 cases use the 3-D octree grid, centered incompressible Navier–Stokes solver, VOF, height-function curvature, and contact-angle boundary from the exact source snapshot in `vof_solver_version.json`.

The official `basilisk.fr` source endpoint returned HTTP 403 through the network proxy. The build therefore used the accessible GitHub mirror at commit `11b5c68731e3facf44b3dfbaaa3b0163f96940b1`. This is a pinned 2019 third-party snapshot; it is not represented as the latest or officially tagged Basilisk release.

## Serial base matrix

On Debian-like systems with GCC, build the pinned Basilisk source and qcc as follows. `-fcommon` is required for this source snapshot with GCC 14.

```sh
git clone https://github.com/tenglongcong/basiliskCFD.git basiliskCFD
cd basiliskCFD
git checkout 11b5c68731e3facf44b3dfbaaa3b0163f96940b1
cd basilisk/src
export BASILISK="$PWD"
gcc -O2 -fcommon -DLIBDIR="\"$BASILISK\"" -c include.c -o include.o
gcc -O2 -fcommon -c rotate.c -o rotate.o
gcc -O2 -fcommon -DLIBDIR="\"$BASILISK\"" \
  -DCC99='"gcc -std=c99 -pedantic -D_GNU_SOURCE=1 -Wno-unused-result"' \
  -DCPP99='""' -DCADNACC='""' -DBASILISK="\"$BASILISK\"" \
  qcc.c include.o rotate.o -o qcc
```

Set `BASILISK_SRC` to this `src` directory and run the serial R/Δx=16/32 matrix plus R/Δx=64 static geometry checks:

```sh
export BASILISK_SRC="$PWD"
python3 run_basilisk_matrix.py
```

The base runner records dynamic exact-curvature spheres, VOF-curvature spheres, 72° sessile drops at R/Δx=16 and 32, time-step repeats at 32, and geometry-only checks at 64. The serial exact-curvature pressure projection is not reliable at R/Δx=64.

## R/Δx=64 dynamic supplement

The dynamic VOF-curvature sphere and sessile-drop cases at R/Δx=64 completed using four Open MPI ranks. Install MPI development headers and runtime (`libopenmpi-dev`, `openmpi-bin`) before building. Set `BASILISK_SRC` to the pinned `src` directory, then run:

```sh
export BASILISK_SRC="/path/to/basiliskCFD/src"
P0_MPI_NP=4 python3 run_basilisk_r64_mpi.py
```

The supplemental runner compiles the same source with qcc's `-D_MPI=1` and `mpicc`, runs the VOF-curvature sphere and 72° sessile drop to t=0.05, stores per-case logs, and upserts their rows into the CSV summaries. `CC99` and `MPIEXEC` can be set by the caller when the compiler or launcher has a non-default name. The R/Δx=32 serial/MPI comparison is stored in `basilisk_cases/diagnostics/mpi_R32_serial_consistency.csv`.

The R/Δx=64 analytical-curvature failure is recorded separately; it is not replaced with the VOF-curvature result. All runs are dimensionless solver verification. They use no Ag material properties and do not predict an Ag melt pool.
