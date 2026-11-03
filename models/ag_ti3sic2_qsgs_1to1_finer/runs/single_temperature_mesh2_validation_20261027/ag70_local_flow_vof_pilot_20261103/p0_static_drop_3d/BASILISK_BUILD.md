# Basilisk build and reproduction

The P0 cases use the 3-D octree grid, centered Navier–Stokes solver, VOF, height-function curvature, and contact-angle boundary from the source snapshot recorded in `vof_solver_version.json`.

The official `basilisk.fr` source endpoint returned HTTP 403 through the current network proxy. The build therefore used the accessible GitHub mirror at commit `11b5c68731e3facf44b3dfbaaa3b0163f96940b1`. This is an explicitly pinned 2019 third-party snapshot; it is not represented as the latest or officially tagged Basilisk release.

On a Debian-like system with GCC and the generated Basilisk C sources (`qcc.c`, `include.c`, and `rotate.c`), build `qcc` as follows. `-fcommon` is needed for this older snapshot with GCC 14.

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

Set `BASILISK_SRC` to the `src` directory and run:

```sh
python3 run_basilisk_matrix.py
```

The script writes one stdout and stderr log per case under `basilisk_cases/`, and the summary tables under `basilisk_results/`. Dynamic cases use `t=0.05`; R/Δx=64 cases are geometry-only at `t=0` because the dynamic projection did not complete within the available run limit. The matrix never uses Ag material properties and is a solver verification, not a molten-Ag prediction.
