#!/usr/bin/env bash
set -euo pipefail

# Reuse an existing validated-version MPI build before downloading anything.
check_existing() {
  local py=/workspace/.venvs/gpaw-mpi/bin/python
  [[ -x "$py" ]] && command -v mpirun >/dev/null &&
  "$py" -c 'import gpaw,ase,gpaw_data; assert gpaw.__version__=="26.7.0" and ase.__version__=="3.29.0"' &&
  /workspace/.venvs/gpaw-mpi/bin/gpaw info | rg -q 'MPI enabled +yes' &&
  /workspace/.venvs/gpaw-mpi/bin/gpaw info | rg -q 'scalapack +yes'
}
if check_existing; then
  echo 'Existing GPAW 26.7.0 MPI/ScaLAPACK environment found; reuse it.'
  exit 0
fi
if [[ "${1:-}" == --check ]]; then
  echo 'Required MPI GPAW environment is not ready; install missing dependencies.' >&2
  exit 2
fi

ROOT=/workspace/-/research/mace-v12-transfer
SETUP=/workspace/.setup
LOCAL=/workspace/.local
VENVS=/workspace/.venvs
export PIP_CACHE_DIR=/workspace/.cache/pip
export XDG_CACHE_HOME=/workspace/.cache
export MPLCONFIGDIR=/workspace/.cache/matplotlib
export FC_CACHEDIR=/workspace/.cache/fontconfig
mkdir -p "$PIP_CACHE_DIR" "$MPLCONFIGDIR" "$FC_CACHEDIR" \
  "$SETUP/apt/lists/partial" "$SETUP/apt/archives/partial" \
  "$SETUP/gpaw" "$SETUP/gpaw-mpi" "$LOCAL/gpaw-libs" \
  "$LOCAL/openblas" "$LOCAL/gpaw-mpi-deps" "$VENVS"

cat > "$SETUP/apt/sources.list" <<'SOURCES'
deb [signed-by=/usr/share/keyrings/debian-archive-keyring.gpg] https://deb.debian.org/debian trixie main
SOURCES
apt_args=(-o Dir::Etc::sourcelist="$SETUP/apt/sources.list"
  -o Dir::Etc::sourceparts=-
  -o Dir::State::lists="$SETUP/apt/lists"
  -o Dir::Cache::archives="$SETUP/apt/archives"
  -o APT::Sandbox::User=agent)
apt-get "${apt_args[@]}" update

deb_stage=$(mktemp -d "$SETUP/apt/downloads.XXXXXX")
trap 'rm -rf "$deb_stage"' EXIT
fetch_and_extract() {
  target=$1
  shift
  mkdir -p "$target"
  mkdir "$deb_stage/group"
  (cd "$deb_stage/group" && apt-get "${apt_args[@]}" download "$@")
  for package in "$deb_stage"/group/*.deb; do
    dpkg-deb -x "$package" "$target"
  done
  rm -rf "$deb_stage/group"
}
fetch_and_extract "$LOCAL/gpaw-libs" \
  libxc-dev=5.2.3-3 libxc9=5.2.3-3 libblas-dev=3.12.1-6 libblas3=3.12.1-6
fetch_and_extract "$LOCAL/openblas" \
  libopenblas-pthread-dev=0.3.29+ds-3 libopenblas0-pthread=0.3.29+ds-3
fetch_and_extract "$LOCAL/gpaw-mpi-deps" \
  libopenmpi-dev=5.0.7-1 libscalapack-openmpi-dev=2.2.2-1 \
  libscalapack-openmpi2.2=2.2.2-1 libpmix-dev=5.0.7-1 \
  libibverbs-dev=56.1-1 libltdl-dev=2.5.4-4 \
  libevent-dev=2.1.12-stable-10+b1 \
  libevent-extra-2.1-7t64=2.1.12-stable-10+b1 \
  libevent-openssl-2.1-7t64=2.1.12-stable-10+b1 \
  libhwloc-dev=2.12.0-4 libnl-3-dev=3.7.0-2 \
  libnl-route-3-dev=3.7.0-2 libnuma-dev=2.0.19-1
command -v mpicxx >/dev/null
command -v mpirun >/dev/null
mpicxx --showme:version
mpirun --version | head -n 1

cat > "$SETUP/gpaw-mpi/siteconfig.py" <<'CONFIG'
mpi = True
scalapack = True
compiler = 'mpicxx'
include_dirs = [
    '/workspace/.local/gpaw-libs/usr/include',
    '/workspace/.local/openblas/usr/include/x86_64-linux-gnu/openblas-pthread',
    '/workspace/.local/gpaw-mpi-deps/usr/lib/x86_64-linux-gnu/openmpi/include',
    '/workspace/.local/gpaw-mpi-deps/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi',
]
library_dirs = [
    '/workspace/.local/gpaw-libs/usr/lib/x86_64-linux-gnu',
    '/workspace/.local/openblas/usr/lib/x86_64-linux-gnu/openblas-pthread',
    '/workspace/.local/gpaw-mpi-deps/usr/lib/x86_64-linux-gnu',
    '/workspace/.local/gpaw-mpi-deps/usr/lib/x86_64-linux-gnu/openmpi/lib',
]
libraries = ['xc', 'openblas', 'scalapack-openmpi']
extra_compile_args = ['-fopenmp']
extra_link_args = [
    '-fopenmp',
    '-Wl,-rpath,/workspace/.local/gpaw-libs/usr/lib/x86_64-linux-gnu',
    '-Wl,-rpath,/workspace/.local/openblas/usr/lib/x86_64-linux-gnu/openblas-pthread',
    '-Wl,-rpath,/workspace/.local/gpaw-mpi-deps/usr/lib/x86_64-linux-gnu',
    '-Wl,-rpath,/workspace/.local/gpaw-mpi-deps/usr/lib/x86_64-linux-gnu/openmpi/lib',
]
CONFIG

GPAW_VENV=$VENVS/gpaw-mpi
if [[ ! -x "$GPAW_VENV/bin/python" ]]; then python3.12 -m venv "$GPAW_VENV"; fi
GPAW_PY=$GPAW_VENV/bin/python
"$GPAW_PY" -m pip install 'numpy==2.5.3' 'ase==3.29.0' 'gpaw-data==1.2.1'
if ! ("$GPAW_VENV/bin/gpaw" info 2>/dev/null | rg -q 'MPI enabled +yes' && \
      "$GPAW_VENV/bin/gpaw" info 2>/dev/null | rg -q 'scalapack +yes' && \
      "$GPAW_PY" -c 'import gpaw; assert gpaw.__version__ == "26.7.0"' 2>/dev/null); then
  export GPAW_CONFIG="$SETUP/gpaw-mpi/siteconfig.py" CC=gcc CXX=g++
  "$GPAW_PY" -m pip install --force-reinstall --no-deps --no-binary=gpaw \
    'gpaw==26.7.0' 2>&1 | tee "$SETUP/gpaw-mpi/build.log"
fi

export GPAW_MPI_BACKEND=cgpaw
MPI_DEPS=$LOCAL/gpaw-mpi-deps
export LD_LIBRARY_PATH="$MPI_DEPS/usr/lib/x86_64-linux-gnu/openmpi/lib:$MPI_DEPS/usr/lib/x86_64-linux-gnu:$LOCAL/gpaw-libs/usr/lib/x86_64-linux-gnu:$LOCAL/openblas/usr/lib/x86_64-linux-gnu/openblas-pthread${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export GPAW_SETUP_PATH=$("$GPAW_PY" -c 'import gpaw_data; print(gpaw_data.datapath())')
"$GPAW_VENV/bin/gpaw" info | tee "$SETUP/gpaw-mpi/gpaw-info.txt"
mpirun --bind-to core --map-by core -n 4 "$GPAW_PY" -c \
  'from gpaw.mpi import world; assert world.size == 4; print(f"GPAW MPI ranks: {world.size}", flush=True)'
cat > "$SETUP/gpaw-mpi/pw_smoke.py" <<'PY'
import os
import numpy as np
import gpaw_data
os.environ['GPAW_SETUP_PATH'] = str(gpaw_data.datapath())
from ase import Atoms
from gpaw import GPAW, PW, FermiDirac
from gpaw.mpi import world
atoms = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]], cell=[5, 5, 5], pbc=True)
atoms.center()
atoms.calc = GPAW(mode=PW(200), xc='PBE', kpts=(1, 1, 1),
                  occupations=FermiDirac(0.05),
                  convergence={'energy': 1e-4, 'density': 1e-4, 'eigenstates': 1e-4},
                  maxiter=50, txt='/workspace/.setup/gpaw-mpi/pw-smoke.log',
                  parallel={'domain': world.size})
energy = float(atoms.get_potential_energy())
forces = np.asarray(atoms.get_forces())
assert world.size == 4 and np.isfinite(energy) and np.isfinite(forces).all()
if world.rank == 0:
    print(f'GPAW PW-PBE smoke PASS: ranks={world.size}, energy={energy:.10f} eV', flush=True)
PY
mpirun --bind-to core --map-by core -n 4 "$GPAW_PY" \
  "$SETUP/gpaw-mpi/pw_smoke.py" > "$SETUP/gpaw-mpi/pw-smoke.stdout" 2>&1

echo "GPAW MPI/ScaLAPACK installation and PW-PBE smoke completed. MACE not required for window B."
