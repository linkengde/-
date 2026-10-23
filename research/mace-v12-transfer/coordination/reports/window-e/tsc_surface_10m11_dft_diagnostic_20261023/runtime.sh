#!/usr/bin/env bash
set -euo pipefail
VENV=/workspace/.venvs/gpaw-mpi
GPAW_BIN="$VENV/bin/gpaw"
export PATH="$VENV/bin:$PATH"
export GPAW_BIN
export GPAW_SETUP_PATH="$($VENV/bin/python -c 'from gpaw_data import datapath; print(datapath())')"
export GPAW_MPI_BACKEND=cgpaw
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
