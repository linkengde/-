#!/bin/bash
set -euo pipefail
cd /workspace/-
PACK=research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v17_rank_recovery_parallel
LABEL=Ti_gap_2p3_rank_recovery_v17_01
RUN=/workspace/mace_v17_rank_recovery_frame0
if [ -e "$RUN/result" ] || [ -e "$PACK/calculations/$LABEL" ]; then
  echo 'Existing result/checkpoint; refusing duplicate launch' >&2
  exit 1
fi
mkdir -p "$RUN"
export GPAW_MPI_BACKEND=cgpaw OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
mpirun -np 4 /workspace/.venvs/gpaw-mpi/bin/python "$PACK/run_pw_reference.py" "$LABEL" "$PACK/inputs/$LABEL.extxyz" "$RUN/result" > "$RUN/mpi.log" 2>&1
/workspace/.venvs/gpaw-mpi/bin/python "$PACK/archive.py" > "$RUN/archive.log" 2>&1
python3 research/mace-v12-transfer/coordination/sync_tasks.py publish window-b
python3 research/mace-v12-transfer/coordination/sync_tasks.py progress window-b --job v17_rank_recovery_frame0_DFT --state completed --iteration 1 --note "Exact frame0 four-rank relabel converged, compact archive verified and published; Gamma kpoint convergence remains unresolved."
