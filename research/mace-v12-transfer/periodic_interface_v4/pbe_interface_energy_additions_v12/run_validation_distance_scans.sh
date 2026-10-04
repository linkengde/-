#!/usr/bin/env bash
set -euo pipefail

WORK_ROOT=/workspace/-/research/mace-v12-transfer/periodic_interface_v4
ADD=$WORK_ROOT/pbe_interface_energy_additions_v12
PY=/workspace/.venvs/gpaw-mpi/bin/python
MPI_DEPS=/workspace/.local/gpaw-mpi-deps
ARCHIVE_ROOT=$ADD/calculations/validation_distance_scans
RUN_ROOT=/tmp/pbe_v12_validation_runs
mkdir -p "$ARCHIVE_ROOT" "$RUN_ROOT"

export GPAW_MPI_BACKEND=cgpaw
export GPAW_SETUP_PATH=$($PY -c 'import gpaw_data; print(gpaw_data.datapath())')
export LD_LIBRARY_PATH="$MPI_DEPS/usr/lib/x86_64-linux-gnu/openmpi/lib:$MPI_DEPS/usr/lib/x86_64-linux-gnu:/workspace/.local/gpaw-libs/usr/lib/x86_64-linux-gnu:/workspace/.local/openblas/usr/lib/x86_64-linux-gnu/openblas-pthread${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

for LABEL in AgC_validation_d2p30 AgSi_validation_d2p60 AgTi_validation_d2p60; do
  INPUT=$ADD/inputs/validation_distance_scans/$LABEL.extxyz
  OUT=$RUN_ROOT/$LABEL
  [[ ! -e "$OUT" ]] || { echo "Refusing to reuse existing run directory: $OUT" >&2; exit 2; }
  [[ ! -e "$ARCHIVE_ROOT/$LABEL" ]] || { echo "Refusing to overwrite existing archive: $ARCHIVE_ROOT/$LABEL" >&2; exit 3; }
  echo "Starting $LABEL from SCF iteration 1 on four MPI ranks."
  mpirun --bind-to core --map-by core -n 4 "$PY" \
    "$ADD/run_validation_pw_reference.py" "$LABEL" "$INPUT" "$OUT" \
    >"$ARCHIVE_ROOT/$LABEL.launcher.log" 2>&1
  "$PY" - "$ARCHIVE_ROOT/$LABEL/summary.json" <<'PY'
import json, sys
from pathlib import Path
s=json.loads(Path(sys.argv[1]).read_text())
assert s['scf_converged'] is True and s['mpi_ranks'] == 4 and s['scf_iterations'] > 0
print(f"{s['label']}: SCF PASS in {s['scf_iterations']} iterations; {s['energy_eV_cell']:.10f} eV")
PY
done

"$PY" "$ADD/archive_validation_pw.py" >"$ARCHIVE_ROOT/archive_validation_pw.log"
echo "All three independent validation PW-PBE labels passed archive verification."
