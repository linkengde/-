#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PY=/workspace/.venvs/gpaw-mpi/bin/python
MPI_DEPS=/workspace/.local/gpaw-mpi-deps
RUN_ROOT=/tmp/pbe_v13_holdouts_runs
mkdir -p "$RUN_ROOT" "$ROOT/calculations"

export GPAW_MPI_BACKEND=cgpaw
export GPAW_SETUP_PATH="$($PY -c 'import gpaw_data; print(gpaw_data.datapath())')"
export LD_LIBRARY_PATH="$MPI_DEPS/usr/lib/x86_64-linux-gnu/openmpi/lib:$MPI_DEPS/usr/lib/x86_64-linux-gnu:/workspace/.local/gpaw-libs/usr/lib/x86_64-linux-gnu:/workspace/.local/openblas/usr/lib/x86_64-linux-gnu/openblas-pthread${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

COORD="$ROOT/../../coordination"
WATCHER_PID=""
stop_watcher() {
  if [[ -n "$WATCHER_PID" ]]; then
    kill "$WATCHER_PID" 2>/dev/null || true
    wait "$WATCHER_PID" 2>/dev/null || true
    WATCHER_PID=""
  fi
}
trap stop_watcher EXIT
"$PY" "$COORD/watch_progress.py" window-a "$RUN_ROOT" "$ROOT/calculations" \
  AgC_lateral_registry_holdout_v13 AgSi_lateral_registry_holdout_v13 --parent-pid "$$" \
  >"$RUN_ROOT/progress-publisher.log" 2>&1 &
WATCHER_PID=$!

for LABEL in AgC_lateral_registry_holdout_v13 AgSi_lateral_registry_holdout_v13; do
  INPUT="$ROOT/inputs/$LABEL.extxyz"
  ARCHIVE="$ROOT/calculations/$LABEL"
  OUT="$RUN_ROOT/$LABEL"
  if [[ -d "$ARCHIVE" ]]; then
    "$PY" - "$ARCHIVE/summary.json" <<'PY'
import json, sys
from pathlib import Path
s=json.loads(Path(sys.argv[1]).read_text())
assert s['scf_converged'] is True and s['mpi_ranks'] == 4 and s['scf_iterations'] > 0
print(f"Skipping archived {s['label']}: SCF PASS in {s['scf_iterations']} iterations.")
PY
    continue
  fi
  [[ ! -e "$OUT" ]] || { echo "Refusing to reuse an unarchived run directory: $OUT" >&2; exit 2; }
  echo "Starting $LABEL on four MPI ranks."
  mpirun --bind-to core --map-by core -n 4 "$PY" \
    "$ROOT/run_independent_pw_reference.py" "$LABEL" "$INPUT" "$OUT" \
    >"$ROOT/calculations/$LABEL.launcher.log" 2>&1
done

stop_watcher
"$PY" "$ROOT/archive_independent_holdouts.py"
echo "Independent v13 holdout archive verification PASS."
"$PY" "$ROOT/publish_independent_holdouts.py"
