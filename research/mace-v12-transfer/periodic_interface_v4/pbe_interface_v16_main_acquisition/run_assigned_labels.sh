#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
COORD="$ROOT/../../coordination/sync_tasks.py"
PY="${GPAW_PYTHON:-/workspace/.venvs/gpaw-mpi/bin/python}"
TASK=window-a
[[ -x "$PY" ]] || { echo "Missing GPAW Python; see coordination/environment_setup.md." >&2; exit 2; }
python3 "$COORD" claim "$TASK"
MPI_DEPS="${MPI_DEPS:-/workspace/.local/gpaw-mpi-deps}"
export GPAW_MPI_BACKEND=cgpaw
export LD_LIBRARY_PATH="$MPI_DEPS/usr/lib/x86_64-linux-gnu/openmpi/lib:$MPI_DEPS/usr/lib/x86_64-linux-gnu:/workspace/.local/gpaw-libs/usr/lib/x86_64-linux-gnu:/workspace/.local/openblas/usr/lib/x86_64-linux-gnu/openblas-pthread${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export GPAW_SETUP_PATH="$("$PY" -c 'import gpaw_data; print(gpaw_data.datapath())')"
"$PY" - <<'PY'
import os
from pathlib import Path
available = len(os.sched_getaffinity(0))
p = Path('/sys/fs/cgroup/cpu.max')
if p.exists():
    quota, period = p.read_text().split()
    if quota != 'max': available = min(available, int(quota) // int(period))
if available < 4: raise SystemExit('Four actual CPU cores are required; do not oversubscribe.')
import gpaw
assert gpaw.__version__ == '26.7.0'
PY
"$(dirname "$PY")/gpaw" info | rg -q 'MPI enabled +yes'
"$(dirname "$PY")/gpaw" info | rg -q 'scalapack +yes'
RUN_ROOT="${DFT_RUN_ROOT:-/workspace/mace_v16_main_window-a}"
mkdir -p "$RUN_ROOT" "$ROOT/calculations"
mapfile -t LABELS < <(python3 -c 'import json,sys; print("\n".join(r["label"] for r in json.load(open(sys.argv[1]))["records"]))' "$ROOT/input_manifest.json")
python3 "$COORD" publish "$TASK"
python3 "$ROOT/../../coordination/watch_progress.py" "$TASK" "$RUN_ROOT" "$ROOT/calculations" "${LABELS[@]}" --parent-pid "$$" >"/tmp/window-a-v16-watcher.log" 2>&1 &
WATCH_PID=$!
cleanup() { kill "$WATCH_PID" 2>/dev/null || true; }
trap cleanup EXIT
for LABEL in "${LABELS[@]}"; do
  LABEL_STATUS="$(python3 - "$ROOT/../../coordination/tasks/$TASK.json" "$LABEL" <<'PY'
import json, sys
task = json.load(open(sys.argv[1]))
print(task.get('labels', {}).get(sys.argv[2], 'queued'))
PY
)"
  if [[ "$LABEL_STATUS" == "completed" ]]; then
    echo "Skipping already published label $LABEL."
    continue
  fi
  ARCHIVE="$ROOT/calculations/$LABEL"
  OUT="$RUN_ROOT/$LABEL"
  if [[ -f "$ARCHIVE/summary.json" ]]; then
    "$PY" "$ROOT/archive_parallel_pw.py" "$LABEL"
    python3 "$COORD" publish "$TASK" --label "$LABEL"
    continue
  fi
  [[ ! -e "$OUT" ]] || { echo "Unarchived run exists: $OUT. Inspect its process/checkpoint; do not start a duplicate." >&2; exit 2; }
  echo "Starting assigned v16 training acquisition $LABEL with four MPI ranks."
  mpirun --bind-to core --map-by core -n 4 "$PY" "$ROOT/run_pw_reference.py" \
    "$LABEL" "$ROOT/inputs/$LABEL.extxyz" "$OUT" >"$ROOT/calculations/$LABEL.launcher.log" 2>&1
  "$PY" "$ROOT/archive_parallel_pw.py" "$LABEL"
  python3 "$COORD" publish "$TASK" --label "$LABEL"
done
cleanup
trap - EXIT
echo "v16 acquisition queue completed and published."
