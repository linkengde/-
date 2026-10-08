#!/bin/bash
set -euo pipefail
# Deliberately unconditional; only A may integrate an approved production entry.
echo 'LAUNCH DISABLED: review-only draft, owner jobs not registered' >&2
exit 2
# Below is the integration design, never executable in this report pack.
PACK="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
TASK="${1:?registered window-a or window-b}"
RUN_ROOT="${2:?explicit local run root}"
PY="${GPAW_PYTHON:-/workspace/.venvs/gpaw-mpi/bin/python}"
export GPAW_MPI_BACKEND=cgpaw OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
# A integration must independently check local identity=owner_instance, per-label
# registration, CPU quota>=4, gpaw MPI/scalapack support and existing processes.
mapfile -t LABELS < <(python3 - "$PACK/input_manifest.json" "$TASK" <<'PY'
import json,sys
print('\n'.join(r['label'] for r in json.load(open(sys.argv[1]))['records'] if r['owner_task']==sys.argv[2]))
PY
)
for LABEL in "${LABELS[@]}"; do
  OUT="$RUN_ROOT/$LABEL"
  if [[ -f "$OUT/summary.json" ]]; then
    "$PY" "$PACK/archive_draft.py" "$LABEL" "$OUT" > "$OUT/verification.new.json"
    # Publish only compact verified files after A production integration.
    continue
  fi
  [[ ! -e "$OUT" ]] || { echo "Existing checkpoint/run $OUT: inspect or resume same state; no duplicate launch" >&2; exit 2; }
  INPUT="$(python3 - "$PACK/input_manifest.json" "$LABEL" <<'PY'
import json,sys
print(next(r['input'] for r in json.load(open(sys.argv[1]))['records'] if r['label']==sys.argv[2]))
PY
)"
  # Register running and attach existing watch_progress.py only after owner check.
  mpirun --bind-to core --map-by core -np 4 "$PY" "$PACK/run_pw_reference_draft.py" "$LABEL" "$PACK/$INPUT" "$OUT"
  "$PY" "$PACK/archive_draft.py" "$LABEL" "$OUT" > "$OUT/verification.new.json"
  # Stop on any failed hash/SCF/finite-label check; publish each then continue.
done
