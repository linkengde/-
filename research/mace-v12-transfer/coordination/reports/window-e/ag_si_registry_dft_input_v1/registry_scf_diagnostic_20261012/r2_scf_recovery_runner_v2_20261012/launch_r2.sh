#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "usage: launch_r2.sh C_AUTH.json C_AUTH_SHA256 ABSOLUTE_WORK_ROOT" >&2
  exit 2
fi

AUTH_FILE=$1
AUTH_SHA=$2
WORK_ROOT=$3
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PYTHON=${PYTHON:-python}
MPIEXEC=${MPIEXEC:-mpiexec}
GPAW=${GPAW:-gpaw}
command -v timeout >/dev/null || { echo "GNU timeout is required for the 4-hour hard deadline" >&2; exit 2; }
command -v "$GPAW" >/dev/null || { echo "GPAW launcher is unavailable" >&2; exit 2; }
[[ -f "$AUTH_FILE" ]] || { echo "C authorization file is missing" >&2; exit 2; }

case "$WORK_ROOT" in
  /*) ;;
  *) echo "work root must be an absolute path outside the repository" >&2; exit 2 ;;
esac
if [[ -e "$WORK_ROOT" ]]; then
  echo "work root already exists; refusing overwrite or resume" >&2
  exit 2
fi
mkdir -p -- "$(dirname -- "$WORK_ROOT")"

export OMP_NUM_THREADS=1
set +e
timeout --signal=TERM --kill-after=1s 14399s \
  "$MPIEXEC" -n 4 "$GPAW" python -u "$HERE/runner.py" \
  --run \
  --authorization "$AUTH_FILE" \
  --authorization-sha256 "$AUTH_SHA" \
  --work-root "$WORK_ROOT" \
  >"${WORK_ROOT}.stdout.log" 2>"${WORK_ROOT}.stderr.log"
runner_rc=$?
set -e

set +e
"$PYTHON" "$HERE/collector.py" \
  --run-dir "$WORK_ROOT" \
  --launcher-stdout "${WORK_ROOT}.stdout.log" \
  --launcher-stderr "${WORK_ROOT}.stderr.log" \
  --output "$WORK_ROOT/RESULTS_COLLECTED.json"
collector_rc=$?
set -e

if [[ $runner_rc -ne 0 ]]; then
  exit "$runner_rc"
fi
exit "$collector_rc"
