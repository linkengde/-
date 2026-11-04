#!/usr/bin/env bash
set -euo pipefail

RUN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS="$RUN_DIR/results"
: "${BASILISK_SRC:?Set BASILISK_SRC to the pinned Basilisk src directory}"
export BASILISK="$BASILISK_SRC"
export OPENGLIBS="${OPENGLIBS:--lfb_tiny}"
export FROZEN_CONTROL_ROOT="${FROZEN_CONTROL_ROOT:-$RESULTS/frozen_dynamic_repro_$(date -u +%Y%m%d_%H%M%S)}"
mkdir -p "$FROZEN_CONTROL_ROOT"
if [[ -n "$(find "$FROZEN_CONTROL_ROOT" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  echo "FROZEN_CONTROL_ROOT must be empty; refusing to overwrite: $FROZEN_CONTROL_ROOT" >&2
  exit 2
fi
FROZEN_CONTROL_ROOT="$(cd "$FROZEN_CONTROL_ROOT" && pwd)"
export FROZEN_CONTROL_ROOT

cd "$RESULTS"
"$BASILISK_SRC/qcc" -disable-dimensions -O2 -Wall -autolink \
  p1f_static_diagnostic.c -o p1f_static_diagnostic -lm \
  > build_frozen_dynamic.stdout.log 2> build_frozen_dynamic.stderr.log

CASE="$FROZEN_CONTROL_ROOT/r16_amr_frozen_dt025_marangoni"
mkdir -p "$CASE"
(
  cd "$CASE"
  set +e
  timeout 1800 "$RESULTS/p1f_static_diagnostic" 16 2 0.25 0 0.50 1 1 \
    > run.stdout.log 2> run.stderr.log
  rc=$?
  printf '%s\n' "$rc" > exit_code.txt
  exit "$rc"
)
FROZEN_CONTROL_CASE="$CASE" python3 "$RESULTS/analyze_frozen_dynamic_control.py"
