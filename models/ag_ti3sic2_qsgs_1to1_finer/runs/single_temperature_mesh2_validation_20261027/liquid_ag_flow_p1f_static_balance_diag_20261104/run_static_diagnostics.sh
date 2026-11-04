#!/usr/bin/env bash
set -euo pipefail

RUN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS="$RUN_DIR/results"
: "${BASILISK_SRC:?Set BASILISK_SRC to the pinned Basilisk src directory}"
export BASILISK="$BASILISK_SRC"
export OPENGLIBS="${OPENGLIBS:--lfb_tiny}"
export STATIC_DIAG_RESULTS_ROOT="${STATIC_DIAG_RESULTS_ROOT:-$RESULTS/reproduction_$(date -u +%Y%m%d_%H%M%S)}"
mkdir -p "$STATIC_DIAG_RESULTS_ROOT"
if [[ -n "$(find "$STATIC_DIAG_RESULTS_ROOT" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  echo "STATIC_DIAG_RESULTS_ROOT must be empty; refusing to overwrite: $STATIC_DIAG_RESULTS_ROOT" >&2
  exit 2
fi
STATIC_DIAG_RESULTS_ROOT="$(cd "$STATIC_DIAG_RESULTS_ROOT" && pwd)"
export STATIC_DIAG_RESULTS_ROOT
cd "$RESULTS"

"$BASILISK_SRC/qcc" -disable-dimensions -O2 -Wall -autolink \
  p1f_static_diagnostic.c -o p1f_static_diagnostic -lm \
  > build.stdout.log 2> build.stderr.log

run_case() {
  local name="$1"
  shift
  local out="$STATIC_DIAG_RESULTS_ROOT/$name"
  mkdir -p "$out"
  (
    cd "$out"
    set +e
    timeout 1200 "$RESULTS/p1f_static_diagnostic" "$@" \
      > run.stdout.log 2> run.stderr.log
    local rc=$?
    printf '%s\n' "$rc" > exit_code.txt
    exit "$rc"
  )
}

# Arguments: cells/radius, AMR mode (0 uniform, 1 dynamic, 2 freeze after
# first adaptation, 3 interface/temperature criteria only), timestep factor,
# pressure-jump sign, end time in t0, and pressure-tolerance scale.
run_case r16_amr_dt025_zero 16 1 0.25 0 0.20 1
run_case r16_fixed_dt025_zero 16 0 0.25 0 0.20 1
run_case r16_amr_frozen_dt025_zero 16 2 0.25 0 0.20 1
run_case r16_amr_dynamic_dt0125_zero 16 1 0.125 0 0.20 1
run_case r32_amr_dynamic_dt025_zero 32 1 0.25 0 0.20 1
run_case r32_amr_frozen_dt025_zero 32 2 0.25 0 0.10 1
run_case r16_amr_dynamic_dt025_plus 16 1 0.25 1 0.20 1
run_case r16_amr_dynamic_dt025_minus 16 1 0.25 -1 0.20 1
run_case r16_amr_dynamic_tol1e-5 16 1 0.25 0 0.20 0.1
run_case r32_amr_dynamic_tol1e-5 32 1 0.25 0 0.10 0.1
run_case r16_amr_interface_only_dt025_zero 16 3 0.25 0 0.20 1

python3 "$RESULTS/summarize_static_diagnostics.py"
