#!/usr/bin/env bash
set -euo pipefail

RUN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS="$RUN_DIR/results"
: "${BASILISK_SRC:?Set BASILISK_SRC to the pinned Basilisk src directory}"
export BASILISK="$BASILISK_SRC"
export OPENGLIBS="${OPENGLIBS:--lfb_tiny}"
export P1F_RESULTS_ROOT="${P1F_RESULTS_ROOT:-$RESULTS/reproduction_$(date -u +%Y%m%d_%H%M%S)}"
mkdir -p "$P1F_RESULTS_ROOT"
if [[ -n "$(find "$P1F_RESULTS_ROOT" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  echo "P1F_RESULTS_ROOT must be empty so existing results are never overwritten: $P1F_RESULTS_ROOT" >&2
  exit 2
fi
P1F_RESULTS_ROOT="$(cd "$P1F_RESULTS_ROOT" && pwd)"
export P1F_RESULTS_ROOT
cd "$RESULTS"

"$BASILISK_SRC/qcc" -disable-dimensions -O2 -Wall -autolink \
  p1f_timestep_sensitivity.c -o p1f_timestep_sensitivity -lm \
  > build_timestep.stdout.log 2> build_timestep.stderr.log
"$BASILISK_SRC/qcc" -disable-dimensions -O2 -Wall -autolink \
  p1f_diffusivity_sensitivity.c -o p1f_diffusivity_sensitivity -lm \
  > build_diffusivity.stdout.log 2> build_diffusivity.stderr.log
"$BASILISK_SRC/qcc" -disable-dimensions -O2 -Wall -autolink \
  p1f_static_balance.c -o p1f_static_balance -lm \
  > build_static.stdout.log 2> build_static.stderr.log

run_case() {
  local case_dir="$1" binary="$2"
  shift 2
  mkdir -p "$P1F_RESULTS_ROOT/$case_dir"
  (
    cd "$P1F_RESULTS_ROOT/$case_dir"
    set +e
    timeout 1800 "$binary" "$@" > run.stdout.log 2> run.stderr.log
    local rc=$?
    printf '%s\n' "$rc" > exit_code.txt
    exit "$rc"
  )
}

BIN_DT="$RESULTS/p1f_timestep_sensitivity"
BIN_PE="$RESULTS/p1f_diffusivity_sensitivity"

# P1F-1: only the max-step cap varies; PeT=0.01 and all other parameters stay fixed.
run_case timestep/dt_1x "$BIN_DT" 1
run_case timestep/dt_0p5x "$BIN_DT" 0.5
run_case timestep/dt_0p25x "$BIN_DT" 0.25
run_case timestep/zero_gradient_dt_0p25x "$BIN_DT" 0.25 0
python3 "$RESULTS/analyze_p1f_timestep.py"

if [[ "$(head -n 1 "$P1F_RESULTS_ROOT/p1f_timestep_decision.txt")" != "PASS" ]]; then
  echo "P1F-1 did not pass; stop before the diffusivity sweep and reduce DT further." >&2
  exit 3
fi

# Supplemental zero-gradient static-drop pressure initialization diagnostics.
# These do not alter the P1F acceptance metrics; they expose local residual flow.
BIN_STATIC="$RESULTS/p1f_static_balance"
run_case timestep/zero_gradient_laplace_pressure "$BIN_STATIC" 1
run_case timestep/zero_gradient_laplace_pressure_negative "$BIN_STATIC" -1
python3 "$RESULTS/analyze_p1f_static.py"

# P1F-2 runs only after P1F-1 passes.  DT cap remains fixed at the P1E baseline.
run_case diffusivity/PeT_0p005 "$BIN_PE" 0.005
run_case diffusivity/PeT_0p01 "$BIN_PE" 0.01
run_case diffusivity/PeT_0p02 "$BIN_PE" 0.02
run_case diffusivity/PeT_0p05 "$BIN_PE" 0.05
python3 "$RESULTS/analyze_p1f_diffusivity.py"
