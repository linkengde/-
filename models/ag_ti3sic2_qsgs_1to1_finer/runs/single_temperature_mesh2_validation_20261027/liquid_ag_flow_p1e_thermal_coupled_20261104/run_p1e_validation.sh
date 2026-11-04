#!/usr/bin/env bash
set -euo pipefail

RUN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS="$RUN_DIR/results"
: "${BASILISK_SRC:?Set BASILISK_SRC to the pinned Basilisk src directory}"
export BASILISK="$BASILISK_SRC"
export OPENGLIBS="${OPENGLIBS:--lfb_tiny}"
cd "$RESULTS"

"$BASILISK_SRC/qcc" -disable-dimensions -O2 -Wall -autolink \
  p1e_thermal_marangoni.c -o p1e_thermal_marangoni -lm \
  > build.stdout.log 2> build.stderr.log

run_case() {
  local name="$1"
  shift
  set +e
  timeout 360 ./p1e_thermal_marangoni "$@" \
    > "${name}.stdout.log" 2> "${name}.stderr.log"
  local rc=$?
  set -e
  printf '%s\n' "$rc" > "${name}.exit_code.txt"
  if [[ $rc -ne 0 ]]; then
    echo "${name} failed with exit code ${rc}" >&2
    return "$rc"
  fi
}

# Arguments: cells/radius, thermal Peclet number, gradient flag,
# end time in t0, adaptive mesh flag.
run_case thermal_R8_Pe001_fixed 8 0.01 1 0.5 0
run_case thermal_R8_Pe001_amr 8 0.01 1 0.5 1
run_case thermal_R8_Pe001_zero_gradient 8 0.01 0 0.5 1
run_case thermal_R16_Pe001_amr_half_t0 16 0.01 1 0.5 1

python3 summarize_p1e_logs.py
