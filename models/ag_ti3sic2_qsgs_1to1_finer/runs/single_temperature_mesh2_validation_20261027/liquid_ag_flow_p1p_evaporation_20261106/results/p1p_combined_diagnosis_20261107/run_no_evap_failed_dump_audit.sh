#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "usage: BASILISK_SRC=/path/to/basilisk/src $0 /path/to/cell_audit.csv" >&2
  exit 2
fi
if [[ -z "${BASILISK_SRC:-}" || ! -x "$BASILISK_SRC/qcc" ]]; then
  echo "Set BASILISK_SRC to the Basilisk src directory containing qcc." >&2
  exit 2
fi

script_dir="$(cd "$(dirname "$0")" && pwd)"
run_dir="$(cd "$script_dir/../.." && pwd)"
diag_tmp="$(mktemp -d)"
source_tmp="$run_dir/p1p/inspect_no_evap_failed_dump_tmp.c"
solver_src_tmp="$run_dir/p1p/p1m_thermal_iteration_diagnostic_audit_tmp.c"
trap 'rm -rf "$diag_tmp"; rm -f "$source_tmp" "$solver_src_tmp"' EXIT

gzip -dc "$run_dir/tests/solver_smoke/original_solver_control/original_p1m_no_evap_failed.dump.gz" > "$diag_tmp/failed.dump"

# qcc rewrites Basilisk .c input into a preprocessed form in-place. Compile
# disposable copies so neither the tracked smoke source nor this reader is
# touched. qcc resolves local solver headers from the p1p source directory.
cp "$run_dir/tests/solver_smoke/original_solver_control/p1m_thermal_iteration_diagnostic.c" \
  "$solver_src_tmp"
sed 's#../../tests/solver_smoke/original_solver_control/p1m_thermal_iteration_diagnostic.c#p1m_thermal_iteration_diagnostic_audit_tmp.c#' \
  "$script_dir/inspect_no_evap_failed_dump.c" > "$source_tmp"

(
  cd "$run_dir/p1p"
  BASILISK="$BASILISK_SRC" OPENGLIBS=-lfb_tiny \
    "$BASILISK_SRC/qcc" -O2 -Wall -autolink -disable-dimensions -fopenmp \
    -DP1H_RECT_NZ=76 -DP1H_TOP_OPEN=1 -I. \
    inspect_no_evap_failed_dump_tmp.c -o "$diag_tmp/inspect_dump" -lm
  "$diag_tmp/inspect_dump" \
    "$diag_tmp/failed.dump" \
    "$script_dir/inputs" \
    "$1"
)
