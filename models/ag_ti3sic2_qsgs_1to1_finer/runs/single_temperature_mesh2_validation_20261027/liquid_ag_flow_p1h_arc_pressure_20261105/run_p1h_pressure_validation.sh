#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
P1G_RUN_ROOT="$(cd "$HERE/.." && pwd)/liquid_ag_flow_p1g_exploratory_meltpool_20261104"
MODE="${1:-test}"
BASILISK_SRC="${BASILISK_SRC:-}"
if [[ -z "$BASILISK_SRC" || ! -x "$BASILISK_SRC/qcc" ]]; then
  echo "Set BASILISK_SRC to the pinned Basilisk src directory (v2026-07-20)." >&2
  exit 2
fi
export BASILISK="$BASILISK_SRC"
export OPENGLIBS="${OPENGLIBS:--lfb_tiny}"

build_solver() {
  (cd "$HERE" && "$BASILISK_SRC/qcc" -disable-dimensions -O2 -Wall -autolink \
    p1h_free_surface_coupled.c -o p1h_free_surface_coupled -lm)
  (cd "$HERE" && "$BASILISK_SRC/qcc" -disable-dimensions -O2 -Wall -autolink \
    p1h_pressure_vof_test.c -o p1h_pressure_vof_test -lm)
}

case "$MODE" in
  build)
    build_solver
    ;;
  test)
    cc -std=c99 -Wall -Wextra -O2 "$HERE/test_p1h_arc_pressure.c" \
      -I"$HERE" -lm -o "$HERE/test_p1h_arc_pressure"
    "$HERE/test_p1h_arc_pressure"
    build_solver
    out="$HERE/results/validation_$(date -u +%Y%m%dT%H%M%SZ)"
    mkdir -p "$out"
    (cd "$out" && P1H_TEST_END=0.04 timeout 300 "$HERE/p1h_pressure_vof_test" uniform >uniform.log 2>&1)
    (cd "$out" && P1H_TEST_END=0.04 P1H_TEST_DT=1e-4 timeout 300 "$HERE/p1h_pressure_vof_test" gaussian >gaussian_dt_1x.log 2>&1)
    (cd "$out" && P1H_TEST_END=0.04 P1H_TEST_DT=5e-5 timeout 300 "$HERE/p1h_pressure_vof_test" gaussian >gaussian_dt_0p5x.log 2>&1)
    (cd "$out" && P1H_TEST_END=0.04 P1H_TEST_DT=2.5e-5 timeout 300 "$HERE/p1h_pressure_vof_test" gaussian >gaussian_dt_0p25x.log 2>&1)
    tail -3 "$out/uniform.log"
    tail -3 "$out/gaussian_dt_1x.log"
    tail -3 "$out/gaussian_dt_0p5x.log"
    tail -3 "$out/gaussian_dt_0p25x.log"
    ;;
  smoke)
    build_solver
    input="${P1H_INPUT_DIR:-$P1G_RUN_ROOT/inputs/basilisk}"
    out="${P1H_OUTPUT_DIR:-$HERE/results/geometry_smoke_$(date -u +%Y%m%dT%H%M%SZ)}"
    max_time="${P1H_MAX_TIME_PS:-0.16}"
    mkdir -p "$out"
    (
      cd "$out"
      "$HERE/p1h_free_surface_coupled" "$input" \
        "${P1H_SIGMA_CASE:-negative}" "${P1H_SOLID_DRAG_S:-1e11}" "$max_time"
    )
    ;;
  *)
    echo "Usage: $0 build|test|smoke" >&2
    exit 2
    ;;
esac
