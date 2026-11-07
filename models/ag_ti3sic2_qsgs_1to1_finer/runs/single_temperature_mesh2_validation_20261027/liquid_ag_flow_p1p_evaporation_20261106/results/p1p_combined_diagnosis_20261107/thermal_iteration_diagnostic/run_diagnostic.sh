#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "usage: BASILISK_SRC=/path/to/basilisk/src $0 INPUT_DIR OUTPUT_DIR" >&2
  exit 2
fi
if [[ -z "${BASILISK_SRC:-}" || ! -x "$BASILISK_SRC/qcc" ]]; then
  echo "Set BASILISK_SRC to the fixed Basilisk src directory containing qcc." >&2
  exit 2
fi
input_dir="$(cd "$1" && pwd)"
output_dir="$2"
mkdir -p "$output_dir"
output_dir="$(cd "$output_dir" && pwd)"
if [[ -e "$output_dir/thermal_nl_iterations.csv" ]]; then
  echo "Refusing to overwrite existing diagnostic output: $output_dir/thermal_nl_iterations.csv" >&2
  exit 2
fi
script_dir="$(cd "$(dirname "$0")" && pwd)"
tmp="$(mktemp -d "${TMPDIR:-/tmp}/p1p-thermal-iteration-diag.XXXXXX")"
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/build" "$tmp/run"
cp "$script_dir/diagnostic_source.c" "$tmp/build/thermal_diag.c"

(
  cd "$tmp/build"
  BASILISK="$BASILISK_SRC" OPENGLIBS=-lfb_tiny \
    "$BASILISK_SRC/qcc" -O2 -Wall -autolink -disable-dimensions -fopenmp \
    -DP1H_RECT_NZ=76 -DP1H_TOP_OPEN=1 -DP1J_ACTIVE_THERMAL_DIAG \
    -DP1H_THERMAL_NL_DIAGNOSTICS -DP1I_FULL_THERMAL_RESIDUAL_DIAGNOSTICS \
    -DP1I_PREFLIGHT -I"$script_dir/local_headers" -I. \
    thermal_diag.c -o "$tmp/thermal_diag" -lm
)
(
  cd "$tmp/run"
  export OMP_NUM_THREADS=2
  export P1M_SOURCE_POWER_W=0.004
  export P1J_DT_MAX_FS=1
  export P1H_MARANGONI=off
  set +e
  timeout 600 "$tmp/thermal_diag" "$input_dir" negative 1e11 0.0004 > run.stdout 2> run.stderr
  status=$?
  set -e
  printf '%s\n' "$status" > exit_status.txt
  cp thermal_nl_iterations.csv "$output_dir/"
  cp run.stdout "$output_dir/thermal_diag_replay.stdout"
  cp run.stderr "$output_dir/thermal_diag_replay.stderr"
  cp exit_status.txt "$output_dir/"
  exit "$status"
)
