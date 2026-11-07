#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 4 || $# -gt 6 ]]; then
  echo "usage: BASILISK_SRC=/path/to/basilisk/src $0 METHOD INPUT_DIR OUTPUT_DIR DT_MAX_FS=1 MAX_TIME_PS=0.0004 OMP_THREADS=1" >&2
  exit 2
fi
method="$1"; input_dir="$(cd "$2" && pwd)"; output_dir="$3"
dt_max_fs="${4:-1}"; max_time_ps="${5:-0.0004}"; omp_threads="${6:-1}"
if [[ -z "${BASILISK_SRC:-}" || ! -x "$BASILISK_SRC/qcc" ]]; then
  echo "Set BASILISK_SRC to the fixed Basilisk src directory containing qcc." >&2
  exit 2
fi
script_dir="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$output_dir"; output_dir="$(cd "$output_dir" && pwd)"
if find "$output_dir" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
  echo "Refusing to overwrite non-empty output directory: $output_dir" >&2
  exit 2
fi
tmp="$(mktemp -d "${TMPDIR:-/tmp}/p1p-solver-trial.XXXXXX")"
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/build" "$tmp/run"
cp "$script_dir/solver_trial_source.c" "$tmp/build/thermal_trial.c"
method_define=()
case "$method" in
  line-search) method_define=(-DP1K_RESIDUAL_LINE_SEARCH) ;;
  enthalpy-residual) method_define=(-DP1K_ENTHALPY_RESIDUAL) ;;
  fixed-point) method_define=(-DP1I_FULL_THERMAL_RESIDUAL_DIAGNOSTICS) ;;
  *) echo "unknown method: $method" >&2; exit 2 ;;
esac
(
  cd "$tmp/build"
  BASILISK="$BASILISK_SRC" OPENGLIBS=-lfb_tiny \
    "$BASILISK_SRC/qcc" -O2 -Wall -autolink -disable-dimensions -fopenmp \
    -DP1H_RECT_NZ=76 -DP1H_TOP_OPEN=1 -DP1J_ACTIVE_THERMAL_DIAG \
    -DP1H_THERMAL_NL_DIAGNOSTICS -DP1I_PREFLIGHT "${method_define[@]}" \
    -I"$script_dir/local_headers" -I. thermal_trial.c -o thermal_trial -lm
)
(
  cd "$tmp/run"
  set +e
  OMP_NUM_THREADS="$omp_threads" P1M_SOURCE_POWER_W=0.004 \
    P1J_DT_MAX_FS="$dt_max_fs" P1H_MARANGONI=off \
    timeout 600 "$tmp/build/thermal_trial" "$input_dir" negative 1e11 "$max_time_ps" \
      > run.stdout 2> run.stderr
  status=$?
  set -e
  cp run.stdout run.stderr "$output_dir/"
  for f in p1k_line_search.csv p1k_enthalpy_iteration.csv thermal_nl_iterations.csv p1h_history.csv; do
    [[ ! -f "$f" ]] || cp "$f" "$output_dir/"
  done
  if grep -q "STOP nonlinear enthalpy solve failed" run.stderr; then
    echo nonlinear_failure > "$output_dir/solver_status.txt"
  elif [[ $status -eq 124 ]]; then
    echo timeout > "$output_dir/solver_status.txt"
  elif [[ $status -eq 0 ]]; then
    echo reached_configured_stop > "$output_dir/solver_status.txt"
  else
    echo "process_exit_$status" > "$output_dir/solver_status.txt"
  fi
)
