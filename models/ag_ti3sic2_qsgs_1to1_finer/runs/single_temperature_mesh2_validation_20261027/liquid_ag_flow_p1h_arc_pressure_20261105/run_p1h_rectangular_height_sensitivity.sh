#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(git -C "$HERE" rev-parse --show-toplevel)"
: "${BASILISK_SRC:?Set BASILISK_SRC to the pinned Basilisk src directory}"
INPUT_DIR="${P1G_INPUT_DIR:-$REPO_ROOT/models/ag_ti3sic2_qsgs_1to1_finer/runs/single_temperature_mesh2_validation_20261027/liquid_ag_flow_p1g_exploratory_meltpool_20261104/inputs/basilisk}"
OUT_ROOT="${P1H_RECT_OUTPUT_DIR:-$HERE/results/rectangular_height_sensitivity_20261105}"
BUILD_ROOT="${P1H_RECT_BUILD_DIR:-/tmp/p1h-rect-build-20261105}"
OMP_NUM_THREADS="${OMP_NUM_THREADS:-5}"
MAX_TIME_PS="${P1H_RECT_MAX_TIME_PS:-0.16}"

[[ -d "$BASILISK_SRC" ]] || { echo "BASILISK_SRC not found: $BASILISK_SRC" >&2; exit 2; }
[[ -d "$INPUT_DIR" ]] || { echo "P1G input directory not found: $INPUT_DIR" >&2; exit 2; }
if [[ -e "$OUT_ROOT/nz64_closed/run_config.env" ]]; then
  echo "Output already exists: $OUT_ROOT (choose a new P1H_RECT_OUTPUT_DIR to preserve prior results)" >&2
  exit 2
fi
mkdir -p "$OUT_ROOT" "$BUILD_ROOT"

compile_case() {
  local nz="$1" boundary="$2" exe="$BUILD_ROOT/p1h_rect_nz${nz}_${boundary}"
  local -a flags=("-DP1H_RECT_NZ=$nz" -disable-dimensions -fopenmp -O2 -Wall -autolink)
  [[ "$boundary" == open ]] && flags+=("-DP1H_TOP_OPEN=1")
  (cd "$HERE" && "$BASILISK_SRC/qcc" "${flags[@]}" p1h_free_surface_rectangular.c -o "$exe" -lm)
  printf '%s\n' "$exe"
}

run_case() {
  local nz="$1" boundary="$2" exe="$3"
  local out="$OUT_ROOT/nz${nz}_${boundary}"
  mkdir -p "$out"
  cat >"$out/run_config.env" <<EOF
grid=64x64x${nz}
sample_cells_z=53
gas_cap_cells=$((nz - 53))
top_boundary=$boundary
marangoni=off
thermal_source=existing Ag70 source, unchanged
pressure=uniform 1.8 MPa software-path smoke only; not calibrated nanoscale arc pressure
max_time_ps=$MAX_TIME_PS
openmp_threads=$OMP_NUM_THREADS
EOF
  set +e
  (cd "$out" && \
    OMP_NUM_THREADS="$OMP_NUM_THREADS" OMP_DYNAMIC=FALSE \
    P1H_MARANGONI=off \
    P1H_ARC_PRESSURE_ENABLED=1 \
    P1H_ARC_PRESSURE_PEAK_PA=1800000 \
    P1H_ARC_PRESSURE_PROFILE=uniform \
    P1H_ARC_PRESSURE_RAMP_PS=0.02 \
    P1H_SOURCE_RAMP_PS=0.02 \
    "$exe" "$INPUT_DIR" negative 1e11 "$MAX_TIME_PS" \
    >run.stdout.log 2>run.stderr.log)
  local rc=$?
  set -e
  if [[ $rc -ne 0 ]] && ! rg -q "MAX_TIME_STOP t_ps=" "$out/run.stderr.log"; then
    echo "case nz=${nz} boundary=${boundary} failed (exit $rc); see $out/run.stderr.log" >&2
    return "$rc"
  fi
  # The solver intentionally returns 1 from its stop event at max_time.
  # Keep the latest restart dump directly usable; compress other artifacts.
  local latest
  latest="$(cat "$out/latest_checkpoint.txt")"
  find "$out" -maxdepth 1 -type f \( -name '*.bin' -o -name '*.dump' \) -print0 |
    while IFS= read -r -d '' f; do
      [[ "$(basename "$f")" == "$latest" ]] && continue
      gzip -n "$f"
    done
  # Archive the latest checkpoint too after the run ends; retain its pointer
  # with the .gz suffix so it can be restored after `gzip -dk`.
  if [[ -f "$out/$latest" ]]; then
    gzip -n "$out/$latest"
    printf '%s.gz\n' "$latest" >"$out/latest_checkpoint.txt"
  fi
  printf 'COMPLETED nz=%s boundary=%s exit=%s out=%s\n' "$nz" "$boundary" "$rc" "$out"
}

declare -A bins
for nz in 64 76 86; do
  for boundary in closed open; do
    bins["${nz}_${boundary}"]="$(compile_case "$nz" "$boundary")"
  done
done

for nz in 64 76 86; do
  for boundary in closed open; do
    run_case "$nz" "$boundary" "${bins[${nz}_${boundary}]}"
  done
done

python3 "$HERE/summarize_p1h_rectangular_height_sensitivity.py" "$OUT_ROOT"
