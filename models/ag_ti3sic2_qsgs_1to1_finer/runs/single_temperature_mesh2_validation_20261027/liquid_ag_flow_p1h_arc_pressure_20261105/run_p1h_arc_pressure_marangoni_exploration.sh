#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(git -C "$HERE" rev-parse --show-toplevel)"
: "${BASILISK_SRC:?Set BASILISK_SRC to the pinned Basilisk src directory}"
INPUT_DIR="${P1H_INPUT_DIR:-$REPO_ROOT/models/ag_ti3sic2_qsgs_1to1_finer/runs/single_temperature_mesh2_validation_20261027/liquid_ag_flow_p1g_exploratory_meltpool_20261104/inputs/basilisk}"
OUT_ROOT="${P1H_EXPLORATION_OUTPUT_DIR:-/workspace/p1h_arc_pressure_exploration_20261105/production_pair}"
BUILD_ROOT="${P1H_EXPLORATION_BUILD_DIR:-/workspace/p1h_arc_pressure_exploration_20261105/build}"
OMP_NUM_THREADS="${OMP_NUM_THREADS:-5}"
SEGMENT_PS="${P1H_SEGMENT_PS:-2}"
MAX_TOTAL_PS="${P1H_MAX_TOTAL_PS:-20}"
PRESSURE_PA="${P1H_EXPLORATION_PRESSURE_PA:-1800000}"
PRESSURE_RADIUS_NM="${P1H_EXPLORATION_PRESSURE_RADIUS_NM:-5}"

[[ -d "$BASILISK_SRC" ]] || { echo "BASILISK_SRC not found: $BASILISK_SRC" >&2; exit 2; }
[[ -d "$INPUT_DIR" ]] || { echo "P1H input directory not found: $INPUT_DIR" >&2; exit 2; }
mkdir -p "$OUT_ROOT" "$BUILD_ROOT"

EXE="$BUILD_ROOT/p1h_nz76_open_marangoni"
(cd "$HERE" && "$BASILISK_SRC/qcc" -DP1H_RECT_NZ=76 -DP1H_TOP_OPEN=1 \
  -DPOOL_REPORT_INTERVAL_PS=0.05 -DPOOL_FIELD_INTERVAL_PS=0.25 \
  -DPOOL_CHECKPOINT_INTERVAL_PS=0.5 -disable-dimensions -fopenmp \
  -O2 -Wall -autolink p1h_free_surface_rectangular.c -o "$EXE" -lm)

compress_completed_artifacts() {
  local out="$1" latest=""
  [[ -f "$out/latest_checkpoint.txt" ]] && latest="$(cat "$out/latest_checkpoint.txt")"
  while IFS= read -r -d '' f; do gzip -n -f "$f"; done \
    < <(find "$out" -maxdepth 1 -type f -name '*.bin' -print0)
  while IFS= read -r -d '' f; do
    local base
    base="$(basename "$f")"
    [[ "$base" == "$latest" ]] && continue
    if [[ -f "$f.complete" ]]; then
      gzip -n -f "$f"
      # A new checkpoint may be published after `latest` was read above. If
      # this compression raced that pointer update, keep the recovery pointer
      # valid by switching it to the completed gzip artifact atomically.
      local latest_now
      latest_now="$(cat "$out/latest_checkpoint.txt" 2>/dev/null || true)"
      if [[ "$latest_now" == "$base" ]]; then
        printf '%s.gz\n' "$base" >"$out/latest_checkpoint.txt.part"
        mv "$out/latest_checkpoint.txt.part" "$out/latest_checkpoint.txt"
      fi
    fi
  done < <(find "$out" -maxdepth 1 -type f -name '*.dump' -print0)
}

compress_while_running() {
  local out="$1"
  while [[ ! -f "$out/.case_finished" ]]; do
    compress_completed_artifacts "$out"
    sleep 5
  done
  compress_completed_artifacts "$out"
  while IFS= read -r -d '' f; do
    [[ -f "$f.complete" ]] && continue
    gzip -n -f "$f"
  done < <(find "$out" -maxdepth 1 -type f -name '*.dump' -print0)
}

run_case() {
  local name="$1" pressure_enabled="$2" profile="$3" out="$OUT_ROOT/$1"
  local configured_peak="$PRESSURE_PA"
  [[ "$pressure_enabled" == 1 ]] || configured_peak=0
  [[ ! -e "$out/run_config.env" ]] || { echo "Refusing to overwrite $out" >&2; return 2; }
  mkdir -p "$out"
  cat >"$out/run_config.env" <<EOF
grid=64x64x76; sample geometry unchanged; 23-cell gas cap
top_boundary=open, ambient gauge pressure 0
heating=existing Ag70 source, 8.93189911185 mW, 0.02 ps smooth ramp
marangoni=on, existing negative-dsigma/dT candidate (exploratory)
pressure_enabled=$pressure_enabled
pressure_profile=$profile
pressure_peak_Pa=$configured_peak
pressure_radius_nm=$PRESSURE_RADIUS_NM
pressure_center=existing heat-source center
pressure_and_heat_off_ramp_ps=0.02
heat_stop=Ti3SiC2 reaches 2300 K OR source/sink mismatch <=0.01 mW for 10 consecutive steps
cooling_stop=two consecutive steps with no liquid Ag after the source and pressure ramps end
dt_max_fs=20, adaptively reduced approaching 2300 K
segment_ps=$SEGMENT_PS
maximum_time_guard_ps=$MAX_TOTAL_PS
EOF
  compress_while_running "$out" &
  local compressor_pid=$!
  local restart_dump="" end_ps="$SEGMENT_PS" chunk=0 done_case=0
  while (( $(awk -v t="$end_ps" -v m="$MAX_TOTAL_PS" 'BEGIN {print (t<=m)?1:0}') )); do
    chunk=$((chunk + 1))
    local stderr="$out/chunk_$(printf '%03d' "$chunk").stderr.log"
    local stdout="$out/chunk_$(printf '%03d' "$chunk").stdout.log"
    set +e
    if [[ -z "$restart_dump" ]]; then
      (cd "$out" && OMP_NUM_THREADS="$OMP_NUM_THREADS" OMP_DYNAMIC=FALSE \
        P1H_MARANGONI=on P1H_ARC_PRESSURE_ENABLED="$pressure_enabled" \
        P1H_ARC_PRESSURE_PEAK_PA="$configured_peak" \
        P1H_ARC_PRESSURE_PROFILE="$profile" \
        P1H_ARC_PRESSURE_RADIUS_NM="$PRESSURE_RADIUS_NM" \
        P1H_ARC_PRESSURE_RAMP_PS=0.02 P1H_SOURCE_RAMP_PS=0.02 \
        "$EXE" "$INPUT_DIR" negative 1e11 "$end_ps" >"$stdout" 2>"$stderr")
      rc=$?
    else
      (cd "$out" && OMP_NUM_THREADS="$OMP_NUM_THREADS" OMP_DYNAMIC=FALSE \
        P1H_MARANGONI=on P1H_ARC_PRESSURE_ENABLED="$pressure_enabled" \
        P1H_ARC_PRESSURE_PEAK_PA="$configured_peak" \
        P1H_ARC_PRESSURE_PROFILE="$profile" \
        P1H_ARC_PRESSURE_RADIUS_NM="$PRESSURE_RADIUS_NM" \
        P1H_ARC_PRESSURE_RAMP_PS=0.02 P1H_SOURCE_RAMP_PS=0.02 \
        P1H_RESTART_DUMP="$restart_dump" \
        "$EXE" "$INPUT_DIR" negative 1e11 "$end_ps" >"$stdout" 2>"$stderr")
      rc=$?
    fi
    set -e
    rg 'INIT |RESTARTED|HEAT_OFF_AT_|HEAT_OFF_PREDICTED_TARGET|DONE_SOLIDIFIED|MAX_TIME_STOP|STOP ' "$stderr" | tail -n 10 || true
    if rg -q 'DONE_SOLIDIFIED' "$stderr"; then
      done_case=1
      break
    fi
    if rg -q 'MAX_TIME_STOP' "$stderr"; then
      local latest
      latest="$(cat "$out/latest_checkpoint.txt")"
      if [[ "$latest" == *.gz ]]; then
        gzip -dk "$out/$latest"
        latest="${latest%.gz}"
        printf '%s\n' "$latest" >"$out/latest_checkpoint.txt"
      fi
      restart_dump="$out/$latest"
      end_ps="$(awk -v t="$end_ps" -v d="$SEGMENT_PS" 'BEGIN {printf "%.9g", t+d}')"
      continue
    fi
    echo "Case $name failed in chunk $chunk with exit=$rc; see $stderr" >&2
    touch "$out/.case_failed"
    break
  done
  touch "$out/.case_finished"
  wait "$compressor_pid"
  if [[ "$done_case" != 1 ]]; then
    echo "Case $name reached the configured maximum-time guard without a solidification event; latest checkpoint retained." >&2
    return 1
  fi
}

if [[ "${P1H_RUN_ZERO_PRESSURE:-1}" == 1 ]]; then
  run_case marangoni_on_zero_pressure 0 uniform
fi
if [[ "${P1H_RUN_GAUSSIAN_PRESSURE:-1}" == 1 ]]; then
  pressure_mpa="$(awk -v p="$PRESSURE_PA" 'BEGIN {printf "%.6g",p/1.e6}')"
  pressure_tag="${pressure_mpa//./p}"
  run_case "marangoni_on_gaussian_${pressure_tag}MPa" 1 periodic_gaussian
fi

echo "Requested exploratory cases reached full solidification. Results: $OUT_ROOT"
