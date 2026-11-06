#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(git -C "$HERE" rev-parse --show-toplevel)"
: "${BASILISK_SRC:?Set BASILISK_SRC to the pinned Basilisk src directory}"
INPUT_DIR="${P1I_INPUT_DIR:-$REPO_ROOT/models/ag_ti3sic2_qsgs_1to1_finer/runs/single_temperature_mesh2_validation_20261027/liquid_ag_flow_p1g_exploratory_meltpool_20261104/inputs/basilisk}"
OUT_DIR="${P1I_RUN_DIR:-/workspace/p1h_arc_pressure_exploration_20261105/p1i_3x_hold3ps_20261106}"
BUILD_DIR="${P1I_BUILD_DIR:-/tmp/p1i-build}"
OMP_THREADS="${P1I_OMP_THREADS:-4}"
MAX_TIME_PS="${P1I_MAX_TIME_PS:-10}"
PRESSURE_PA=501000000
PRESSURE_RADIUS_NM=5
MAX_SAVED_CHECKPOINTS=2
SOURCE_POWER_MW="${P1I_SOURCE_POWER_MW:-8.93189911185}"
SOURCE_POWER_W="$(awk -v p="$SOURCE_POWER_MW" 'BEGIN { printf "%.15g", p/1000 }')"

[[ -d "$BASILISK_SRC" ]] || { echo "BASILISK_SRC not found: $BASILISK_SRC" >&2; exit 2; }
[[ -d "$INPUT_DIR" ]] || { echo "P1G input directory not found: $INPUT_DIR" >&2; exit 2; }
[[ ! -e "$OUT_DIR/run_config.env" ]] || { echo "Refusing to overwrite existing run: $OUT_DIR" >&2; exit 2; }
mkdir -p "$OUT_DIR" "$BUILD_DIR"
EXE="$BUILD_DIR/p1i_nz76_open_marangoni3x"
(cd "$HERE" && "$BASILISK_SRC/qcc" -DP1H_RECT_NZ=76 -DP1H_TOP_OPEN=1 \
  -DP1I_SOURCE_POWER_W="$SOURCE_POWER_W" \
  -DPOOL_REPORT_INTERVAL_PS=0.05 -DPOOL_FIELD_INTERVAL_PS=20.0 \
  -DPOOL_CHECKPOINT_INTERVAL_PS=0.5 -disable-dimensions -fopenmp \
  -O2 -Wall -autolink p1i_free_surface_rectangular.c -o "$EXE" -lm)

cat >"$OUT_DIR/run_config.env" <<CFG
git_branch=$(git -C "$REPO_ROOT" branch --show-current)
git_head=$(git -C "$REPO_ROOT" rev-parse HEAD)
grid=64x64x76; sample geometry unchanged; Ag70; 23-cell gas cap
top_boundary=open, ambient gauge pressure 0
thermal_model=single-temperature enthalpy; existing composite heat source; source power=${SOURCE_POWER_MW}mW
source_template_integral=8.93189911185mW; VOF-mapped source is renormalized to the requested target
marangoni=on; tangential multiplier=3; normal capillary stress unchanged
arc_pressure=periodic_gaussian; configured peak=501MPa; radius=5nm; center=existing source center
heating_target=Ti3SiC2 max 2100K; then full 3ps hold; heat and pressure off together after hold
pre_target_balance_stop=abs(input-bottom)<=0.01mW for 10 consecutive steps
safety=hard 2300K table ceiling with predictive off-ramp; divergence guard retained
cooling_stop=two successive steps with no liquid Ag after source and pressure ramps end
max_time_guard_ps=$MAX_TIME_PS; dt_max=20fs; source and pressure ramp=0.02ps
openmp_threads=$OMP_THREADS
CFG
sha256sum "$HERE/p1i_free_surface_rectangular.c" "$HERE/p1h_arc_pressure.h" \
  "$HERE/two-phase-clsvof-p1g.h" >"$OUT_DIR/code_hashes.sha256"

compress_completed_artifacts() {
  local latest latest_base latest_index="0"
  latest="$(cat "$OUT_DIR/latest_checkpoint.txt" 2>/dev/null || true)"
  latest_base="${latest%.gz}"
  if [[ "$latest_base" =~ restart_i([0-9]+)_ ]]; then latest_index="${BASH_REMATCH[1]}"; fi

  while IFS= read -r -d '' f; do gzip -1 -n -f "$f"; done \
    < <(find "$OUT_DIR" -maxdepth 1 -type f -name 'p1h_*.bin' -print0)

  while IFS= read -r -d '' f; do
    local base idx
    base="$(basename "$f")"
    [[ -f "$f.complete" ]] || continue
    [[ "$base" == "$latest_base" ]] && continue
    if [[ "$base" =~ restart_i([0-9]+)_ ]]; then
      idx="${BASH_REMATCH[1]}"
      # Wait for the solver to publish this as latest before compressing it.
      (( 10#$idx < 10#$latest_index )) || continue
    fi
    gzip -1 -n -f "$f"
    local latest_now
    latest_now="$(cat "$OUT_DIR/latest_checkpoint.txt" 2>/dev/null || true)"
    if [[ "$latest_now" == "$base" ]]; then
      printf '%s.gz\n' "$base" >"$OUT_DIR/latest_checkpoint.txt.part"
      mv "$OUT_DIR/latest_checkpoint.txt.part" "$OUT_DIR/latest_checkpoint.txt"
    fi
  done < <(find "$OUT_DIR" -maxdepth 1 -type f -name 'p1h_restart*.dump' -print0)

  # Keep only the newest and immediately previous complete restart. Field
  # snapshots remain available at the heating, hold, heat-off, and solidified events.
  mapfile -t checkpoints < <(find "$OUT_DIR" -maxdepth 1 -type f -name 'p1h_restart*.dump.gz' -printf '%f\n' | sort -rV)
  for ((n=MAX_SAVED_CHECKPOINTS; n<${#checkpoints[@]}; n++)); do
    local old="${checkpoints[$n]}"
    [[ -f "$OUT_DIR/${old%.gz}.complete" ]] || continue
    rm -f -- "$OUT_DIR/$old" "$OUT_DIR/${old%.gz}.state" "$OUT_DIR/${old%.gz}.complete"
    printf 'PRUNED_OLD_CHECKPOINT %s\n' "$old" >>"$OUT_DIR/checkpoint_cleanup.log"
  done
}

compress_all_after_exit() {
  while IFS= read -r -d '' f; do gzip -1 -n -f "$f"; done \
    < <(find "$OUT_DIR" -maxdepth 1 -type f \( -name 'p1h_*.bin' -o -name 'p1h_*.dump' \) -print0)
  local latest
  latest="$(find "$OUT_DIR" -maxdepth 1 -type f -name 'p1h_restart*.dump.gz' -printf '%f\n' | sort -rV | head -1)"
  if [[ -n "$latest" ]]; then
    printf '%s\n' "$latest" >"$OUT_DIR/latest_checkpoint.txt.part"
    mv "$OUT_DIR/latest_checkpoint.txt.part" "$OUT_DIR/latest_checkpoint.txt"
  fi
  compress_completed_artifacts
}

export OMP_NUM_THREADS="$OMP_THREADS" OMP_DYNAMIC=FALSE
export P1H_MARANGONI=on P1H_MARANGONI_MULTIPLIER=3
export P1H_ARC_PRESSURE_ENABLED=1 P1H_ARC_PRESSURE_PEAK_PA="$PRESSURE_PA"
export P1H_ARC_PRESSURE_PROFILE=periodic_gaussian P1H_ARC_PRESSURE_RADIUS_NM="$PRESSURE_RADIUS_NM"
export P1H_ARC_PRESSURE_RAMP_PS=0.02 P1H_SOURCE_RAMP_PS=0.02

(cd "$OUT_DIR" && "$EXE" "$INPUT_DIR" negative 1e11 "$MAX_TIME_PS" >run.stdout.log 2>run.stderr.log) &
solver_pid=$!
compress_completed_artifacts
while state="$(ps -o stat= -p "$solver_pid" 2>/dev/null || true)" &&
      [[ -n "$state" && "$state" != Z* ]]; do
  compress_completed_artifacts
  sleep 3
done
set +e
wait "$solver_pid"
solver_status=$?
set -e
printf '%s\n' "$solver_status" >"$OUT_DIR/exit_code.txt"
compress_all_after_exit
if rg -q 'DONE_SOLIDIFIED' "$OUT_DIR/run.stderr.log"; then
  touch "$OUT_DIR/.case_finished"
  echo "P1I completed through solidification: $OUT_DIR"
elif rg -q 'MAX_TIME_STOP' "$OUT_DIR/run.stderr.log"; then
  touch "$OUT_DIR/.case_time_guard"
  echo "P1I reached the configured time guard; latest valid checkpoint retained: $OUT_DIR" >&2
  exit 3
else
  touch "$OUT_DIR/.case_failed"
  echo "P1I stopped abnormally (exit=$solver_status); inspect $OUT_DIR/run.stderr.log" >&2
  exit "${solver_status:-1}"
fi
