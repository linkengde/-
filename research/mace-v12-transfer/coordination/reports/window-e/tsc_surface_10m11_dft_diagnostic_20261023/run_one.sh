#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../../../../../" && pwd)"
source "$HERE/runtime.sh"
TASK_ID="${1:?pass an exact structure_id}"
RESUME_FLAG="${2:-}"
case "$TASK_ID" in
  TSC69_GEOMHOLD_PYR_10m11_COD9009647_TITI_REF|TSC69_GEOMHOLD_PYR_10m11_COD9009647_TITI_RESP_TIZ020) ;;
  *) echo "Unknown task ID" >&2; exit 2 ;;
esac
CONFIG="$HERE/DFT_DIAGNOSTIC_CONFIG.json"
LOCK="$HERE/EXECUTION_LOCK.json"
INPUT="$REPO/research/mace-v12-transfer/coordination/reports/window-e/tsc_surface_geometry_dev_candidate_10m11_20261023/inputs/$TASK_ID.extxyz"
RUN_DIR="$HERE/results/$TASK_ID"
mkdir -p "$RUN_DIR"
CONFIG_SHA="$(sha256sum "$CONFIG" | awk '{print $1}')"
if [[ "$RESUME_FLAG" == "--resume" ]]; then
  RESUME=(--resume)
elif [[ -z "$RESUME_FLAG" ]]; then
  RESUME=()
else
  echo "Only optional flag is --resume" >&2; exit 2
fi
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
exec mpiexec -n 4 "$GPAW_BIN" python "$HERE/run_fixed_geometry.py" \
  --structure-id "$TASK_ID" --input "$INPUT" --run-dir "$RUN_DIR" \
  --config "$CONFIG" --config-sha256 "$CONFIG_SHA" --lock "$LOCK" "${RESUME[@]}" \
  >"$RUN_DIR/stdout_$STAMP.log" 2>"$RUN_DIR/stderr_$STAMP.log"
