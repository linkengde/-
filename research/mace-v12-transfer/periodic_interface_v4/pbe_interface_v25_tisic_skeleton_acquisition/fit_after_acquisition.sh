#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd -- "$(dirname -- "$0")" && pwd)"
REPO="/workspace/jeitschko-main-run"
PYTHON="/workspace/.venvs/mace/bin/python"
STATUS="$HERE/sequence_status.json"
LOCK="/workspace/.v25-skeleton-fit.lock"
LOG="$HERE/calculations/postprocess_supervisor.log"
mkdir -p "$(dirname -- "$LOCK")" "$HERE/calculations"
exec 9>"$LOCK"
flock -n 9 || { echo 'Another v25 skeleton postprocessor is already active.' >&2; exit 2; }
exec > >(tee -a "$LOG") 2>&1

while true; do
  state="$("$PYTHON" -c 'import json,sys; print(json.load(open(sys.argv[1])).get("status", "missing"))' "$STATUS")"
  case "$state" in
    all_skeleton_archives_verified) break ;;
    running) sleep 30 ;;
    paused_by_user|failed_stop_no_retry|failed*)
      echo "Acquisition ended in state $state; refusing to fit partial data."
      exit 3
      ;;
    *)
      echo "Unexpected acquisition state: $state"
      exit 4
      ;;
  esac
done

echo "All DFT archives verified; assembling and fitting the Ti-Si-C skeleton candidate."
cd "$REPO"
"$PYTHON" "$HERE/prepare_dataset.py"
"$PYTHON" "$HERE/preflight_train_skeleton.py"
bash "$HERE/run_train_skeleton.sh"
echo "Skeleton fit verification finished: $HERE/run/results/skeleton_fit_verification.json"
