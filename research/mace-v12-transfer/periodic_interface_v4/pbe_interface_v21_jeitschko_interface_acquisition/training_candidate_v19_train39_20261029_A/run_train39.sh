#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd -- "$(dirname -- "$0")" && pwd)"
DATA="$HERE/dataset"
OUT="$HERE/run"
FOUNDATION="${MACE_FOUNDATION_MODEL:-/workspace/-/research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v12_interface_energy/foundation_models/mace-mp-0b3-medium.model}"
PYTHON="${MACE_PYTHON:-/workspace/.venvs/mace/bin/python}"
TRAIN="${MACE_RUN_TRAIN:-/workspace/.venvs/mace/bin/mace_run_train}"
THREADS="${MACE_NUM_THREADS:-4}"
[[ -x "$PYTHON" && -x "$TRAIN" ]] || { echo 'MACE runtime missing.' >&2; exit 2; }
[[ -s "$FOUNDATION" ]] || { echo 'Foundation model missing.' >&2; exit 2; }
[[ "$THREADS" =~ ^[1-4]$ ]] || { echo 'MACE_NUM_THREADS must be 1-4.' >&2; exit 2; }
"$PYTHON" "$HERE/preflight_train39.py"
mkdir -p "$OUT"/{logs,models,checkpoints,results}
export XDG_CACHE_HOME="$OUT/cache" MPLCONFIGDIR="$OUT/mpl-cache" FC_CACHEDIR="$OUT/cache/fontconfig"
export OMP_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS=1
printf 'CPU training threads: %s\n' "$THREADS"
"$TRAIN" --name MACE_periodic_v19_residual_response --work_dir "$OUT" \
 --log_dir "$OUT/logs" --model_dir "$OUT/models" --checkpoints_dir "$OUT/checkpoints" --results_dir "$OUT/results" \
 --train_file "$DATA/train.extxyz" --valid_file "$DATA/valid.extxyz" \
 --energy_key REF_energy --forces_key REF_forces --foundation_model "$FOUNDATION" --E0s estimated \
 --seed 45 --device cpu --batch_size 1 --max_num_epochs 80 --lr 0.0001 --weight_decay 5e-7 \
 --energy_weight 100 --forces_weight 1000 --save_all_checkpoints --keep_checkpoints 2>&1 | tee "$OUT/train.stdout"
[[ -s "$OUT/checkpoints/MACE_periodic_v19_residual_response_run-45.model" ]] || { echo 'Expected selected checkpoint absent; preserve outputs.' >&2; exit 3; }
[[ -s "$OUT/models/MACE_periodic_v19_residual_response.model" ]] || { echo 'Expected final model absent; preserve outputs.' >&2; exit 3; }
grep -Eq 'Loaded Stage one model from epoch 79 for evaluation' "$OUT/train.stdout" || { echo 'Expected epoch-79 export marker absent; preserve outputs.' >&2; exit 3; }
"$PYTHON" - "$OUT" <<'PY'
import datetime, json, sys
from pathlib import Path
(Path(sys.argv[1])/'launcher_completion.json').write_text(json.dumps({
    'training_command_exit_code':0,
    'export_epoch79_marker_verified':True,
    'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'status':'TRAIN39_EXPERIMENT_TRAINING_COMPLETE_PENDING_VERIFICATION'
},indent=2)+'\n')
PY
