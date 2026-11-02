#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd -- "$(dirname -- "$0")" && pwd)"
DATA="$HERE/dataset"
OUT="$HERE/run"
FOUNDATION="${MACE_FOUNDATION_MODEL:-/workspace/-/research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v12_interface_energy/foundation_models/mace-mp-0b3-medium.model}"
PYTHON="${MACE_PYTHON:-/workspace/.venvs/mace/bin/python}"
TRAIN="${MACE_RUN_TRAIN:-/workspace/.venvs/mace/bin/mace_run_train}"
THREADS="${MACE_NUM_THREADS:-4}"
[[ -x "$PYTHON" && -x "$TRAIN" && -s "$FOUNDATION" ]] || { echo 'MACE runtime or foundation model missing.' >&2; exit 2; }
[[ "$THREADS" =~ ^[1-4]$ ]] || { echo 'MACE_NUM_THREADS must be 1-4.' >&2; exit 2; }
"$PYTHON" "$HERE/preflight_train_skeleton.py"
RESUME_ARG=()
TEE_APPEND=()
if [[ -d "$OUT/checkpoints" ]] && compgen -G "$OUT/checkpoints/MACE_TiSiC_skeleton_v25_run-45_epoch-*.pt" >/dev/null; then
  [[ ! -s "$OUT/launcher_completion.json" ]] || { echo 'skeleton fit already completed; refusing to rerun.' >&2; exit 2; }
  RESUME_ARG=(--restart_latest); TEE_APPEND=(-a)
else
  if [[ -d "$OUT" ]] && find "$OUT" -mindepth 1 -print -quit | rg -q .; then
    echo 'Partial fit outputs have no resumable epoch checkpoint; preserving them.' >&2; exit 2
  fi
fi
mkdir -p "$OUT"/{logs,models,checkpoints,results}
export XDG_CACHE_HOME="$OUT/cache" MPLCONFIGDIR="$OUT/mpl-cache" FC_CACHEDIR="$OUT/cache/fontconfig"
export OMP_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS=1
"$TRAIN" --name MACE_TiSiC_skeleton_v25 --work_dir "$OUT" \
 --log_dir "$OUT/logs" --model_dir "$OUT/models" --checkpoints_dir "$OUT/checkpoints" --results_dir "$OUT/results" \
 --train_file "$DATA/train.extxyz" --valid_file "$DATA/valid.extxyz" \
 --energy_key REF_energy --forces_key REF_forces --foundation_model "$FOUNDATION" --E0s estimated \
 --seed 45 --device cpu --batch_size 1 --max_num_epochs 80 --lr 0.0001 --weight_decay 5e-7 \
 --energy_weight 100 --forces_weight 1000 --save_all_checkpoints --keep_checkpoints "${RESUME_ARG[@]}" 2>&1 | tee "${TEE_APPEND[@]}" "$OUT/train.stdout"
[[ -s "$OUT/checkpoints/MACE_TiSiC_skeleton_v25_run-45.model" ]] || { echo 'Selected checkpoint absent.' >&2; exit 3; }
[[ -s "$OUT/models/MACE_TiSiC_skeleton_v25.model" ]] || { echo 'Exported model absent.' >&2; exit 3; }
"$PYTHON" - "$OUT" <<'PY'
import datetime,json,sys
from pathlib import Path
(Path(sys.argv[1])/'launcher_completion.json').write_text(json.dumps({
 'training_command_exit_code':0,
 'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
 'status':'SKELETON_TRAINING_COMPLETE_PENDING_VERIFICATION'
},indent=2)+'\n')
PY
"$PYTHON" "$HERE/verify_skeleton_fit.py" > "$OUT/results/skeleton_fit_verification.stdout"
