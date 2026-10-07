#!/usr/bin/env bash
set -euo pipefail
[[ $# == 4 && "$1" == --repo-root && "$3" == --output-dir ]] || { echo 'usage: run_v17_training.sh --repo-root REPO --output-dir NEW_RUN_CHILD' >&2; exit 2; }
REPO="$(realpath -- "$2")"
OUT="$(realpath -m -- "$4")"
ENTRY="$REPO/research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v17_interface_energy"
[[ "$OUT" == "$ENTRY/"* && "$OUT" != "$ENTRY" ]] || { echo 'Output must be a new child of the V17 run directory.' >&2; exit 2; }
[[ ! -e "$OUT" || ( -d "$OUT" && -z "$(find "$OUT" -mindepth 1 -print -quit)" ) ]] || { echo 'Refusing nonempty output/existing model.' >&2; exit 2; }
PYTHON="${MACE_PYTHON:-/workspace/.venvs/mace-v12/bin/python}"
TRAIN="${MACE_RUN_TRAIN:-/workspace/.venvs/mace-v12/bin/mace_run_train}"
[[ -x "$PYTHON" && -x "$TRAIN" ]] || { echo 'MACE runtime missing.' >&2; exit 2; }
"$PYTHON" -B "$ENTRY/preflight_v17.py" --repo-root "$REPO"
FOUNDATION="$REPO/research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v12_interface_energy/foundation_models/mace-mp-0b3-medium.model"
DATA="$ENTRY/data"
mkdir -p "$OUT"/{logs,models,checkpoints,results}
cd "$OUT"
export XDG_CACHE_HOME="$OUT/cache" MPLCONFIGDIR="$OUT/mpl-cache" FC_CACHEDIR="$OUT/cache/fontconfig"
THREADS="${MACE_NUM_THREADS:-4}"
[[ "$THREADS" =~ ^[1-4]$ ]] || { echo 'MACE_NUM_THREADS must be 1-4 on this five-core instance.' >&2; exit 2; }
export OMP_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS=1
echo "CPU training threads: $THREADS"
"$TRAIN" --name MACE_periodic_v17_provisional_screening --work_dir "$OUT" \
 --log_dir "$OUT/logs" --model_dir "$OUT/models" --checkpoints_dir "$OUT/checkpoints" --results_dir "$OUT/results" \
 --train_file "$DATA/train.extxyz" --valid_file "$DATA/valid.extxyz" \
 --energy_key REF_energy --forces_key REF_forces --foundation_model "$FOUNDATION" --E0s estimated \
 --seed 45 --device cpu --batch_size 1 --max_num_epochs 80 --lr 0.0001 --weight_decay 5e-7 \
 --energy_weight 100 --forces_weight 1000 2>&1 | tee "$OUT/train.stdout"
[[ -s "$OUT/checkpoints/MACE_periodic_v17_provisional_screening_run-45.model" ]] || { echo 'Expected selected checkpoint absent; preserve outputs.' >&2; exit 3; }
[[ -s "$OUT/models/MACE_periodic_v17_provisional_screening.model" ]] || { echo 'Expected final model absent; preserve outputs.' >&2; exit 3; }
echo 'Training process finished. Inspect all 80 epochs and final checkpoint before writing a frozen selection record. V17 remains provisional screening; no blind labels or test_file were used.'
