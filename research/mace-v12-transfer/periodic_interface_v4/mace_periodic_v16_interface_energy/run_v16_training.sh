#!/usr/bin/env bash
set -euo pipefail
[[ $# == 4 && "$1" == --repo-root && "$3" == --output-dir ]] || { echo 'usage: run_v15_training.sh --repo-root REPO --output-dir NEW_REPORT_CHILD' >&2; exit 2; }
REPO="$(realpath -- "$2")"
OUT="$(realpath -m -- "$4")"
ENTRY="$REPO/research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v16_interface_energy"
[[ "$OUT" == "$ENTRY/"* && "$OUT" != "$ENTRY" ]] || { echo 'Output must be a new child of the V15 run directory.' >&2; exit 2; }
[[ ! -e "$OUT" || ( -d "$OUT" && -z "$(find "$OUT" -mindepth 1 -print -quit)" ) ]] || { echo 'Refusing nonempty output/existing model.' >&2; exit 2; }
PYTHON="${MACE_PYTHON:-/workspace/.venvs/mace-v12/bin/python}"
TRAIN="${MACE_RUN_TRAIN:-/workspace/.venvs/mace-v12/bin/mace_run_train}"
[[ -x "$PYTHON" && -x "$TRAIN" ]] || { echo 'MACE runtime missing; A must prepare it.' >&2; exit 2; }
# A runtime only: confirms all three blind archives; no blind values select the checkpoint.
"$PYTHON" -B "$ENTRY/preflight.py" --repo-root "$REPO" --check-blind-archives
FOUNDATION="$REPO/research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v12_interface_energy/foundation_models/mace-mp-0b3-medium.model"
DATA="$REPO/research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v16_interface_energy/data"
mkdir -p "$OUT"/{logs,models,checkpoints,results}
cd "$OUT"
export XDG_CACHE_HOME="$OUT/cache" MPLCONFIGDIR="$OUT/mpl-cache" FC_CACHEDIR="$OUT/cache/fontconfig"
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
"$TRAIN" --name MACE_periodic_v16_interface_energy --work_dir "$OUT" \
 --log_dir "$OUT/logs" --model_dir "$OUT/models" --checkpoints_dir "$OUT/checkpoints" --results_dir "$OUT/results" \
 --train_file "$DATA/train.extxyz" --valid_file "$DATA/valid.extxyz" --test_file "$DATA/test.extxyz" \
 --energy_key REF_energy --forces_key REF_forces --foundation_model "$FOUNDATION" --E0s estimated \
 --seed 45 --device cpu --batch_size 1 --max_num_epochs 80 --lr 0.0001 --weight_decay 5e-7 \
 --energy_weight 100 --forces_weight 1000 2>&1 | tee "$OUT/train.stdout"
# Success is process exit only; A must inspect 80-epoch and selected-checkpoint evidence.
[[ -f "$OUT/checkpoints/MACE_periodic_v16_interface_energy_run-45.model" ]] || { echo 'Expected checkpoint absent; preserve outputs.' >&2; exit 3; }
echo 'Training process finished; A must review epoch/checkpoint evidence and freeze selection before evaluation. No model PASS is assigned.'
