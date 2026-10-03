#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
FOUNDATION_MODEL="${FOUNDATION_MODEL:-$ROOT/foundation_models/mace-mp-0b3-medium.model}"
[[ -f "$FOUNDATION_MODEL" ]] || { echo "Foundation model not found: $FOUNDATION_MODEL" >&2; exit 2; }
cd "$ROOT"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-/tmp/mace-v12-cache}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/mpl-v12-cache}"
export FC_CACHEDIR="${FC_CACHEDIR:-$XDG_CACHE_HOME/fontconfig}"
mkdir -p logs models checkpoints results "$XDG_CACHE_HOME" "$MPLCONFIGDIR" "$FC_CACHEDIR"
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
mace_run_train \
  --name MACE_periodic_v12_interface_energy \
  --work_dir "$ROOT" --log_dir "$ROOT/logs" --model_dir "$ROOT/models" \
  --checkpoints_dir "$ROOT/checkpoints" --results_dir "$ROOT/results" \
  --train_file "$ROOT/data/train.extxyz" --valid_file "$ROOT/data/valid.extxyz" \
  --test_file "$ROOT/data/test.extxyz" --energy_key REF_energy --forces_key REF_forces \
  --foundation_model "$FOUNDATION_MODEL" --E0s estimated \
  --seed 45 --device cpu --batch_size 1 --max_num_epochs 80 \
  --lr 0.0001 --weight_decay 5e-7 --energy_weight 100 --forces_weight 1000 \
  2>&1 | tee "$ROOT/train.stdout"
