#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
FOUNDATION_MODEL="${FOUNDATION_MODEL:-$ROOT/../mace_periodic_v12_interface_energy/foundation_models/mace-mp-0b3-medium.model}"
[[ -f "$FOUNDATION_MODEL" ]] || { echo "Foundation model not found: $FOUNDATION_MODEL" >&2; exit 2; }
if [[ -f "$ROOT/checkpoints/MACE_periodic_v13_interface_energy_run-45.model" && "${ALLOW_RETRAIN:-0}" != 1 ]]; then
  echo "v13 is already trained. Read its assessment and continue new DFT checks; ALLOW_RETRAIN=1 explicitly permits retraining." >&2
  exit 2
fi
MACE_TRAIN="${MACE_RUN_TRAIN:-mace_run_train}"
if ! command -v "$MACE_TRAIN" >/dev/null 2>&1; then
  MACE_TRAIN="${MACE_VENV_BIN:-/workspace/.venvs/mace-v12/bin}/mace_run_train"
fi
[[ -x "$MACE_TRAIN" ]] || command -v "$MACE_TRAIN" >/dev/null 2>&1 || { echo "MACE training command unavailable; set MACE_RUN_TRAIN or MACE_VENV_BIN." >&2; exit 2; }
cd "$ROOT"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-/tmp/mace-v13-cache}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/mpl-v13-cache}"
export FC_CACHEDIR="${FC_CACHEDIR:-$XDG_CACHE_HOME/fontconfig}"
mkdir -p logs models checkpoints results "$XDG_CACHE_HOME" "$MPLCONFIGDIR" "$FC_CACHEDIR"
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
"$MACE_TRAIN" \
  --name MACE_periodic_v13_interface_energy \
  --work_dir "$ROOT" --log_dir "$ROOT/logs" --model_dir "$ROOT/models" \
  --checkpoints_dir "$ROOT/checkpoints" --results_dir "$ROOT/results" \
  --train_file "$ROOT/data/train.extxyz" --valid_file "$ROOT/data/valid.extxyz" \
  --test_file "$ROOT/data/test.extxyz" --energy_key REF_energy --forces_key REF_forces \
  --foundation_model "$FOUNDATION_MODEL" --E0s estimated \
  --seed 45 --device cpu --batch_size 1 --max_num_epochs 80 \
  --lr 0.0001 --weight_decay 5e-7 --energy_weight 100 --forces_weight 1000 \
  2>&1 | tee "$ROOT/train.stdout"
