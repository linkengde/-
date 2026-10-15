#!/usr/bin/env bash
set -euo pipefail
ROOT=/workspace/.local/lammps-ttm/run_1ps_diagnostic_20261001/potential_screening/periodic_interface_v4/mace_periodic_v10_energycal
cd "$ROOT"
mkdir -p logs models checkpoints results /tmp/mace-v10-cache /tmp/mpl-v10-cache
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
export XDG_CACHE_HOME=/tmp/mace-v10-cache MPLCONFIGDIR=/tmp/mpl-v10-cache
/workspace/.local/venvs/mace/bin/mace_run_train \
  --name MACE_periodic_v10_energycal \
  --work_dir "$ROOT" --log_dir "$ROOT/logs" --model_dir "$ROOT/models" \
  --checkpoints_dir "$ROOT/checkpoints" --results_dir "$ROOT/results" \
  --train_file "$ROOT/data/train.extxyz" --valid_file "$ROOT/data/valid.extxyz" \
  --test_file "$ROOT/data/test.extxyz" --energy_key REF_energy --forces_key REF_forces \
  --foundation_model /workspace/.local/cache/mace/mace-mp-0b3-medium.model --E0s estimated \
  --seed 43 --device cpu --batch_size 1 --max_num_epochs 80 \
  --lr 0.0001 --weight_decay 5e-7 --energy_weight 1 --forces_weight 1000 \
  2>&1 | tee "$ROOT/train.stdout"
