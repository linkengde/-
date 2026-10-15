#!/usr/bin/env bash
set -euo pipefail
cd /workspace/.local/lammps-ttm/run_1ps_diagnostic_20261001/potential_screening/periodic_interface_v4/mace_periodic_v9_forceonly
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1
export XDG_CACHE_HOME=/tmp/mace-v9-cache MPLCONFIGDIR=/tmp/mpl-v9-cache
mkdir -p "$XDG_CACHE_HOME" "$MPLCONFIGDIR"
/workspace/.local/venvs/mace/bin/mace_run_train \
  --name Ag_Ti_Si_C_MACE_periodic_v9_forceonly \
  --work_dir "$PWD" --log_dir "$PWD/logs" --model_dir "$PWD/models" \
  --checkpoints_dir "$PWD/checkpoints" --results_dir "$PWD/results" \
  --train_file "$PWD/data/train.extxyz" --valid_file "$PWD/data/valid.extxyz" \
  --test_file "$PWD/data/test.extxyz" --energy_key REF_energy --forces_key REF_forces \
  --foundation_model /workspace/.local/cache/mace/mace-mp-0b3-medium.model --E0s estimated \
  --seed 42 --device cpu --batch_size 1 --max_num_epochs 80 \
  --lr 0.0001 --weight_decay 5e-7 --energy_weight 1 --forces_weight 1000 \
  2>&1 | tee "$PWD/train.stdout"
