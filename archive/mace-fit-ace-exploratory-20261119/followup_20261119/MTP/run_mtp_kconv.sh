#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="${HERE%/dft_kpoint_recheck_20261119}"
DATA="$HERE/fit_candidate"
OUT="$HERE/mtp_kconv_refits"
MLP="/workspace/potential-comparison/software/mlip-2-build/mlp"
PYTHON="/workspace/.venvs/gpaw-diagnosis/bin/python"
SCORE="$BASE/score_mlip_cfg.py"
INIT="$BASE/MTP10_AgCSiTi_v19v25_v24_untrained.mtp"
mkdir -p "$OUT"

for force_weight in 0.05 0.10 0.20; do
  tag="F${force_weight/./}"
  model="$OUT/MTP10_TiGap1x6x3_${tag}.mtp"
  "$MLP" train "$INIT" "$DATA/train.cfg" \
    --energy-weight=1 --force-weight="$force_weight" --max-iter=1000 \
    --valid-cfgs="$DATA/dev.cfg" --trained-pot-name="$model" \
    > "$OUT/train_${tag}.log" 2>&1
  "$MLP" calc-efs "$model" "$DATA/train.cfg" "$OUT/train_${tag}_pred.cfg" \
    > "$OUT/calc_train_${tag}.log" 2>&1
  "$MLP" calc-efs "$model" "$DATA/dev.cfg" "$OUT/dev_${tag}_pred.cfg" \
    > "$OUT/calc_dev_${tag}.log" 2>&1
  "$PYTHON" "$SCORE" "$DATA/train.cfg" "$OUT/train_${tag}_pred.cfg" "$OUT/train_${tag}_metrics.json"
  "$PYTHON" "$SCORE" "$DATA/dev.cfg" "$OUT/dev_${tag}_pred.cfg" "$OUT/dev_${tag}_metrics.json"
done
