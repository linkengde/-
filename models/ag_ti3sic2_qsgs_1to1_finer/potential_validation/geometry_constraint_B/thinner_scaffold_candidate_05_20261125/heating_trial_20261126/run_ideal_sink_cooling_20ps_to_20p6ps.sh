#!/usr/bin/env bash
# 从20 ps加热检查点停热，进行0.6 ps底部300 K理想热沉冷却。
set -euo pipefail
TRIAL="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AG_TTM_CHECKPOINT_DIR="$TRIAL/continued_19p5_to_20p0ps" \
AG_TTM_COOL_OUTPUT_DIR="$TRIAL/cooling_from20p0_to20p6ps_reproduction" \
AG_TTM_COOL_STEPS=1200 \
AG_TTM_COOL_COUPLING_SUBSTEPS=5 \
AG_TTM_BOTTOM_SINK_Z_MAX_A=50 \
AG_TTM_BOTTOM_SINK_TEMPERATURE_K=300 \
AG_TTM_COOL_BOX_ZHI_A=500 \
/workspace/.venvs/qsgs/bin/python "$TRIAL/cool_with_ideal_bottom_sink_candidate05.py"
