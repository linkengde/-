#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "$0")" && pwd)"
REPO_ROOT="$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)"
MODEL_DIR="$REPO_ROOT/models/ag_ti3sic2_qsgs_1to1_finer"
PHASE_DIR="$MODEL_DIR/phase_field_qsgs"
PYTHONPATH="$PHASE_DIR${PYTHONPATH:+:$PYTHONPATH}" python "$SCRIPT_DIR/analyze_existing_run.py" \
  --model-dir "$MODEL_DIR" --diag-dir "$SCRIPT_DIR/source_data" --out-dir "$SCRIPT_DIR"
python "$SCRIPT_DIR/run_sensitivity.py"
