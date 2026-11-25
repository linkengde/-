#!/usr/bin/env bash
# 从3.0 ps双检查点连续加热至20.0 ps，每0.5 ps保存一组原子/TTM检查点。
set -euo pipefail
TRIAL="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="/workspace/.venvs/qsgs/bin/python"
SCRIPT="$TRIAL/continue_heat_candidate05.py"
PREVIOUS="continued_2p5_to_3p0ps"
for START_TENTHS in $(seq 30 5 195); do
    END_TENTHS=$((START_TENTHS + 5))
    START_WHOLE=$((START_TENTHS / 10)); START_FRAC=$((START_TENTHS % 10))
    END_WHOLE=$((END_TENTHS / 10)); END_FRAC=$((END_TENTHS % 10))
    START_LABEL="${START_WHOLE}p${START_FRAC}"; END_LABEL="${END_WHOLE}p${END_FRAC}"
    NEXT="continued_${START_LABEL}_to_${END_LABEL}ps"
    if [[ -e "$TRIAL/$NEXT/continuation_summary.json" ]]; then
        echo "SKIP completed segment $NEXT"
        PREVIOUS="$NEXT"
        continue
    fi
    mkdir -p "$TRIAL/$NEXT"
    echo "START ${START_LABEL} -> ${END_LABEL} ps; checkpoint=$PREVIOUS"
    AG_TTM_CHECKPOINT_DIR="$TRIAL/$PREVIOUS" \
      AG_TTM_CONTINUE_OUTPUT_DIR="$TRIAL/$NEXT" \
      AG_TTM_CONTINUE_STEPS=1000 \
      "$PYTHON" "$SCRIPT" >"$TRIAL/$NEXT/continuation.log" 2>&1
    python3 - "$TRIAL/$NEXT/continuation_summary.json" "$END_WHOLE.$END_FRAC" <<'PY'
import json, sys
from pathlib import Path
summary = json.loads(Path(sys.argv[1]).read_text())
expected = float(sys.argv[2])
if summary.get("status") != "COMPLETED" or abs(summary.get("final_time_ps", -1)-expected)>1e-10:
    raise SystemExit(f"segment did not complete at {expected} ps: {summary}")
print("DONE", f"t={expected:.1f} ps", f"Fmax={summary['mobile_ag_force_max_seen_eV_A']:.4f} eV/A",
      f"ledger={summary['energy_ledger_residual_eV']:.6g} eV", flush=True)
PY
    test -s "$TRIAL/$NEXT/lammps.restart"
    test -s "$TRIAL/$NEXT/ttm_state.npz"
    test -s "$TRIAL/$NEXT/CHECKPOINT.json"
    PREVIOUS="$NEXT"
done
echo "Heating reached 20.0 ps. Final checkpoint: $TRIAL/$PREVIOUS"
