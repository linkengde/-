#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
KDIR="$HERE/k1x6x3"
DFT_PID="144943"
GPAW_PY="/workspace/.venvs/gpaw-diagnosis/bin/python"
ACE_PY="/workspace/.venvs/ace/bin/python"
PIPELINE_LOG="$HERE/fit_pipeline_after_dft.log"

exec > >(tee -a "$PIPELINE_LOG") 2>&1
echo "[$(date -u +%FT%TZ)] Waiting for GPAW parent PID $DFT_PID to finish."
while :; do
  state="$(ps -p "$DFT_PID" -o stat= 2>/dev/null || true)"
  state="${state%% *}"
  if [[ -z "$state" || "$state" == Z* ]]; then
    break
  fi
  sleep 30
done

if [[ ! -f "$KDIR/result_resume.json" ]]; then
  echo "Missing DFT result_resume.json; preserving outputs and stopping."
  exit 10
fi
status="$("$GPAW_PY" -c 'import json,sys; print(json.load(open(sys.argv[1])).get("status", "MISSING"))' "$KDIR/result_resume.json")"
if [[ "$status" != "SCF_CONVERGED" ]]; then
  echo "Dense-mesh calculation status: $status; skip label-dependent refits."
  exit 11
fi

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
echo "[$(date -u +%FT%TZ)] Comparing k-point labels."
"$GPAW_PY" "$KDIR/compare_ti_gap_kpoint.py" | tee "$KDIR/comparison.stdout"
echo "[$(date -u +%FT%TZ)] Preparing candidate datasets."
"$GPAW_PY" "$HERE/prepare_fit_candidate.py" | tee "$HERE/fit_candidate.stdout"
echo "[$(date -u +%FT%TZ)] Running ACE refits."
"$ACE_PY" "$HERE/fit_ace_kconv.py" > "$HERE/ace_kconv_refits.stdout" 2>&1
echo "[$(date -u +%FT%TZ)] Running MTP force-weight sweep."
bash "$HERE/run_mtp_kconv.sh" > "$HERE/mtp_kconv_refits.stdout" 2>&1
echo "[$(date -u +%FT%TZ)] ACE and MTP refits completed."
