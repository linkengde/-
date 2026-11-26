#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../../../../.." && pwd)"
PYTHON="${AGW_PYTHON:-/workspace/.venvs/ag-qsgs/bin/python}"
TMP="$(mktemp -d)"
OUT="${AGW_OUTPUT_DIR:-$HERE/reruns/run_$(date +%Y%m%d_%H%M%S)}"
mkdir -p "$OUT/snapshots"
trap 'rm -rf "$TMP"' EXIT

gzip -dc "$HERE/run_inputs/v15_relaxed.data.gz" > "$TMP/v15_relaxed.data"
gzip -dc "$HERE/run_inputs/v15_300K_after_0p2ps.dump.gz" > "$TMP/v15_300K_after_0p2ps.dump"
export LD_LIBRARY_PATH="/workspace/.venvs/ag-qsgs/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

"$PYTHON" "$REPO/research/ag_w_reference_pack/models/ag_w_qsgs_200A/ttm/agw_ttm_md_pilot.py" \
  --data "$TMP/v15_relaxed.data" \
  --audit "$HERE/v15_exploratory_manifest.json" \
  --exploratory \
  --restart-dump "$TMP/v15_300K_after_0p2ps.dump" \
  --restart-step 400 \
  --output "$OUT/v15_ttm_scale80_fwhm0p5ps_cooling5ps.json" \
  --snapshot-dir "$OUT/snapshots" \
  --snapshot-times-ps 0.5 1 2 4 7.5 \
  --grid-A 10 \
  --geometry-fractions-npz "$HERE/v15_qsgs_phase_fractions_5A_10A_20A.npz" \
  --dt-ps 0.0005 \
  --pulse-fwhm-ps 0.5 \
  --pulse-center-ps 0.5 \
  --cooling-ps 5 \
  --source-scale 80 \
  --interface-h-W-m2-K 1e9 \
  --initial-electron-temperature-from-lattice \
  --seed 20261124 \
  --omp-threads 1

echo "Result and snapshots retained under $OUT"
