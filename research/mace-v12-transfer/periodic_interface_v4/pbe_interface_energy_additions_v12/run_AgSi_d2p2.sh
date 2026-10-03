#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${GPAW_PYTHON:-python}"
RUN_ROOT="${GPAW_RUN_ROOT:-/tmp/pbe_interface_v12_runs}"
OUT="$RUN_ROOT/AgSi_d2p2"
[[ ! -e "$OUT" ]] || { echo "Refusing to overwrite existing run directory: $OUT" >&2; exit 2; }
"$PYTHON" "$HERE/run_contact_pw_reference.py" \
  AgSi_d2p2 "$HERE/inputs/AgSi_d2p2.extxyz" "$OUT"
