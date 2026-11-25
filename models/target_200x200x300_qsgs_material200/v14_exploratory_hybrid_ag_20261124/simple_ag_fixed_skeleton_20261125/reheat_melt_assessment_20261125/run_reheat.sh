#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(git -C "$HERE" rev-parse --show-toplevel)"
LMP_BIN="${LMP_BIN:-lmp}"
OUT="${1:-$HERE/reproduced_$(date -u +%Y%m%dT%H%M%SZ)}"
mkdir -p "$OUT"
"$LMP_BIN" \
  -var initial_data "$HERE/data/ag_equilibrated_300K.data" \
  -var ag_eam_file "$REPO_ROOT/models/ag_ti3sic2_qsgs_1to1_finer/potentials/Ag_u3.eam" \
  -var output_dump "$OUT/ag_reheat_q6_ptm_msd.dump" \
  -var final_dump "$OUT/ag_reheat_q6_ptm_msd_final.dump" \
  -in "$HERE/inputs/in.reheat.lmp" \
  -log "$OUT/log.reheat_q6_ptm_msd.lammps" \
  -screen "$OUT/screen.reheat_q6_ptm_msd.out"
