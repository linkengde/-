#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(git -C "$HERE" rev-parse --show-toplevel)"
LMP_BIN="${LMP_BIN:-lmp}"
OVITO_PYTHON="${OVITO_PYTHON:-python}"
TRAJ="$(realpath "${1:?usage: run_analysis.sh trajectory.dump [reheat-screen.out] [output-dir]}")"
THERMO_SCREEN="${2:-$(dirname "$TRAJ")/screen.reheat_q6_ptm_msd.out}"
OUT="${3:-$HERE/analysis_$(date -u +%Y%m%dT%H%M%SZ)}"
mkdir -p "$OUT"
"$OVITO_PYTHON" "$HERE/analysis/analyze_ptm.py" "$TRAJ" "$OUT/ptm_per_frame.csv"
"$LMP_BIN" \
  -var initial_data "$HERE/data/ag_equilibrated_300K.data" \
  -var ag_eam_file "$REPO_ROOT/models/ag_ti3sic2_qsgs_1to1_finer/potentials/Ag_u3.eam" \
  -var input_dump "$TRAJ" \
  -in "$HERE/inputs/in.reanalyze_q6.lmp" \
  -log "$OUT/log.reanalyze_q6.lammps" \
  -screen "$OUT/screen.reanalyze_q6.out"
"$OVITO_PYTHON" "$HERE/analysis/assemble_summary.py" \
  "$OUT/ptm_per_frame.csv" "$OUT/screen.reanalyze_q6.out" \
  "$THERMO_SCREEN" "$OUT/reheat_q6_ptm_msd_summary.csv"
