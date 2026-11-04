#!/usr/bin/env bash
set -euo pipefail
: "${BASILISK_SRC:?Set BASILISK_SRC to the pinned Basilisk src directory}"
export BASILISK="$BASILISK_SRC"
export OPENGLIBS=-lfb_tiny
"$BASILISK_SRC/qcc" -disable-dimensions -O2 -Wall -autolink p1c_vof_dynamic.c -o p1c_vof_dynamic -lm
"$BASILISK_SRC/qcc" -disable-dimensions -O2 -Wall -autolink p1c_axisymmetric_force_balance.c -o p1c_axisymmetric_force_balance -lm
run_case() {
  local name="$1"; shift
  timeout 300 ./p1c_vof_dynamic "$@" > "$name.stdout.log" 2> "$name.stderr.log"
}
run_case combined_amr_r8 8 0 1 1 1.0
run_case combined_amr_r16 16 0 1 1 1.0
run_case combined_amr_r32 32 0 1 1 1.0
run_case combined_amr_r16_dt05 16 0 1 0.5 1.0
run_case normal_only_r16 16 1 1 1 1.0
run_case tangent_only_r16 16 2 1 1 1.0
run_case combined_fixed_r8 8 0 0 1 1.0
run_case balanced_tangent_r8 8 0 1 1 0.9792
run_case balanced_tangent_r16 16 0 1 1 0.9792
run_case balanced_tangent_r32 32 0 1 1 0.9792
for n in 8 16 32; do
  if [[ "$n" == 8 ]]; then
    ./p1c_axisymmetric_force_balance "$n" > p1c_axisymmetric_force_balance.csv
  else
    ./p1c_axisymmetric_force_balance "$n" | tail -n 1 >> p1c_axisymmetric_force_balance.csv
  fi
done
