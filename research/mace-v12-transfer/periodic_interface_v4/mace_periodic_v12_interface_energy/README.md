# MACE v12 AgSi/AgC energy additions

This model is staged but must not be trained until the four v12 PW-PBE single points have converged and their archive manifest passes review. It replaces the matching LCAO-only AgSi and AgC force-scan frames with same-geometry PW-PBE energy+force labels, preserves all current validation/test structures, and keeps the AgSi/AgC midpoint DFT structures external.

Training objective matches the successful v11 energy-focused trial (`energy_weight=100`, `forces_weight=1000`). Do not run long MD or TTM from v12 without passing the same frozen and external checks.
