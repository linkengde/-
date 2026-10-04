# MACE v12 AgSi/AgC energy additions

The four AgSi/AgC PW-PBE labels converged, passed archive checks, and were used to build the 23/2/3 train/validation/test split. MACE v12 completed 80 epochs; the selected checkpoint is from epoch 76.

The current validation status is **UNDETERMINED**. V12 improves energy errors over v11 on the three-frame frozen test and two external structures, while external force results are mixed. The MACE training summary also reports a 1800.8 meV/atom energy RMSE that does not match the separate frame-by-frame audit (41.0 meV/atom). Resolve this mismatch and add more independent interface checks before using v12 on candidate structures, long MD, or TTM.

See `results/v12_validation_assessment.json` for the stage-by-stage status and evidence, `results/v9_v10_v11_v12_candidate_comparison.json` for per-structure metrics and model hashes, and `results/v12_training_frame_audit.json` for the training-set audit.

Training used the v11 energy-focused objective (`energy_weight=100`, `forces_weight=1000`) and the bundled MACE-MP-0b3-medium model.
