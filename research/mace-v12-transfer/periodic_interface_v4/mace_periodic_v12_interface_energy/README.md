# MACE v12 AgSi/AgC energy additions

The four AgSi/AgC PW-PBE labels converged, passed archive checks, and were used to build the 23/2/3 train/validation/test split. MACE v12 completed 80 epochs; the selected checkpoint is from epoch 76.

The current validation status is **UNDETERMINED**. Three added independent fixed-distance PW-PBE labels (Ag-C 2.30 A, Ag-Si 2.60 A, and Ag-Ti 2.60 A) converged and passed hash, geometry, and label checks. On three frozen test frames and five external fixed geometries, v12 substantially improves energy errors over v11: external energy MAE is 9.95 vs 127.46 meV/atom. External aggregate force-vector RMSE is effectively unchanged (0.08553 vs 0.08557 eV/A), with mixed contact-level results: Ag-C forces worsen, Ag-Ti improves, and Ag-Si vector RMSE improves while its separating force is underpredicted. The new scan points are distance perturbations of existing midpoint motifs, so more independent interface geometries are needed before v12 can support candidate structures, long MD, or TTM.

The MACE table's 1800.8 meV/atom training energy RMSE is explained by four force-only frames that have no energy labels but are included as zero-energy targets in its unmasked summary metric; the selected checkpoint's RMSE on the 19 actually energy-labelled training frames is 41.0 meV/atom. The next model iteration should add the three newly verified DFT labels to a separate v13 dataset, preserve the v12 model and frozen splits, and target further Ag-C/Ag-Si force checks.

See `results/v12_validation_assessment.json` for the stage-by-stage status and evidence, `results/v9_v10_v11_v12_candidate_comparison.json` for per-structure metrics and model hashes, and `results/v12_training_frame_audit.json` for the training-set audit.

Training used the v11 energy-focused objective (`energy_weight=100`, `forces_weight=1000`) and the bundled MACE-MP-0b3-medium model.
