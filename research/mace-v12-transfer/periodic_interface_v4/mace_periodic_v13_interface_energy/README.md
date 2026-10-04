# MACE v13 Ag/Ti/Si/C interface update

The v13 training set preserves all v12 training frames and appends the three verified PW-PBE distance labels. The v12 validation and test files are byte-identical. Training completed 80 epochs on CPU; the selected checkpoint is epoch 72.

The current validation result is **UNDETERMINED**. On the three frozen test frames, v13 is close to v12 overall (force-vector RMSE 0.05716 vs 0.05821 eV/A), although the C-gap frame is worse (0.09288 vs 0.08492 eV/A). On the two held-out Ag-C/Ag-Si midpoint structures, the combined force RMSE falls from 0.09224 to 0.02412 eV/A. Ag-C midpoint RMSE falls from 0.11659 to 0.02546 eV/A; its separating-force error falls from 0.279 to 0.110 eV/A. Ag-Si also improves, with a remaining 0.143 eV/A separating-force error.

Those midpoint checks share source geometries with the distance scans added to v13 training, so they test interpolation within related interface motifs. They do not establish transfer to new registries, thermal interface configurations, liquid Ag, or the production candidate. Keep long MD and TTM on hold pending independent interface geometries. The trainer's 1697.4 meV/atom train energy RMSE includes four force-only frames with no energy labels as zero-valued targets; an audit reproduces that value and finds 7.45 meV/atom unweighted energy RMSE on the 22 actually labelled frames (5.41 meV/atom when atom weighted). See `results/v13_training_frame_audit.json`.

See `results/v9_v10_v11_v12_v13_comparison.json` for per-structure metrics and model hashes, `results/v13_validation_assessment.json` for the staged assessment, and `data/dataset_manifest.json` for dataset hashes and split provenance.
