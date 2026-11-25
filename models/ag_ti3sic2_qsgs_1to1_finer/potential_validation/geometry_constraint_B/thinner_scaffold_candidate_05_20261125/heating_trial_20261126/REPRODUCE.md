# Reproduce / continue candidate 05 heating

This is an exploratory ideal-support run, not a calibrated four-element contact potential. The Ti/Si/C atoms are fixed and Ag contacts them through a WCA repulsive wall.

The repo is expected at `/workspace/-`; Python/LAMMPS is the existing `/workspace/.venvs/qsgs` environment. The model root already contains `ag_ttm_reference.py`, `potentials/Ag_u3.eam`, the candidate geometry, `potential.geometry_wca_fixed_skeleton.inc`, and TTM parameter CSV.

Run the initial 0.5 ps segment from the model root:

```bash
cd /workspace/-/models/ag_ti3sic2_qsgs_1to1_finer
AG_TTM_TRIAL_STEPS=1000 AG_TTM_TRIAL_OUTPUT_DIR=/workspace/-/models/ag_ti3sic2_qsgs_1to1_finer/potential_validation/geometry_constraint_B/thinner_scaffold_candidate_05_20261125/heating_trial_20261126/ttm_initial_0p5ps /workspace/.venvs/qsgs/bin/python potential_validation/geometry_constraint_B/thinner_scaffold_candidate_05_20261125/heating_trial_20261126/ttm_heat_trial_candidate05.py
```

Each continuation segment is 1000 steps = 0.5 ps at 0.5 fs. Set `AG_TTM_CHECKPOINT_DIR` to the preceding segment folder and `AG_TTM_CONTINUE_OUTPUT_DIR` to a new output folder, then run `continue_heat_candidate05.py`. The 3.0 ps LAMMPS restart and TTM grid are paired in `continued_2p5_to_3p0ps/`; do not separate or overwrite them.

Recompute the trajectory analysis in OVITO's Python environment by placing/using the frame symlinks in `analysis_q6_ptm_msd_0_to_3ps/` and running `analyze_q6_ptm_msd.py`. Export per-particle values with `export_ovito_frame_3ps.py`, then render the saved arrays with `MPLCONFIGDIR=/tmp/mpl MPLBACKEND=Agg python3 render_3ps_software.py`. The image is rendered from OVITO-computed data with Matplotlib because the Qt offscreen Vulkan integration is unavailable on this host.
