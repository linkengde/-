# Round7 v12k archive

This is an isolated POD/SNAP refit extending v12j with one converged Ti–C local force label. It retains the 128 training / 12 validation split and keeps Si01 center atom 104359 outside training.

- `candidate/v12k_ticpair_candidate/`: complete fitter input data, POD and SNAP coefficients, validation summaries, holdout result, and local center comparison.
- `labels/tic_pair5_v12j_cold_10fs_pw350_gamma/`: extracted crop, converged GPAW label, SCF output, and metadata.
- `full_model/`: exact LAMMPS scripts and logs for run-0 static POD/SNAP screens and POD 10 fs no-laser cold screen. The 458,248-atom final dump is compressed; the corresponding final `.data` file is omitted because positions/velocities are already represented by the dump and the archived starting model is present under `models/target_200x200x300_qsgs_material200/`.
- `ROUND7_STATUS_V12K_20261118.md`: interpretation and decision.

The 12-frame validation RMSE meets the exploratory 0.1 eV/atom and 3 eV/Å targets. This is not full-model MD approval: the POD cold screen rises from 300 K to about 707 K in 10 fs with no laser. Do not treat the static Fmax cutoff of 500 eV/Å as a stability criterion.
