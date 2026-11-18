# Round7 v12l POD/SNAP candidate and cold-screen archive

This package extends v12k with a converged local Ti-C DFT force label, refits POD/SNAP, and records independent holdout plus 458,248-atom screens. It retains the 128/12 split and keeps Si01 center 104359 outside training.

- `candidate/v12l_ticpair_candidate/`: full fitter inputs/outputs, coefficients, validation, Si01 holdout, and local-force diagnostics.
- `labels/`: pair6 Ti-C candidate and converged GPAW label.
- `full_model/`: exact static and cold-run scripts/logs; full-size dumps/data are gzip-compressed. The starting 458,248-atom data file is reused from `models/target_200x200x300_qsgs_material200/Ag_Ti3SiC2_bicontinuous_200x200x300.data.gz` already in the parent archive.
- `ROUND7_STATUS_V12L_20261118.md`: full metrics and decision.

Both models meet the small 12-frame RMSE target, but neither passes the full-model cold screen. POD v12l crosses 500 eV/A at cumulative 37 fs by a local Ag-C-C collapse; SNAP v12l reaches 3330 K and creates 1989 Ti-C pairs shorter than 1.8 A in 10 fs. No heating/cooling/solidification run is approved.
