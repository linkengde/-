# POD/SNAP v12g Si04 candidate

This directory contains a separate candidate refit. It does not overwrite v12f.

## Data and role split

- 128 baseline configurations remain in training.
- 14 local force configurations are in the training group, including the new Si hotspot centered on atom ID 244669.
- 12 ordinary validation configurations remain validation data.
- Si hotspot 01, center ID 104359, remains excluded from training and is reported separately as a development holdout.
- The new Si04 DFT calculation used GPAW PBE/PW350/Gamma with 0.10 eV smearing and converged in 50 SCF cycles. Its DFT center force is `[-1.3300, 0.6769, -6.3098] eV/Å`.

The rich local source frames are kept in `dataset/local_train/local_dft_targets.xyz`. Fitter copies retain only positions and named forces, in that order. The marked DFT target rows use the DFT values; surrounding rows use POD v10 force anchors and the configured local force weight. `local_fit_input_audit.json` confirms the column order and values.

The first POD attempt was stopped before being used because its custom-property order was unsafe for FitPOD's positional extxyz reader. Its partial log is retained as `fitpod_v12g_aborted_bad_property_order.log`. The fitted v12g model uses the corrected input copy and has a separate log.

## Current results

| Check | POD v12g | SNAP v12g |
|---|---:|---:|
| 12-frame ordinary validation force RMSE | 0.605 eV/Å | 0.627 eV/Å |
| Si01 held-out center vector error | 7.19 eV/Å | 54.50 eV/Å |
| Si01 held-out direction cosine | 0.816 | -0.169 |

For context, v12f's corresponding Si01 vector errors were 13.16 eV/Å for POD and 48.01 eV/Å for SNAP. Adding Si04 improves POD's Si01 local error but does not bring it below the exploratory 3 eV/Å force-error target. It does not improve SNAP on this holdout.

These local vector errors are single-center checks, not global RMSE. The ordinary validation RMSE is below 3 eV/Å, while the Si01 hotspot remains difficult. Neither v12g potential is cleared for heating from these results alone.

## Next checks

1. The full-model v12f FIRE and v12g static/cold screens are complete; see `ROUND7_CURRENT_STATUS_20261118.md` at the archive root.
2. The independent 458,248-atom v12g FIRE audit is in progress locally; its result will be added after completion.
3. Do not start heating until short contacts and rapid cold energy release are resolved.
