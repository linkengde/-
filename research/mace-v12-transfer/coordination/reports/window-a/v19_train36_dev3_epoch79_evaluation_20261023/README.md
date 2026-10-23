# V19 train36/dev3 fixed run and epoch-79 evaluation

## Result

The registered CPU fit completed normally with exit code 0. The log contains epochs 0–79, all 80 epoch checkpoints were present during verification, and the exported model is tensor-identical to the epoch-79 checkpoint across 99 state entries. The checkpoint and model hashes are recorded in `training_verification.json`.

The candidate **fails the development force and separation-force gates**. Do not use it for production, MD, or impact simulations. No test or sealed labels were opened.

| Development frame | Energy absolute error (meV/atom) | Force-vector RMSE (eV/Å) | Separation-force absolute error (eV/Å) | Result |
|---|---:|---:|---:|---|
| AgC | 672.275 | 2.299662 | 0.105592 | Fail both force gates |
| AgSi | 651.006 | 2.749850 | 1.554593 | Fail both force gates |
| AgTi | 667.965 | 2.505232 | 0.652634 | Fail both force gates |

Thresholds were 0.05 eV/Å per-frame force-vector RMSE and 0.10 eV/Å separation-force absolute error. AgC misses the latter narrowly. Pooled per-element force-vector RMSE is Ag 0.3794, Si 0.7121, C 2.4578, and Ti 3.0973 eV/Å; Ti and C dominate the residuals.

## Run and artifact verification

The fixed recipe used seed 45, CPU, four threads, batch size 1, 80 epochs, learning rate 1e-4, weight decay 5e-7, energy weight 100, and force weight 1000. It loaded the fixed foundation model and estimated elemental reference energies from train36. The verification script confirmed the epoch log/checkpoint coverage and exact epoch-79 export identity. `training_epoch_metrics.csv` contains all 80 validation rows.

The exported model and epoch-79 checkpoint are retained under the run directory in `mace_periodic_v19_residual_response/runs/v19_train36_dev3_20261023_A/`. Intermediate checkpoints 0–78 and transient caches are not part of the published bundle; the run log, complete epoch metrics, final checkpoint, and exported model are retained.

## Train/dev diagnostic

`train_dev_diagnostic.json` compares the foundation model and epoch-79 model using only approved train36/dev3. Train36 fit improves to 20.95 meV/atom energy RMSE and 0.11721 eV/Å force-component RMSE (0.20301 eV/Å force-vector RMSE). On dev3, the foundation model has 0.38679 eV/Å force-vector RMSE; epoch 79 rises to 2.52496 eV/Å. This is a large generalization gap, not evidence that the development labels are wrong.

The input audit shows a concrete domain and method mismatch to resolve before another fit:

- All 36 training frames use fully periodic `T T T` cells. Metadata records Gamma-point PW-PBE 500 eV for all; 27 frames have a generic Gamma single-point method string and 9 explicitly record Gamma with Fermi smearing 0.1 eV. Most cells are roughly 28 Å cubic.
- All 3 development frames are the same formula, `C16Ag4Si8Ti24`, absent from train36. They use `T T F` slabs with approximately 6.115 × 5.296 × 34 Å cells, a 6×6×1 mesh, sigma 0.10 eV, and xy dipole correction.
- The train composition matrix has rank 4, but its singular values span 152.759 to 0.397 (condition ratio about 385). The dev composition is outside the exact train composition list. The trainer estimated element reference corrections from this limited composition set; its log reports corrections of −16.709 eV for C, +9.876 eV for Si, +6.428 eV for Ti, and +4.748 eV for Ag. This makes composition extrapolation a plausible contributor to the nearly uniform dev energy offset, but does not prove it is the only cause.
- The inherited train rows carry an energy-reference caveat in the dataset manifest. Newer training labels and dev labels preserve native extrapolated energy and free energy separately, but the old rows' original energy-key lineage should be reconciled before using energy metrics to tune the next model.

These differences confound source-geometry generalization with periodic-boundary, k-point, smearing, and energy-reference differences. The three dev structures also share one Ti₃SiC₂ crystallographic source family, so they are not three independent source-family tests. Slab mesh/smearing/vacuum convergence is unresolved.

## Next data work

1. Audit the original source energy keys for every train36 row and document native extrapolated energy versus free energy, units, k-points, smearing, cell, and PBC. Preserve the existing train36/dev3 files while correcting any proven inconsistency in a versioned dataset.
2. Reuse only archive-verified Ag and Ti₃SiC₂ pure-phase/backbone labels whose settings and geometries are appropriate. Add force-rich Ti/C backbone configurations if coverage is missing; remote Ti/C sites dominate the current dev error.
3. Find a second, source-verified Ti₃SiC₂ parent CIF family. Build AgC, AgSi, and AgTi interface candidates with recorded structure/source hashes and separated data roles. Keep any structures used for fitting out of the held-out source-family development group.
4. Resolve or explicitly freeze slab numerical settings, then label only registered, hash-verified structures. Have B independently audit the new data entry before A approves a versioned train/dev split and runs the fixed preflight/training recipe again.

No new DFT or training was started as part of this evaluation and diagnosis. The B audit's shared-source limitation remains in force; test/sealed data remain unopened.

## Reproduction files

- `training_verification.json`: all-epoch and epoch-79 identity checks plus dev3 metrics.
- `training_epoch_metrics.csv`: per-epoch validation curve, epochs 0–79.
- `train_dev_diagnostic.json`: foundation-versus-epoch79 diagnostics on train36/dev3 only.
- `training_fit_detail.txt`: the trainer's complete structured fit/evaluation output.
- `training_stdout.log` and `training_metrics.log`: original compact run logs.
- `SHA256SUMS.txt`: hashes for this report, scripts, logs, final checkpoint/model, and approved data manifest.
- Reproduction scripts are in `research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v19_residual_response/`.
