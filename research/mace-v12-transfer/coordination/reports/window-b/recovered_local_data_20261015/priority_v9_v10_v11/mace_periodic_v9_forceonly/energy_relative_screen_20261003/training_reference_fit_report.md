# v9 training-only energy reference diagnostic

Only the v9 training energy labels were used to fit reference-energy correction coefficients. Validation/test and external PW-PBE frames were not used in the fit. No model artifact was changed.

- Energy-labeled training frames: 9
- Species composition matrix rank: 3/4
- Species terms (eV/atom): C=+0.02953140, Si=+0.12352657, Ti=-0.04570446, Ag=+0.02562796
- Training energy RMSE before reference correction: 17.374 meV/atom
- Training energy RMSE after the training-only reference fit: 3.394 meV/atom

| Split | Energy frames | Raw RMSE (meV/atom) | After train-only ref fit (meV/atom) |
|---|---:|---:|---:|
| train | 9 | 17.374 | 3.394 |
| valid | 2 | 19.744 | 7.283 |
| test | 2 | 24.926 | 7.483 |

This diagnostic tests whether the observed absolute-energy bias is mostly a reference term. If the design-matrix rank is below four, individual element references cannot be uniquely inferred from the current training set; add energy-labeled structures with independent compositions before using per-element corrections. A small corrected RMSE on train alone is not validation; valid/test and untouched external structures control that claim.
