# V19 train39 Jeitschko interface experiment

## Result

The experiment-only refit completed all 80 epochs with exit code 0. The epoch-79 exported model is tensor-identical to the epoch-79 checkpoint. The frozen development screen failed: aggregate force-vector RMSE is 0.20958 eV/Å, above the 0.05 eV/Å per-interface limit. The official V19 train36/dev3 files and old V19 model remain unchanged; no test or sealed data were opened.

| Interface | Force-vector RMSE (eV/Å) | Limit | Separation-force error (eV/Å) | Limit | Energy absolute error (meV/atom) |
|---|---:|---:|---:|---:|---:|
| AgC | 0.196795 | 0.05 | 0.082977 | 0.10 | 6.140 |
| AgSi | 0.231677 | 0.05 | 0.081472 | 0.10 | 38.765 |
| AgTi | 0.198404 | 0.05 | 0.155916 | 0.10 | 2.307 |

All three force-vector gates fail; AgTi also fails the separation-force gate. The old V19 force-vector RMSEs were 2.300, 2.750 and 2.505 eV/Å for AgC, AgSi and AgTi respectively, so this experiment substantially improved the screen but did not qualify the model. Species-resolved force-vector RMSE is 0.0592 eV/Å for Ag, 0.2608 for C, 0.0331 for Si and 0.2211 for Ti; the largest residuals are on C and Ti.

## Inputs and recipe

The training file has 39 frames: the byte-identical frozen V19 train36 plus three independently calculated Ag/Ti₃SiC₂ interface labels. The validation file is the exact frozen 3-frame dev3. The added structures share the Jeitschko–Nowotny 1967 parent family, while dev3 belongs to the distinct COD 9009647/Kisi 1998 family. Their labels map GPAW native extrapolated energy to `REF_energy`, sigma 0.10 eV forces to `REF_forces`, and retain free energy separately. The interface geometries are fixed single-point acquisition structures, not relaxed equilibrium configurations.

The fixed V19 recipe was seed 45, CPU with four threads, batch size 1, 80 epochs, learning rate 1e-4, weight decay 5e-7, energy weight 100, force weight 1000, and the same foundation model. The verifier evaluates only the frozen dev3 and confirms the epoch-79 export matches the selected checkpoint. The evaluation script records the recipe, input hashes, all 80 epoch validation metrics, per-species and per-interface errors, and the pass/fail thresholds in `run/results/train39_verification.json`.

## Files and use

- `dataset/`: frozen train39/dev3 inputs, manifest and checksums.
- `run/results/train39_verification.json`: complete verification and development-screen result.
- `run/results/train39_epoch_metrics.csv`: validation metrics for epochs 0–79.
- `run/logs/` and `run/train.stdout`: training logs and epoch-79 export record.
- `run/checkpoints/MACE_periodic_v19_residual_response_run-45_epoch-79.pt`: final checkpoint.
- `run/models/MACE_periodic_v19_residual_response.model`: exported epoch-79 model.
- `run/checkpoints/MACE_periodic_v19_residual_response_run-45.model`: checkpoint model export used for identity verification.

All 80 epoch checkpoints were present and verified before archiving. Intermediate checkpoints and the derived compiled model were removed to keep the repository package compact; the complete epoch coverage remains in the verified logs and result JSON. The retained model is an experiment artifact only. It has not passed the development gates and must not be treated as a qualified production or molecular-dynamics model. Slab k-point, smearing, and vacuum convergence are also not fully established.

## Status

`TRAIN39_EXPERIMENT_DEV3_SCREEN_FAIL`. This experiment is separate from the official dataset and model. The result does not resolve the V19 development-set failure. Residuals are largest on C and Ti, suggesting that the next acquisition should add audited force-rich Ti/C bulk configurations before another interface refit; any follow-up should remain a separately versioned experiment with the same frozen dev3 screen.
