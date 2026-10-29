# V19 train43 Ti/C backbone augmentation experiment

This separate experiment appends four verified Ti₃SiC₂ internal-coordinate displacement labels to the completed experiment-only train39 data. Train39's largest force residuals on frozen dev3 were C and Ti, so these paired Ti(4f)/C(4f) displacement labels directly increase backbone force-response coverage. The frozen 3-frame dev3 and old V19 model remain unchanged; no test or sealed data are included.

## Added source labels

The four calculations are from Jeitschko and Nowotny's 1967 single-crystal source family, reconstructed from full-text Fourier coordinates with ideal full occupancy. Their source records are pinned to repository commit `847676691c89faa78750802d2217bea9d8702f85`. Each source result passed its archived result verifier and the per-result SHA256 inventory before inclusion. The user waived a separate owner-side B review for this experiment; the labels remain experiment-only and are not official training assignments.

| Displaced site | Displacement (Å) | Native energy (eV/cell) | Displaced-site force z (eV/Å) | SCF iterations |
|---|---:|---:|---:|---:|
| Ti(4f) | -0.02 | -97.40285644 | -0.090502 | 33 |
| Ti(4f) | +0.02 | -97.38402897 | -0.809206 | 32 |
| C(4f) | -0.02 | -97.37629875 | 1.187084 | 30 |
| C(4f) | +0.02 | -97.41150932 | 0.561830 | 30 |

The labels use PBE, a 500 eV plane-wave cutoff, 10×10×2 k points, and 0.10 eV Fermi smearing. Native extrapolated energy is mapped to `REF_energy`, GPAW forces to `REF_forces`, and finite-width free energy is retained separately. Paired free-energy finite-difference force checks agree with the directly calculated Ti and C baseline forces within 0.0004 eV/Å. Train39 also contains Gamma/TTT and 6×6×1/TTF labels, so the combined training set retains a k-point and boundary-condition mixture.

## Dataset and fixed fit

The candidate dataset contains 43 training frames (the byte-identical 39-frame train39 prefix plus these four bulk labels) and an exact byte-copy of frozen dev3. Its composition matrix has rank 4; the raw condition number decreases from 40.71 for train39 to 38.53 for train43. This is an experiment-only set, not a replacement for official V19 train36/dev3.

The planned fit uses the same V19 recipe and foundation model as train39: seed 45, CPU, four threads, batch size 1, 80 epochs, learning rate 1e-4, weight decay 5e-7, energy weight 100, force weight 1000, and `REF_energy`/`REF_forces`. The screen evaluates only the frozen dev3 and retains the existing per-interface force-vector RMSE ≤0.05 eV/Å and separation-force error ≤0.10 eV/Å gates.

## Fit result

The fixed 80-epoch fit completed successfully. All epochs and checkpoints were verified before archival, and the exported epoch-79 model is tensor-identical to its checkpoint. On the frozen dev3, aggregate force-vector RMSE fell from 0.20958 eV/Å for train39 to 0.06003 eV/Å for train43. This is a large improvement, but the 0.05 eV/Å per-interface force gate still fails for all three interfaces, and AgTi still fails its 0.10 eV/Å separation-force gate.

| Interface | Force-vector RMSE (eV/Å) | Train39 (eV/Å) | Separation-force error (eV/Å) | Train39 (eV/Å) | Energy absolute error (meV/atom) |
|---|---:|---:|---:|---:|---:|
| AgC | 0.054053 | 0.196795 | 0.035325 | 0.082977 | 1.355 |
| AgSi | 0.071126 | 0.231677 | 0.045248 | 0.081472 | 19.107 |
| AgTi | 0.053195 | 0.198404 | 0.158689 | 0.155916 | 6.279 |

The final epoch's validation energy RMSE is 11.64 meV/atom and force-component RMSE is 34.66 meV/Å. Species force-vector RMSE is 0.0620 eV/Å for Ag, 0.0697 for C, 0.0259 for Si, and 0.0609 for Ti. The complete epoch metrics and gate results are in `run/results/train43_verification.json` and `run/results/train43_epoch_metrics.csv`.

Status: `TRAIN43_EXPERIMENT_DEV3_SCREEN_FAIL`. All 80 epoch checkpoints were verified before archiving; only epoch 79 is retained with the model exports, logs, and validation record to keep the package compact. Official train36/dev3, train39, and the prior V19 model are unchanged. This candidate is not a qualified production or molecular-dynamics model. The bulk source assumes full occupancy despite the paper's discussion of possible nonstoichiometry. Source-family review is not a universal supercell-equivalence audit, and slab k-point, smearing, and vacuum convergence remain incomplete.
