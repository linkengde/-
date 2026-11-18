# Round7 v12f follow-up review (2026-11-18)

## Archive

The completed v12f artifacts and pre-existing v11 SI hotspot labels were archived to GitHub branch `codex/round7-v12f-pod-snap-artifacts-20261118`, commit `38e9f4b7097cb01d2acfe60970ed70e03c77099a`. The archive contains the v12f coefficients/configurations, input datasets, local DFT labels, evaluation outputs, run summaries, and SHA-256 manifest. The SI split remains SI00/SI02/SI03 training and SI01 development holdout. SI01 is not in the training set.

## SI01 holdout result

The center atom is Si ID 104359 in a 38-atom crop. Its GPAW reference force is `[1.9764, 1.8379, 12.1415] eV/Å`.

| Model | Predicted center force (eV/Å) | Vector error (eV/Å) | Direction cosine |
|---|---:|---:|---:|
| POD v12f | `[-7.7108, 6.1393, 4.3436]` | 13.1586 | 0.3641 |
| SNAP v12f | `[22.4870, -37.8937, -5.3358]` | 48.0077 | -0.1630 |

These are single-center vector-error norms, not dataset RMSEs. They exceed the project’s exploratory force-RMSE target of about 3 eV/Å and the force directions are poor. The all-atom component RMSEs on the local crop are POD 2.675 and SNAP 6.290 eV/Å; treat these as descriptive only because local crop edge forces are not the fit target.

For context, the ordinary 12-configuration validation set has aggregate force RMSE 0.613 eV/Å for POD and 0.652 eV/Å for SNAP (412 atoms total). Those overall values meet the exploratory target, but do not cover the distinctive SI01 compressed Si–Ag/Ti environment; SI01 therefore remains a separate local failure.

The SI01 center has a distinctive compressed mixed environment: one Ti at 1.909 Å and three Ag at 1.970, 2.066, and 2.130 Å. The three training Si centers do not have this same short Ti–Si plus Ag shell. A separate training candidate, Si ID 244669 (38 atoms), has a similar local shell and is extracted under `active_si_hotspot_expansion/`; DFT labeling is pending. SI01 remains held out.

The poor SI01 result is not only a crop-boundary artifact. In the full 347,464-atom static calculations, the POD prediction for ID 104359 was `[-7.91281, 4.00959, 4.72619] eV/Å`, with a DFT vector error of 12.5499 eV/Å and direction cosine 0.3929. SNAP predicted `[21.0504, -37.7574, -2.65663] eV/Å`, with a 46.3745 eV/Å error and direction cosine -0.1115. A run-0 check on the 38-atom full-cell crop exactly reproduced the recorded local-crop POD force, confirming that the crop evaluation is repeatable; it does not replace the full-environment result.

## Full-model geometry and stability findings

The current 347,464-atom `start_wrapped.data` has a 200 × 200 × 300 Å simulation box, but its atoms occupy only `z = -0.258 ... 150.53 Å`. There are no atoms above `z=200 Å`. The model builder explicitly sets a 200 × 200 × 150 Å material slab plus 150 Å vacuum. This differs from interpreting the full target as material through `z=200 Å`; confirm geometry scope before any heating run.

A separate corrected-height QSGS geometry has now been built at `models/target_200x200x300_qsgs_material200/`. It keeps the 200 × 200 × 300 Å box, fills material through z≈198.36 Å, and leaves z=200–300 Å as vacuum. It contains 458,248 atoms with Ag:Ti₃SiC₂ mass ratio 1.000003:1; both phases percolate across x/y and span z. Its minimum nearest-neighbor distance is 2.079 Å (Ti–C); minimum Ag–skeleton distance is 2.293 Å; no atoms have a nearest neighbor below 1.8 Å. A POD v12f `run 0` on this corrected geometry gives Fmax 22.3546 eV/Å (Ag 22.3546, C 11.4358, Si 13.4359, Ti 14.8887), with no atom above 500 eV/Å. This is a static force screen only; a full FIRE run is in progress and the new geometry has not passed cold stability or MD validation.

The initial model contains severe short contacts: the minimum is 1.441 Å (C–Ag, IDs 141572/311294), followed by 1.567 Å (Ag–Ag, IDs 294310/294311). 1,321 atoms have a nearest-neighbor distance below 1.8 Å. Full static POD forces reach 67.86 eV/Å.

The full 1 fs no-laser run lasted 0.01 ps: mobile-region temperature rose from 300 K to 2054 K, potential energy fell by about 69,051 eV, maximum force fell from 67.86 to 23.23 eV/Å, and maximum displacement reached 0.628 Å. This is rapid structural relaxation; it does not pass a 300 K cold-stability check.

Original FIRE logs stopped at step 60 after reporting 347,462 or 347,463 atoms. The padded run was not a clean control: it changed the z boundary, used four rather than five MPI ranks, and fixed only the bottom ceramic skeleton, while the original fixed every atom in z=0–20 Å.

An exact-settings replay (same start structure, potential, `p p m` boundary, fixed group, five MPI ranks, and FIRE settings; only the lost-atom policy, iteration limit, and diagnostic dump were changed) completed 62 iterations. Frames 60–62 each contained all 347,464 unique IDs, with finite coordinates. No atom crossed the nonperiodic z bounds; the few x/y coordinates outside [0,200] Å were at most 0.032 Å beyond a periodic boundary. Thus the reported loss was not reproduced as a physical atom escaping the cell, but the original `Lost atoms` error is still unexplained and not resolved.

The replay was not converged: FIRE stopped at its iteration limit. The force-norm maximum was 49.70 eV/Å at step 60 and 21.47 eV/Å at step 62; the energy fell from -2,061,013.31 to -2,241,271.96 eV (about -0.519 eV/atom). The original run’s step-60 loss error therefore occurred during severe structural relaxation, not after a completed minimization. The existing cold no-laser run also failed the 300 K stability check: mobile temperature rose to about 2,054 K in 0.01 ps as potential energy dropped sharply.

## Current decision

Do not start heating from this v12f model. The exact replay rules out an observed z-boundary escape in that run, but does not explain the original LAMMPS atom-count error, and FIRE remains unconverged. The independent Si candidate is now being labeled by GPAW; if it converges, audit the label, fit a separate v12g copy while keeping SI01 held out, then rerun full static/FIRE and cold-stability checks. The material-height discrepancy must also be resolved before claiming validation of the requested full geometry.
