# Round7 POD/SNAP v12g validation status — 2026-11-18

## Decision

Do not start heating. The ordinary 12-frame validation errors are below the exploratory force-RMSE target of 3 eV/Å, but the current potentials fail transfer/stability checks on the 458,248-atom bicontinuous model.

## Metrics and what they mean

- POD v12g: energy RMSE 0.0795 eV/atom; force RMSE 0.6046 eV/Å on the ordinary validation set.
- SNAP v12g: energy RMSE 0.0570 eV/atom; force RMSE 0.6273 eV/Å on the ordinary validation set.
- On the separate held-out Si01 center (ID 104359), the vector force errors are 7.192 eV/Å for POD and 54.500 eV/Å for SNAP. These are one-atom local errors, not global RMSE. Si01 remains excluded from training.
- A max force of 500 eV/Å is only an emergency ceiling. It is not a pass criterion for relaxation or MD; the present structure shows severe collapse at much lower Fmax.

## Full-model checks

- Model: 458,248 atoms, corrected 200 × 200 × 300 Å geometry.
- POD v12g static screen: Fmax 18.222 eV/Å, F99 6.158 eV/Å; no force exceeded 500 eV/Å.
- Short cold screen: 10 steps at dt=1 fs, 300 K velocity initialization, bottom skeleton fixed and bath layer at 300 K, no laser/TTM heat. All 458,248 atoms retained, but mobile temperature rose to 977.9 K. Potential energy fell by 37,333 eV while kinetic energy rose by 37,011 eV; the temperature increase is from rapid structural energy release, not applied heating. There are 53 unique Ti–C pairs below 1.8 Å (106 atoms); shortest pair is IDs 62864/62871 at 1.68265 Å. This is not cold stability.
- Full-model FIRE: 100 iterations reached the iteration cap, not the minimizer convergence criterion. All atoms retained; final Fmax 17.844 eV/Å and F99 4.621 eV/Å. The potential energy fell by 0.3943 eV/atom, with 4,418 nearest-neighbor pairs below 1.8 Å, including 4,277 Ti–C pairs; the shortest contact is 1.37239 Å. The minimized structure is unusable for MD.

## Local DFT follow-up

- The initial 36-atom crop around Ti 62864–C 62871 converged with GPAW PBE/PW350/Gamma in 34 SCF iterations. v12g center-force errors are 2.625 eV/Å (Ti) and 10.171 eV/Å (C).
- DFT labels for that pair in the 10 fs cold snapshot and for an independent short Ti–C pair (IDs 324028/324034) are in progress. GPAW MPI failed due launcher/library incompatibility and a parallel matrix-distribution error; the active retries use single-process BLAS threads. No incomplete result will enter training.
- If labels converge, prepare a separate v12h candidate from the 128 baseline train / 12 ordinary validation split, the established local-force set, and the new targeted Ti–C rows. Keep Si01 center ID 104359 held out.

## Heating gate

No heating, cooling, or solidification run has been started. A candidate must first improve the local validation/stability behavior and pass a short cold full-model check without continuing Ti–C collapse or unphysical temperature growth.
