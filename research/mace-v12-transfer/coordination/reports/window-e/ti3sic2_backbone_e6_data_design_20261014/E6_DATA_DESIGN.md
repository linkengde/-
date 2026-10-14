# Stage E6 data design and evidence review

## Scope and integrity

E2, E3, E4 and E5 checksum manifests were checked successfully. The E5 package has 868 files including its manifest; all 867 manifest entries pass SHA-256. E5 is archived on `cloud-e-b1-dft` at commit `195436e25d8bf8e4a524ce035b2c35c519f0bdc9`. The original TSCBOP file and E5 archived source copy both hash to `729273a6488d75756ed7927432c77c7347a1b115cecbaefae4981bafcec80fd7`.

No input data, DFT label, structure, or potential was edited in this stage. The E5 force Jacobian is reused as-is. `role_freeze.json` remains TRAIN/DEV/TEST empty.

## What the current force data cover

E5 contains a `1062 × 8` finite-difference force response matrix built from 354 atom-force vectors (six geometries: five cases from one 48-atom bulk parent family and one 114-atom TSC-only slab derived from an Ag interface). Each parameter column perturbs one tied unordered-pair/scalar block; all six shared triplet rows for that block were changed together. The response is ` [F(1.01p)-F(0.99p)]/0.02`, an approximate derivative with respect to `ln(p)`.

The six active columns have algebraic rank 6. Their unit-column-normalized condition number is `1104.743`, so the matrix is badly conditioned even before energy, stress, elastic, or independent-family tests. E5 found strong correlations including `B_Ti-C` / `lambda1_Ti-C` (`+0.999881`), `A_Ti-Si` / `B_Ti-Si` (`-0.998186`) and `A_Ti-C` / `lambda1_Ti-C` (`-0.998078`). These values support reducing simultaneous degrees of freedom; they do not identify a unique physical cause or fitted value.

## Why E5 cannot constrain Si–C A/B

Every tested geometry has all periodic minimum Si–C distances outside the Tersoff taper: `3.4982–3.5965 Å` versus the outer cutoff `R+D=3.211505 Å`. Thus the Si–C radial block contributes exactly zero to the tested forces, giving zero columns for `A_Si-C` and `B_Si-C`. This is a geometric activation gap, not evidence that those parameters are physically unimportant or erroneous.

The smallest *activation test* is a physically relevant, symmetry-audited local environment containing at least one Si–C pair inside the nonzero taper, with converged DFT forces on the same frozen geometry. Do not obtain it by an arbitrary compression solely to activate the code. A clean surface reconstruction, step/edge, or point-defect motif is a candidate only if it falls within the intended solid-backbone domain and DFT confirms that geometry. First screen candidate geometries with the unchanged potential's distance/taper logic. If no intended operating structure activates Si–C, freeze the Si–C blocks at their source values and exclude them from the fitted scope. If both Si–C A and B are to be fitted, one active snapshot is not enough evidence: require distinct Si–C distance/coordination environments and an independent held-out family, then recompute rank and conditioning.

## New DFT versus extending `run 0`

For an existing geometry that already has verified DFT forces, additional LAMMPS `run 0` variants can extend its force-response matrix without new DFT. E5 has already done this for the eight requested blocks on the six existing geometries. `run 0` cannot create reference forces, stresses, energy differences, or independent validation data. The missing Si response, volume/shear stress, clean-surface labels, and independent DEV/TEST evidence therefore require new DFT labels unless an already archived, method-matched result is found and audited. A candidate Si–C-active geometry can be screened for pair activation before spending DFT time; it still needs DFT labels before it can constrain a fit objective.

## E5 effective-rank threshold convention

Recalculation from the archived `PARAMETER_RESPONSE_MATRIX.csv`, after dropping the two zero columns and normalizing each active column to unit norm, gives singular-value ratios `1, 0.27017, 0.05047, 0.01742, 0.003777, 0.0009052`: effective rank 5 at relative threshold `1e-3`, and 4 at `1e-2`. `JACOBIAN_RANK.json` reports raw, unnormalized singular values, for which a `1e-2` threshold gives rank 3; `SVD_ANALYSIS.md` reports the normalized result. This is a scaling-convention difference, not a disagreement in algebraic rank 6 or in the normalized condition number. Future reports must state the column scaling before quoting thresholded effective rank.

## Main design consequence

The force-only E5 matrix supports a cautious two-parameter candidate for a later, separately authorized exploration. It does not support an eight-parameter fit, a Si–C fit, a claim of independent validation, or a production potential. Elastic, energy, surface, and split-family data must be added and audited before any fit gate can be reconsidered.
