# Stage E6 fit-readiness decision

## Decision

- `START_TRIAL_FIT=NO`
- `PRODUCTION_READY=NO`

## Why the fit gate remains closed

1. The current `1062×8` force Jacobian is only locally rank 6; two Si–C columns are identically zero and the six active directions have normalized condition number `1104.743`.
2. Best-conditioned current subsets of size 2/3/4 have conditions `3.288 / 14.156 / 46.718`; none has been tested on the proposed new data or on a fit objective containing energies, stress, and independent families.
3. Existing geometries provide no Si–C response; Si–C cannot be estimated from zero columns. Freeze those blocks unless target-relevant, DFT-labelled configurations activate them.
4. Bulk stress, volume/shear DFT response, Si displacement labels and independent clean-surface references are missing.
5. There are no frozen independent DEV/TEST source families; current role freeze is empty.
6. The 48-atom CIF-derived bulk structure is not a verified PBE equilibrium cell; the 114-atom slab is interface-derived and highly forced.
7. Numerical parameter bounds and prior constraints cannot be inferred from ±1% probes. E3 reports external paper-parameter comparison while local access/independent verification remained unverified. Do not turn source parameter values into fit bounds without checking units and author constraints.

## Smallest path to reconsideration

1. Keep source potential SHA and all DFT sources immutable; retain E5 tied-block mapping.
2. Audit/reuse existing DFT labels and decide whether Si–C interactions are in the intended solid operating domain.
3. Obtain the four volume/shear points, two Si displacement points, and two clean-surface reference/response families, unless equivalent method-matched DFT files are found and independently checked.
4. Freeze family-level TRAIN/DEV/TEST in a C-approved manifest; no common parent lineage may cross roles.
5. Recompute response Jacobians on the selected, actual DFT-labelled training geometries; report column scaling, rank, singular values, condition, and family/element response. Keep any candidate set to 2 blocks initially unless evidence supports more.
6. Only then submit a separate fit protocol for approval, with a frozen objective, family-balanced weights, source-backed bounds, regularization, resource estimate, and a test-sealed evaluation plan.

This stage did not optimize parameters, change a potential, assign data roles, generate a production file, or start DFT/MD/TTM.
