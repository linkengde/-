# TSCBOP 48-atom P0 six-point publication and data-role audit

Date: 2026-10-17

## Decision and role status

This release publishes the six already completed fixed-geometry DFT points and the E-side audit. No DFT, LAMMPS run 0, fitting, MD, role reassignment, or potential modification was performed.

TRAIN=[], DEV=[], and TEST=[] remain unchanged. All six points belong to the single bulk_parent_48 family and remain DIAGNOSTIC_ONLY pending C label qualification and separate explicit use approval. C's review permits considering this family for a future TRAIN replay, but does not assign it to TRAIN or authorize new DFT or fitting.

C clarifies that frozen-potential diagnostic run 0 does not require TRAIN status, but such results cannot be claimed as independent validation. No new run 0 was requested or run in this publication task.

## C review binding

- C review commit: de24ee62e38b33bae87db0cfb4c6198e84cd4e95
- Review path: research/mace-v12-transfer/coordination/reports/window-c/tsc48_family_role_review_20261017/C_ROLE_REVIEW.md
- Versioned family manifest: research/mace-v12-transfer/coordination/reports/window-c/tsc48_family_role_review_20261017/FAMILY_USAGE_MANIFEST_v1.json
- Manifest SHA-256: dbfdca1be0579bfeb72fd4694cd5a42af084ccca0baaa9db4dc3f6c492eced05
- C manifest records TRAIN/DEV/TEST empty and role_change_authorized=false.
- Each published input.extxyz SHA matches its corresponding record in the C manifest.

## Existing result verification

Each of the six result directories passed its own SHA256SUMS.txt. Input hashes match the C manifest and RESULT.json; labelled.extxyz hashes match the output hash in RESULT.json. All six report converged SCF, 48 finite force vectors, finite stress, and the same fixed-geometry method: GPAW 26.7.0, ASE 3.29.0, PBE/PAW, PW500, Fermi-Dirac 0.10 eV, Monkhorst-Pack 6x6x2, MPI=4, OMP=1. Detailed values and per-file hashes are in P0_SIX_POINT_PUBLICATION_MANIFEST.json; the raw inputs, labelled outputs, logs, statuses and results are published at its listed paths.

The source batch summary is research/mace-v12-transfer/coordination/reports/window-e/dft_supplement_results_20261016/P0_FIXED_GEOMETRY_DFT_SUMMARY_20261016.md (SHA-256 d0fdf5304288579080c6ceb43bdae6cad0819ba783e0e7becce3ec8c01dbe3d9). These six descendants are one parent family, not six independent validation families. Label qualification remains pending C review.

C's “8 additional DFT points” is a conditional budget only, not an execution authorization. No additional points were run.

## Stop state

- Existing six-point results published; no recalculation.
- TRAIN/DEV/TEST unchanged.
- No DFT, run 0, fitting, MD or TTM in this release.
- Label qualification: pending C review.

## Launcher status note

The P0_VOL_m02 supervisor recorded exit_code=0. For the other five records, RUN_STATUS.json has a null exit_code because the batch supervisor did not flush it or the parent launcher exit code was unavailable. Those five are not represented as verified process exit codes. Their archived RESULT.json reports SUCCESS and SCF convergence; the 48-force/stress outputs, input/output SHA bindings, and per-run SHA256SUMS were verified. The result payload and checksum evidence support the recorded DFT label status, while the missing launcher exit codes remain a documented provenance limitation.
