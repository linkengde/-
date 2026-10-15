# C review — Ti₃SiC₂ P0 six-point DFT input package

Review date: 2026-10-16 (Asia/Shanghai)

## Decision

```text
INPUT_INTEGRITY_PASS=YES
GEOMETRY_REVIEW_PASS=YES
ROLE_DECISION=DIAGNOSTIC_ONLY
TRAIN=NO
DEV=NO
TEST=NO
SCIENTIFIC_METHOD_REVIEW=PASS_WITH_EXPLICIT_LOCK_VALUES
EXECUTION_READY=NO
EXECUTION_AUTHORIZED=NO
```

This is a read-only C review. It does not authorize DFT and does not reclassify any historical label.

## Binding and evidence

- Reviewed E branch commit: `f3e09b7a3b595d21c1fd02cd72a6dda9cb74b8c1`.
- E package manifest SHA-256: `421499e6594eb090fdd7076ef6ae2ba6fc246a133283336a024609f6e7f0b29f`.
- E package `SHA256SUMS.txt`: all 16 listed files passed after fetching the commit.
- Source parent is 48 atoms, Ti24Si8C16, fully periodic; its file SHA-256 `157c1dab02a28e6140dd5191187ccd914fe558950e2ac6c136866bcc7805c3e4` matches the published baseline input at source commit `499ddf0870235594444abad780582ea4fcc89d61`.
- Existing baseline result `Ti3SiC2_baseline_k10x10x4_v19` is SCF-converged and has archived verification status `PASS`; its energy/force output geometry matches this parent. The associated `role_freeze.json` SHA-256 is `709d84b6f7205d3d437bae92dc1aaecd81278311f807b805e06c1b7d22acc27f`.

The six candidates preserve atom count, composition, atom ordering/IDs, and PBC. Independent transform checks reproduced the requested volume scaling, basal shear, and Si-ID-7 z displacements within the extxyz rounding tolerance. Their measured minimum periodic separations are 2.074–2.102 Å for volume variants, 2.085 Å for shear variants, and 2.088 Å for Si displacement variants. No geometry defect was found.

| Task ID | Transform | Structure SHA-256 |
|---|---|---|
| `TSC48_VOL_m02` | volume factor 0.98 | `1e8aeed059783a36e471873a3b3821369750677ae93ec5965d1d5760189ed683` |
| `TSC48_VOL_p02` | volume factor 1.02 | `a86ee4560676e1a9679729cd6c9721b4b4a875d41cc80862d6287a056ddb0db1` |
| `TSC48_BASAL_XY_SHEAR_m005` | `x′=x−0.005y` | `72be50f6e2fc3c7fcb610fadf520a7be3eadb8b14c2c7f27416da2fe7881ed7b` |
| `TSC48_BASAL_XY_SHEAR_p005` | `x′=x+0.005y` | `cd0e3d09e413d25a54dbaa38868106934b71e2d95c3c879c8346dbd0cff8f3a1` |
| `TSC48_SI_ID7_Z_m020` | Si ID 7, z −0.02 Å | `77ab3ae7dc59a41adea035003543fe19ae9c8deb8dd2e1413807ed25164e91f4` |
| `TSC48_SI_ID7_Z_p020` | Si ID 7, z +0.02 Å | `8222f1538719c859765f1ab48db483e4599d4f64ac2835f67dc5709ee860002d` |

## Data role decision

All six are children of the same 48-atom parent family as the existing baseline and Ti/C displacement controls. E2 freezes that entire family as `numerical_pure_phase_control_not_training_or_independent_validation`, with TRAIN/DEV/TEST empty. I retain that role unchanged: these six labels, if later authorized, are diagnostic-only response controls. They are not training labels, independent DEV, or TEST. The six variants cannot be split across roles. Any later request to train on this family requires a separate, explicit family-level role/split review; it must not silently alter the E2 freeze.

## Scientific method freeze for a future execution package

The proposed method is consistent with the verified archived 48-atom baseline. Its actual GPAW log reports `MonkhorstPack(size=[10,10,4])`, Monkhorst-Pack shift `[0,0,0]`, and PPCG `niter=2`. The final lock must state this exact convention; do not relabel the even mesh as Gamma-containing.

For a separately authorized fixed-geometry run, bind these settings to every input:

- GPAW 26.7.0, ASE 3.29.0, `gpaw-data` 1.2.1, PBE, PAW, plane-wave cutoff 500 eV.
- Fermi–Dirac width 0.10 eV; energy, density, and eigenstate convergence thresholds each `1e-5`; Pulay beta 0.05, `nmaxold=8`, weight 100; maximum 160 SCF iterations.
- k mesh `(10,10,4)`, GPAW Monkhorst-Pack integer mesh, zero shift as reported by the archived run, and no added offset or Gamma override.
- PPCG eigensolver with `niter=2`, matching the archived log.
- Four MPI ranks and one OpenMP thread, matching the archived baseline execution and E runtime.
- Fixed input atoms and cell; no relaxation. Record native total energy and force-consistent free energy separately, all atomic force vectors, stress tensor, SCF convergence/iterations, and input/result hashes. Recording stress for all six is preferred for a uniform output schema; at minimum it is required for the volume and shear variants.

These are proposed/frozen scientific values for E to encode in a final immutable lock. They do not constitute permission to execute.

## Execution blockers

1. `DFT_WORK_ROOT` is unset and the expected durable E work mount is absent. The reported `/workspace` capacity is temporary and cannot alone protect a long run or its only valid results.
2. No final signed archive/restore procedure is supplied for each calculation's output and any restart file. Before execution, demonstrate a fresh remote restore and SHA verification; any checkpoint exceeding the repository's 50 MiB conservative per-file limit needs an approved alternative store.
3. Per-structure wall-time limit, checkpoint cadence/failure policy, and maximum concurrency are not frozen. The archived parent run took 11,722 s at four ranks; this is a planning reference, not a guaranteed runtime for each strained structure.
4. E must publish one immutable final lock binding the six structure hashes and manifest hash, with the explicit method values above, role `DIAGNOSTIC_ONLY`, output contract, work-root/archive rules, and bounded resource policy. The current proposal remains `PROPOSAL_NOT_FROZEN` and has `execution_authorized=false`.

Until those items are resolved and a separate execution authorization binds the final lock, the batch remains `HOLD`. No DFT, fitting, MD, or TTM was run in this review.
