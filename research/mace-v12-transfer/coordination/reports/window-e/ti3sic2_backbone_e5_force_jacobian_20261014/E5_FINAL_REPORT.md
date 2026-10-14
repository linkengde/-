# Stage E5 final report — Ti₃SiC₂ TSCBOP force-response Jacobian

## Findings

- Eight tied scalar blocks were tested separately at −1%/+1%; 102 single-process LAMMPS `run 0` evaluations include common baselines. Each block changes only the six rows sharing its unordered pair/field.
- E4 DFT-force labels were hash-verified and contain 354 atoms; all six new baselines reproduce E4 forces exactly. The source TSCBOP SHA remains `729273a6488d75756ed7927432c77c7347a1b115cecbaefae4981bafcec80fd7`.
- Full Jacobian shape: `1062×8`; rank **6/8**. `A_Si-C` and `B_Si-C` are zero-response columns because every available Si–C distance exceeds their outer cutoff. Active six-column rank is 6, but normalized condition is 1104.7; effective rank falls to 5 at 1e-3 and 4 at 1e-2.
- Strong coupling prevents robust simultaneous use of the four Ti–C radial blocks and Ti–Si A/B. The best-conditioned active two-column slice is `A_Ti-C, B_Ti-Si` (condition 3.288); this is a local response-screening result only.
- The five 48-atom cases are one parent family; the one 114-atom slab is derived. Thus there is no independent validation family and no claim of general parameter identifiability or predictive performance.

## Required decision fields

| Field | Result | Basis |
|---|---|---|
| `PARAMETER_BLOCK_MAPPING_PASS` | YES | E3 shared-pair map, original file, exact tied-row edits, and LAMMPS implementation mapping |
| `RAW_FORCE_DATA_AVAILABLE` | YES | Six converged E4 geometries, actual hashed force labels |
| `RUN0_REPRODUCED` | YES | All six E5 source baselines match E4 per component exactly |
| `FORCE_RESPONSE_COMPUTED` | YES | Eight columns, 354 atom vectors, 1062 force components |
| `JACOBIAN_IDENTIFIABILITY` | PARTIAL / ILL-CONDITIONED | Six active algebraic directions, two inactive columns, high correlations and no independent validation family |
| `IDENTIFIABLE_PARAMETER_SET` | `A_Ti-C, B_Ti-C, lambda1_Ti-C, lambda2_Ti-C, A_Ti-Si, B_Ti-Si (local force-active only); Si–C A/B not constrained` | Rank and distance/taper audit |
| `TRIAL_FIT_RECOMMENDATION` | NO | Leave all fitting on hold pending data coverage and independent family roles |
| `PRODUCTION_READY` | NO | No fit and no independent validation |

The output package includes raw per-atom forces for every baseline/variant, parameter response matrix, per-element/per-structure summaries, correlations, SVD analysis, finite-difference stability, pair activation checks, and SHA-256 manifest. No fit, DFT, MD, TTM, or MACE training was run.
