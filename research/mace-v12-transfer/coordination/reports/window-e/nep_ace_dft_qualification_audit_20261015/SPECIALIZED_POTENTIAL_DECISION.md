# Specialized Ag–Ti–Si–C potential data decision

## Decision fields

| Field | Decision | Evidence / limitation |
|---|---|---|
| `DFT_LABEL_INTEGRITY` | `INSUFFICIENT_DATA` | 71 exact geometries have verified converged energy+forces; 21 remain incomplete; method and energy conventions vary; stress/virial coverage is 5/71. |
| `AG_SOLID_DATA_READY` | `INSUFFICIENT_DATA` | Seven Ag FCC controls cover baseline, volume, shear and displacement, but only a narrow solid family. |
| `AG_LIQUID_DATA_READY` | `FAIL` | Zero confirmed liquid Ag DFT structures or solid-liquid coexistence labels. |
| `TSC_BACKBONE_DATA_READY` | `INSUFFICIENT_DATA` | Five bulk force geometries and one small Ti/C displacement family; no confirmed Si displacement, bulk stress response or independent clean surface family. |
| `INTERFACE_DATA_READY` | `INSUFFICIENT_DATA` | 54 confirmed Ag–Ti/Si/C family geometries exist, but parent/method diversity is incompletely verified and liquid-contact/independent test coverage is absent. |
| `INDEPENDENT_TEST_READY` | `FAIL` | V17/V18 reuse the same three dev geometries and have no test set; parent IDs are incomplete. |
| `NEP_TRAINING_DATA_READY` | `FAIL` | Physical coverage, virial coverage, parent-independent test and GPUMD schema/runtime checks are insufficient. |
| `ACE_TRAINING_DATA_READY` | `FAIL` | Same label/coverage gaps; ACE/LAMMPS fitting/deployment benchmark absent. |
| `ADDITIONAL_DFT_REQUIRED` | `YES` | Liquid Ag, mixed liquid-solid contact, backbone elastic/Si, and independent surface/interface family are absent. |
| `RECOMMENDED_FIRST_ROUTE` | `ACE, conditional` | Best fit to the present CPU/LAMMPS-oriented path without an exposed GPU; this is a route to evaluate after data qualification, not a claim that ACE is production-fast at 700k atoms. NEP remains a candidate if an approved GPU environment is available and passes benchmark. |

## Practical next gate

1. Recover raw SCF/run metadata for the 21 incomplete geometries and normalize energy/free-energy and units without changing frozen roles.
2. Freeze parent-family TRAIN/DEV/TEST roles; obtain at least one genuinely independent sealed family.
3. Complete the prioritized DFT additions in `MINIMUM_ADDITIONAL_DFT_PLAN.md`.
4. Validate an extxyz conversion with the actual NEP/ACE readers and run a small target-hardware benchmark.
5. Only then compare exploratory ACE and NEP fits. Neither is training-ready now, and neither is production-ready for approximately 700,000 atoms.

`RECOMMENDED_FIRST_ROUTE = ACE (conditional after the gates above)`.
`PRODUCTION_READY = NO`.

No new DFT, training, relaxation, MD, TTM or potential modification was performed in this audit.
