# PHASE2 Collector and Runner Test Report

- Status: `PREPARED_NOT_ISSUED`; execution authorization remains false.
- Scope: `R1`, `R2` only; fixed-geometry single points; `8×8×1`; PBE/PAW; MPI 4; OMP 1; role `REGISTRY_SCREENING_ONLY`.
- DFT started: **NO**. No structure or source input was changed.

## Validation

- Mock JSON result parsed total energy, free energy, a 3×3 stress tensor, forces, SCF convergence and iteration count.
- Mock GPAW log parsed the `8×8×1` mesh, convergence in 9 iterations, total/free energy and forces.
- A `6×6×1` result was rejected for Phase2 and its labels withheld.
- Input SHA mismatch was rejected and result labels withheld.
- Authorization gate rejected a correctly hashed `PREPARED_NOT_ISSUED` auth document.
- Request bindings matched the current runner, collector, parameter-lock, input and structure SHA-256 values.
- Summary JSON contains the required Phase2 fields and only `R2`, `R1`.
- Test suite: **7 passed**.
- Source structure/input hash preflight: **PASS** for R1 and R2.
- Collector contains no `ASECalculator.scf` access; SCF status and step count are parsed from GPAW logs or explicit result metadata.

## SHA-256 bindings

| Object | SHA-256 |
|---|---|
| Phase2 runner | `2f5f0a449bce86320275f3d18df785686327ca0ab49664cc05e484a6b79d916a` |
| Result collector | `15592f6ac93755a1b0870574d49d9cac5d856c52ec48ac42a17bd82df375469f` |
| Parameter lock | `d97de3fe5a26ab6529d5e6e8ceef15adc1f3a2f548f03b59f9f092361c2c0c78` |
| Authorization request | `8aa15d1f948ddbbd3cfa6b64033659ab106c4cfe838b39c24cc02eba0d2a12b4` |
| Results summary | `7e2f834770f7ef3567469f671da82b727651e08140b83e73e036de1d66e377ee` |

Per-registry immutable source bindings:

| ID | Input SHA-256 | Structure SHA-256 |
|---|---|---|
| R1 | `d3f7b5542f071cbc432ade3277a22b0af5277cfd539d1eb48a0d006e7250ca6b` | `b758b39571a6ce8b24f687e9b341c4328336c108d761b4adaa25e3fb51932890` |
| R2 | `7c1327f403c94240e6181da990fe7b078f420cc32f92a68bee3282792b0bcc66` | `cb2e16376887f7cdeb05c7ddf4aacfacacbea0a000837a3888a0f3ea91785f70` |

## Current collection state

The collector reports `PARTIAL_OR_PENDING`. No Phase2 `result.json` or GPAW log is present for either allowed registry, so energies, free energies, forces, stress and SCF fields remain unavailable; `config_sha` remains null until a C-issued phase-2 execution record exists. Phase1 6×6×1 results are deliberately not substituted.

The Phase2 runner defaults to read-only preflight and will continue to HOLD until C issues `REGISTRY_PHASE2_EXEC_AUTH.json` with status `ISSUED` and the exact bindings in `PHASE2_EXEC_AUTH_REQUEST.json`.
