# Ag-Si Registry Screening v2: 6x6x1 Execution Report

- Scope: R1-R4 only; fixed-geometry GPAW/PBE/PAW single points.
- Frozen settings: 500 eV plane-wave cutoff, 0.10 eV Fermi-Dirac smearing, Monkhorst-Pack 6x6x1, no atomic or cell relaxation, 4 MPI ranks, 1 OpenMP thread.
- Configuration SHA-256: 515cd81f16e491d8d521426ba63f17946cb47dd5fa0361fa9a8d946e07f0d035.
- Authorization SHA-256: 86b5058f12745e1484680bfc926659b3a58fe426a873076514ef614467fc6b9f.
- Role: REGISTRY_SCREENING_ONLY.
- 8x8x1 was not run; no automatic follow-on was launched.
- Result: no registry produced a complete result.json; no energy ranking is valid. Ranking status is INSUFFICIENT_FOR_8X8.

| Registry | input.py SHA-256 | structure.xyz SHA-256 | Status | SCF status | Exit |
|---|---|---|---|---|---:|
| R1 | d3f7b5542f071cbc432ade3277a22b0af5277cfd539d1eb48a0d006e7250ca6b | b758b39571a6ce8b24f687e9b341c4328336c108d761b4adaa25e3fb51932890 | FAILED_SCF | FAILED_SCF | 42 |
| R2 | 7c1327f403c94240e6181da990fe7b078f420cc32f92a68bee3282792b0bcc66 | cb2e16376887f7cdeb05c7ddf4aacfacacbea0a000837a3888a0f3ea91785f70 | FAILED_POSTPROCESS | CONVERGENCE_UNVERIFIED | 42 |
| R3 | c9b731143aafe2fb7ba367093cad72ad1489a4a812aa63b4c4ccba53256f33e2 | edc7a00668be7fb0f76616463d3eacf2ac53d0ea8837b21fb9fb3c66eafa6f76 | FAILED_SCF | FAILED_SCF | 42 |
| R4 | 1a8fca555c0dd7a7028e8f223a54770ef5183549f0a9c230869ed88dbc4198e8 | 6465b1e2ea8c675a4d368a6f385cc31302a5c4b017d9b39842e82c93c196210c | FAILED_POSTPROCESS | CONVERGENCE_UNVERIFIED | 42 |

Failure details:

- R1 and R3 reached the frozen 160-iteration limit and raised GPAW KohnShamConvergenceError. Their complete SCF logs and partial output are retained.
- R2 and R4 reached force/stress evaluation, then the frozen input failed at its convergence-status check because the GPAW/ASE calculator exposed as ASECalculator has no .scf attribute. The traceback does not establish SCF convergence; no result JSON was written, so these are classified FAILED_POSTPROCESS / CONVERGENCE_UNVERIFIED. Their complete logs and partial output are retained.
- No input, structure, cell, k-point, smearing, mixing, cutoff, or role was modified. No rerun or convergence repair was performed.

Each registry's stdout.log, stderr.log, results/coarse/gpaw.log, exit code, and status JSON are retained. SHA256SUMS.txt covers this execution package.
