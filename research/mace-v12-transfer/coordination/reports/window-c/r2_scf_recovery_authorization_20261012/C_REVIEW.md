# C review — E R2 single-point SCF recovery

**Decision: APPROVED. Checkpoint test status: ACCEPT_LIMITATION.** This authorization covers one R2 fixed-geometry recovery single point only. No DFT was run by C.

## Fresh-fetch and source identity

C fetched E branch head `4d97374aefc338a714411003aa06ec35aa894611` from GitHub. Runner, collector, and config hashes match the requested SHA-256 values. The proposal SHA is `2d16458e0d3139722609746c067cb65a9577533a29db04282213d69da4a0bf2f`. The Phase2 baseline commit is `0dd8628e4a91e88f31aba3fd7bbce1c8a40a2aa2`. Comparing that baseline with E's final commit, the R2 `input.py`, frozen R2 `structure.xyz`, Phase2 parameter lock, and archived Phase2 R2 result have identical SHA-256 values. The runner also rechecks source input, structure, geometry, PBC, and lock at launch.

## Method and scope

The recovery keeps PBE/PAW, 500 eV, 8×8×1, Fermi-Dirac 0.10 eV, the 1e-5 convergence thresholds, and the Phase2 `Mixer(0.05, 8, 100)`. The authorized recovery settings are PPCG `niter=4`, max SCF 200, four MPI ranks, one OMP thread, and a 14,400-second wall limit with orderly stop requested at 12,600 seconds. Geometry and cell remain fixed. The scope and collector role remain `R2_ONLY_SINGLE_DIAGNOSTIC` / `REGISTRY_SCREENING_ONLY_DIAGNOSTIC_COPY`; no R1, other registry, d0 scan, geometry edit, retry, or training is authorized.

## Checkpoint limitation and failure handling

`PASS_WITH_LIMITATION` is accepted. A serial real-GPAW synthetic full-wavefunction GPW write/readback passed, the four-rank write-error synchronization test passed, and one synthetic four-rank readback completed. A later optional synthetic four-rank initialization repeat did not reach the hook; it is inconclusive rather than a failed checkpoint. This leaves some risk that a successful four-rank checkpoint in the real R2 run could fail, but the runner's failure path synchronizes MPI ranks, preserves earlier verified checkpoints and logs, stops without retry, and the collector withholds all scientific labels from incomplete/unconverged output.

The collector only emits `SUCCESS` after matching a GPAW convergence marker, runner and raw success status, converged iteration counts, registry/role, authorization/config/source hashes, finite energy and free energy, and finite force/stress arrays of the expected shapes. SCF nonconvergence is marked `FAILED_SCF`; incomplete/error output cannot be promoted to a label.

## Disk gate

The recorded preflight reports 27,939,389,440 bytes free versus a 23,353,522,176-byte budget, leaving 4,585,867,264 bytes above the required reserve; inode, CPU, and memory checks also passed. That snapshot was produced on the prior package base commit, so the runner must repeat its live preflight immediately before starting. It refuses to start if the current space is insufficient and rechecks before each checkpoint, stopping if one estimated checkpoint plus the reserve no longer fits.

The C authorization binds E's exact commit, proposal, config, runner, collector, Phase2 parameter lock, R2 input and structure, resource report, test report, and manifest. It authorizes no Phase2 modification. E must fetch this C commit, verify the authorization blob SHA, and run its normal preflight/authorization gate before launching the single point.
