# D Ag64 ASE segment integration candidate

This C-only directory contains an ASE 3.29.0 restart adapter and tests for the single bounded continuation `C-D-AGL-1250K-rho097-seed01`, global step 1000 to 2000.

Current status is `INTEGRATION_TESTS_PASS / RELEASE_UNSIGNED / REMOTE_GITHUB_GATE_HOLD`. The adapter, tests, and candidate lock are new. They do not modify the old `run_sampling_codex.py`, the D step-1000 checkpoint, or its archived files. No Ag potential, coordinates from the D checkpoint, or trajectory data are copied here.

The synthetic integration tests use ASE 3.29.0 with a marked virtual harmonic calculator. A separate read-only probe reconstructs the D checkpoint in memory and verifies exact state equality; it creates no calculator and evaluates no energy or force.

`SOURCE_RESTORE_PROBE.json` records the checkpoint hashes and exact restoration results. `ASE_API_SEMANTICS.json` records the inspected ASE implementation hashes. `HASH_SCOPE_AUDIT.md` explains the historical identity digest. `CANDIDATE_RELEASE_LOCK_UNSIGNED.json` is explicitly not authorization. `LAUNCH_COMMAND_DRAFT.md` is a blocked command template, not an instruction to D.

The local fresh-extract archive tests are synthetic. A real GitHub push, fresh fetch, restore, and hash verification for a new step-1100–2000 archive has not been performed, so production execution remains HOLD.
