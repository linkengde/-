# D Ag64 step 1300 continuation release

This package contains a hash-bound continuation adapter and C release config for exactly one task: `C-D-AGL-1250K-rho097-seed01`, continuing its archived state from global step 1300 through a hard stop at step 2000. It does not authorize step 2001, any other task, or continuation to 120000.

The source checkpoint was fetched from D's `cloud-d-b1-dft` branch at pinned head `9aea23179dd73da749d67f7be0b3a3b86e5df085`. Its archive SHA-256 is `e75322b9bdd4b6bc2d485225a76c889536219ac22b6fe299858f96c1411839ee`; `state.json` SHA-256 is `3b3950729735765dc5c4d205d16a3654b4e9f8ba193c58bf8421685f36d66b9e`; `state_arrays.npz` SHA-256 is `6eb03c97a33f6c2f3b5d9521c98435a77e5ee9228286090602733c3d22a65a2a`. The arrays, all three PCG64 stream states, 15 thermal rows, and four trajectory frames were restored read-only. Coordinates, momenta, velocities, cell, PBC, atom order, and RNG are carried forward without reinitialization.

The original scientific settings remain pinned to scientific-parameters SHA `8f792d6206c1b7bae7b900e7ea81445d6cd3e837de2b5979c88bb977650ab7b2`, and the Ag potential must match SHA-256 `958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447`. The signed walltime cap remains 600 seconds. If it expires before step 2000, the runner preserves the latest verified archive and stops; it does not retry or continue automatically. A further continuation needs a new review and release.

The archive's `execution_status` field is retained verbatim as `SYNTHETIC_OR_UNSIGNED_CANDIDATE_ONLY`. Its provenance separately records the prior C-authorized step-1000 release SHA and `synthetic=false`, and the user reports that D's step-1000-to-1300 execution audit passed. This mismatch is disclosed and is not rewritten in the old archive.

Before any D execution, the adapter must fresh-clone D's branch and find the exact pinned head and source archive, verify all archive and internal hashes, confirm the local potential hash, confirm the signed output directory is absent and outside the Git checkout, then construct ASE's EAM calculator. Each new 100-step checkpoint must upload to D's branch, fresh-fetch, restore, and compare positions, momenta, velocities, and RNG state before another integration step. The potential file is never archived.

`test_offline_release.py` performed read-only source restoration, config/hash and identity refusal checks, and a step-2000 guard check that refuses before calculator/force access. It executed zero MD steps. No D continuation, DFT, or MACE training was started by C.

`SIGNED_CONTINUATION_CONFIG.json` is the release artifact; verify it against `SIGNED_CONTINUATION_CONFIG.sha256` before use. Its C signature is procedural and audit-bound to the immutable Git commit and SHA, not a cryptographic identity signature.
