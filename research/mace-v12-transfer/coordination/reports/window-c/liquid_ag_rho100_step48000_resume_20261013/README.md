# Ag64 rho100 step 48000 recovery entry

This isolated runtime entry restores the immutable, fresh-fetched step-48000 archive and preserves its positions, momenta, three PCG64 streams, and Stage A state. It continues only the user-authorized Stage A rule: PENDING continues at 1600 K, PASS transitions to the frozen cooling stage, QA_HOLD stops, and global step 60000 is the hard cap.

No scientific lock, potential, density, cell, integrator, thermostat, timestep, RNG state, or QA threshold was changed. The packaged execution config remains `execution_authorized=false`; the separate SHA-bound authorization records the user's current-session instruction. Each internal segment remains limited to 1000 steps and 3600 seconds, and the production runner requires checkpoint archive, fresh-fetch, and restore validation before continuing.

The step-48000 archive was independently fresh-fetched from D branch commit `d081d78f07c799fbd8b79102461104f6f9177ae2`; all 18 files in its SHA256SUMS passed. The checkpoint SHA is `f21165280a3524612ab51702f7b28ee56643aace598d484f2966d0ce25291b57`.
