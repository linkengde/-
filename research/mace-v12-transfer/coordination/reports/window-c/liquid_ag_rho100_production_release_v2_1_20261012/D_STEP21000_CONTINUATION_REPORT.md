# D Stage A continuation report: step21000 → global cap 60000

Status: prepared and validated; execution resumes only from the immutable, remotely verified step21000 archive. No historical step is rerun.

- D branch base before publishing this package: `66fa2bf1ac27939098042f3d5ab277f22ac654e5`
- Task: `C-D-AGL-1250K-rho100-followup01`
- Step21000 checkpoint SHA-256: `d9cdf98a08eedaaaabd740f5d23b6951d33d644612feeba836615a650e7d5efa`
- Archive commit: `66fa2bf1ac27939098042f3d5ab277f22ac654e5`
- New runner SHA-256: `a4f17d3732bc9b5d98bc3057934a6bd11293616db889c7b5ca5282d2cab51d84`
- Config SHA-256: `ba45d443ed7b67b759e130203428e756818fe69f220c787d16d1f10d88c0ae1b`
- Step21000 binding SHA-256: `e6b1e0a4f364ca1207b26601a9881635d243a1d268fbd7f59ff43fa5a2c7581b`
- Direct user continuation authorization SHA-256: `acfd9fb4987f6fd308e3ec3cf43e8ea9146328774b787927a6afae9a9613f6b3`
- Environment SHA-256: `e71e04b8edc989001a702ffb308b1b897d6b521d0ebfdcb41b43592ae2bff1fb`
- EAM SHA-256: `958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447`
- Frozen QA v2.1 lock SHA-256: `848f8db72d83cea91dcf4cd42ee9fbbf49a80e7b5e3759cb4f535be8cc5873d2`

## Stop policy

- Stage A at 1600 K: `PASS` transitions immediately to the existing frozen cooling stage; `PENDING` continues at 1600 K; `QA_HOLD` stops immediately.
- Global safety maximum: step 60000.
- Each internal chunk remains at most 1000 steps and 3600 seconds, followed by checkpoint, upload, fresh-fetch, restore, and SHA checks.
- EAM, density, structure, RNG, timestep, temperature, QA thresholds, and QA implementation are unchanged.

## Verification and interruption handling

- Six non-integrating tests pass, including exact restoration of positions, momenta, and all RNG streams; Stage A PENDING/PASS/HOLD transitions; the step60000 bound; and prevention of duplicate writes to an immutable archive boundary after a graceful signal.
- Live preflight: `PREPARED_HOLD`. Code/config/environment/EAM hashes, ASE EAM parsing, source and resume checkpoint restores, historical Q6 regression, resources, and fresh output-directory isolation pass. No integration occurred during this preflight.
- The earlier step20000→21000 segment completed 1000 steps and produced a fresh-fetch-verified remote archive at step21000 (checkpoint SHA above). QA was `PENDING` due insufficient RDF uncertainty data.
- A graceful signal arrived at the boundary after the step21000 archive was already pushed. A repeated attempt to write the same immutable archive path was rejected; the original remote archive remains intact. The new runner avoids retrying a boundary that is already archived. Recovery uses the remote step21000 checkpoint, not the local post-stop status rewrite.

The config keeps `execution_authorized=false`; the separate SHA-bound authorization records the user's current stop policy and controls execution. No QA threshold is changed or treated as PASS by default.
