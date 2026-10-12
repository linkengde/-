# D step20000 → 40000 continuation entry report

Status: **prepared, not executed**. The new isolated entry restores the immutable step20000 checkpoint and caps Stage A at global step 40000 (at most 20000 added steps). No MD, DFT, or training was started.

## Bound inputs

- Task: `C-D-AGL-1250K-rho100-followup01`
- Owner: `D`
- D branch base tested before packaging: `969ef314c51c7e82c20312275ca5d27e95da6bc5`
- Step20000 checkpoint SHA-256: `826f7d663506dbe7e94237c9955696a3cb35f9cd3da9d8297845c39347526849`
- Step20000 archive SHA256SUMS SHA-256: `ffe911677b44e6467273d1eea9f3057f979b6943228a317018b2663ebce04d28`
- Ag EAM SHA-256: `958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447`
- Scientific lock SHA-256: `e1adaf12d542184c709d2eb2a5363cbfe9578c168ebc5e7ad260ee39a4faf09b`
- QA v2.1 lock SHA-256: `848f8db72d83cea91dcf4cd42ee9fbbf49a80e7b5e3759cb4f535be8cc5873d2`
- New runner SHA-256: `a6da5ac7b020aa7dbbe61936606b75e12af07cb7458c64a4a0cd7629cb61b8e9`
- Config SHA-256: `5280b4402849f892d0299bcf1b735b73f4bfae56048023da559cc4a4db73d14b`
- Step20000 binding SHA-256: `d117016da18b3e6ef3531c30a1b855d9d59a02ea823822c9798e98362dcbf092`
- User bounded-extension instruction record SHA-256: `34c44b8fcd784e0714b49f73d47bb1cfb43740fcabf6803a1e0198296ea85079`
- Runtime environment SHA-256: `e71e04b8edc989001a702ffb308b1b897d6b521d0ebfdcb41b43592ae2bff1fb`

The original Stage A temperature (1600 K), EAM, density, timestep, coordinates, velocities, RNG streams and QA criteria are bound to the existing checkpoint/configuration and were not edited. The new continuation config has `execution_authorized=false`; this is an entry-preparation artifact and does not authorize or start integration.

## Verification

- Unit tests: **4/4 PASS** (exact bounds, Stage A/global cap, bitwise state restore, restore does not create new structure/integrator state).
- Live preflight: **PREPARED_HOLD**; EAM parser, code/config/environment/QA lock hashes, ASE/GPAW/MPI runtime, checkpoint member hash, positions/momenta/RNG restore, step boundary, output isolation, and resources all pass. `formal_md_started=false`, `formal_steps=0`.
- Historical global-Q6 regression: **PASS**, 11 frames; local Q6 remains report-only.
- Fresh GitHub clone/fetch from `969ef314c51c7e82c20312275ca5d27e95da6bc5`: **PASS**. Step20000 checkpoint archive member SHA, positions, velocities and all RNG streams restored; stage restored as `A_HOT_1600K`; no MD started.
- Current observed storage: `28198395904` free bytes; memory limit `34359738368` bytes; CPU quota `4.0` cores.

## QA state and execution boundary

The source checkpoint's QA is **PENDING**; its earlier report also records a coordination-number result above the frozen upper bound in one evaluated window. This extension preserves the frozen gates and fail-closed behavior. The package does not convert PENDING/FAIL to PASS, relax thresholds, or run MD. The extension entry is ready for review/use under the stated cap, with execution still held by `execution_authorized=false`.
