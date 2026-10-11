# Ag64 step2000 review-gap validation report

## Scope

Candidate config SHA-256: `a6e794f208cc9b874b3dfa8e9c1ce0359c134f2babe9e768b6249d7d17fa8ac5`  
Task: `C-D-AGL-1250K-rho097-seed01`; target segment `2000→3000`; `max_added_steps=1000`.  
Branch snapshot used for environment identity: `ag64-step2000-extension-dev` at `538f23e6fc2b3b3461693085487e065325229e5e`.

No real MD, production calculator, production force or energy calculation, or production potential load occurred. The timeout probe used only virtual integer steps and constant-valued atom test doubles. The source checkpoint, signed step1900→2000 package, sampling settings, and potential bytes were not modified. No execution-authorized config was created.

## Gate fields

```text
ENVIRONMENT_EVIDENCE_PASS=false
BACKEND_V2_REMOTE_RESTORE_PASS=false
TIMEOUT_RECOVERY_PASS=false
```

The three requested evidence artifacts are attached beside this report. Their subchecks and blockers are separated below.

## 1. Environment evidence

`ENVIRONMENT_EVIDENCE.json` was generated and its normalized payload SHA-256 is `a751d639e67fe282ff65bf1d11813be79b284f7e35150ebfe09cf4534ab7512a` (canonical UTF-8 JSON, sorted keys, compact separators, excluding the digest field itself).

Observed environment: Python `3.12.14`, NumPy `2.3.5`, ASE `MISSING` (required `3.29.0`). The backend imports, while adapter/core import is blocked because ASE is absent. The candidate pins adapter/core SHA `581fea79…` / `29fd90b0…`; the target branch snapshot’s files are `4dfc0485…` / `adcc939f…`, so the code binding does not match. Candidate `ENVIRONMENT_SHA=14084cfa…` cannot be reproduced from this workspace. Evidence payload hashing passes; candidate environment validation does not.

## 2. Backend v2 remote restore

The pinned backend SHA is `0a8d8b71…`. `GitHubArchiveBackend.upload()` was called with the 2000→3000 v2 binding, but the call stopped at `git push`: this runtime has no GitHub CLI write credential (`could not read Username`). Therefore the requested full `upload() → verify() → fresh-clone fetch() → restore()` backend path is **not passed**.

To still exercise the real remote read/restore path, the same transport fixture and v2 receipt were committed through the connected GitHub Git Data API to isolated branch `codex-ag64-v2-restore-test-20261011-a6e794f2`. The backend’s `verify()` and `restore()` then passed from fresh Git clones. `RESTORE_STATUS=PASS`; receipt SHA-256 `50fb96e3200d00892ad8c5e9ff4fa11ff211a6268afd6a9490d2700f517efc58`, manifest SHA-256 `743b068ca781f03a8bf4413e7c7dfcb8074532a42a6f9429c15ad15e1c87af86`, restored synthetic state SHA-256 `29598cb45bd49e7baaefb139d4d07fb6eb09e20efd8f3d258c518eb23c22170b`, arrays SHA-256 `0f7830eb…`, and RNG SHA-256 `43030996…` all verified. The arrays and RNG bytes match the source step2000 copy; the receipt binds source checkpoint SHA `5d03864f…`. The synthetic step3000 metadata exists only to meet the backend v2 range guard; it is not a computed output. The test branch was reset to the D base head after the test.

## 3. Timeout recovery

The exact candidate core SHA `29fd90b0…` was run with the existing `advance_segment_with_outputs()` runner and virtual integer-only dynamics. At the simulated wall-time stop, the runner halted at step 2500 with `SEGMENT_WALLTIME_BUDGET_EXCEEDED` and created **zero checkpoint files**. This reproduces the C blocker.

A separate manual supervisor save using the existing checkpoint bundle writer at step 2500 restored step, arrays, and all RNG streams; resuming the same runner advanced virtual steps `2501..3000`, produced one final bundle, and restored step 3000. The total tested sequence `2001..3000` was contiguous with no repeated or out-of-range steps. This proves the bundle restore/resume path when an external actor saves, but it does not make the runner’s timeout path recoverable by itself. Therefore `TIMEOUT_RECOVERY_PASS=false`.

## Current blockers and C review status

- ASE 3.29.0 is not available in the current runtime, and the claimed environment SHA has no matching evidence payload.
- The backend’s direct `upload()` call could not authenticate; only real GitHub API seeding followed by backend fresh-clone verify/restore passed.
- The runner does not create a checkpoint on wall-time timeout before the 1000-step archive boundary.
- Candidate adapter/core pins do not match the code currently on the development branch.

**Do not request C signing yet.** The evidence identifies the remaining gates; this report does not authorize execution. No source-code or config files were modified by this validation.
