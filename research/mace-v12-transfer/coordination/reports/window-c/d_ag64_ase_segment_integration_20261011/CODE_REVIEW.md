# Code review and current release decision

Review result: the offline restart adapter is internally consistent for ASE 3.29.0 and the bounded segment. The unit suite exercises the actual ASE `Langevin.step()` implementation, not an independently counted wrapper. It proves the stop occurs before the integrator can consume RNG or calculate another force at global step 2000.

The source checkpoint restore validates the outer archive manifest, all 18 archived members, D task/owner/step identity, checkpoint and array SHA, source-manifest SHA, atom sequence, composition, cell, PBC, masses, positions, momenta, and all three PCG64 states. It does not attach a calculator. The separate source probe confirms byte-exact restoration.

The bounded step path restores RNG and atom state after an injected nonfinite calculator response. Timeout checks occur before the next step. A checkpoint/archive callback failure terminates the loop and leaves the step-1000 source untouched. `run()` and `irun()` are disabled.

Two blockers prevent production use:

1. The tests use a local synthetic archive/fresh-extract callback. There is no verified GitHub push/fresh-fetch for new segment output from D, so the lock validator requires `remote_archive_verified=true` and refuses the current candidate.
2. The CLI has no real GitHub archive backend. Even a future signed lock cannot run through the current CLI; it explicitly returns `VERIFIED_GITHUB_ARCHIVE_ORCHESTRATION_NOT_IMPLEMENTED`. That backend and its D-node readback test must be added and reviewed before any C signing.

There is also a source-data label defect: `state_arrays.npz` names one array `velocities_A_per_fs`, while the stored numbers equal raw ASE `get_velocities()` and `momenta/masses`. The adapter restores momenta as the authoritative state and proves the reconstructed velocity exactly matches the stored bytes. New bundles use the unambiguous key `velocities_ASE_internal_units`; no unit conversion is silently applied.

Decision: `CODE_REVIEW = PASS_FOR_SYNTHETIC_ADAPTER`; `RELEASE_CANDIDATE = PREPARED_UNSIGNED`; `NEW_MD_RELEASE = NO`; `D_CONTINUE = HOLD`. The code is frozen for review only, not ready for D to execute.
