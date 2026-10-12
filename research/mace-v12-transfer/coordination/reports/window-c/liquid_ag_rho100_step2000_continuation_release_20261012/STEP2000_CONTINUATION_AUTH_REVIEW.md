# Ag64 rho100 step-2000 continuation authorization review

**Decision:** `APPROVED` for the exact, hash-bound step-2000 checkpoint continuation. C did not start MD.

## Source and code identity

The D branch was fresh-fetched from GitHub at `4418a3f498803332b521e4870f4c559db1bd462c`, matching the live remote branch head at review. The code-repair commit is `2fb55ac472b2175f39bbe4f18edcd426077e413d`.

The runner, recovery config, criteria, and checkpoint-binding file matched the submitted SHA-256 values. The step-2000 archive commit `d549ef76b1e9f6a208744bea229fad9c581de2e6` is in the D history. Every member in the archived step-2000 `SHA256SUMS.txt` passed verification.

## Runner repair and checkpoint recovery

The fix repairs manifest/checksum writing and adds a hash-checked resume loader. It also changes the runner's initialization path to restore from the signed archive and sets the remaining step budget from the resume config. Review of the diff found no change to the per-step ASE `Langevin.run(1)` integration, force evaluation sequence, timestep, friction, temperature schedule, or locked criteria. The new path adjusts only an in-memory control-stage marker after independently checking the archived QA monitor and thermal log; the checkpoint itself remains immutable.

C independently invoked the recovery loader against the fresh-fetched archive without attaching a calculator or integrating. It restored step 2000 / 1000.0 fs, all 64 atom positions, momenta, masses, cell, PBC, element and atom-ID order, and all three PCG64 states exactly. The archived trajectory's final frame, cell, and PBC match the checkpoint; the thermal log ends at step 2000 with the 1600 K target; the archived QA state is `PENDING` in `A_HOT_1600K`. The recovery changes no physical state and does not treat `PENDING` as a pass.

## Original C authorization boundary

The science, strategy, QA, potential, environment, and resource-limit sections in the step-2000 config match the original locked configuration. The execution-bound changes are the required checkpoint continuation point (`1000` → `2000`), remaining added-step budget (`59000` → `58000`), the new runner hash, and the residual 1600 K hot-stage maximum (`19000` → `18000` steps). That stage-cap reduction accounts for resuming 1000 steps later; the global cap remains step 60000. Chunks remain at most 1000 steps and 3600 seconds. The physical temperature path and automatic 1600 K → 1250 K stage plan are unchanged. The new authorization supersedes the old step-1000 resume pointer only for this continuation and does not widen the original plan.

The authorization requires D to rerun live preflight at the bound D branch head before launch. D's recorded preflight at code-release commit `2fb55ac472b2175f39bbe4f18edcd426077e413d` was `PREPARED_HOLD`, with execution intentionally unauthorized; it reported the locked runtime and resource checks passing. Since the current release request is at D head `4418a3f498803332b521e4870f4c559db1bd462c`, D must refresh that preflight at the bound head before starting. The runner will fail closed if the live environment, branch head, or any input binding has changed.

## Validation note

D's committed recovery report records 20/20 tests passing in its locked environment, including exact checkpoint recovery and fresh-fetch/restore. C independently ran the same suite: 19 tests passed; the OVITO PTM test could not import `ovito` in C's local environment. D's committed environment record identifies OVITO 3.16.1, and D must repeat the runner's live preflight immediately before execution. This C-side dependency difference does not alter the checkpoint or authorize skipping D's preflight gate.

## Signed scope

Only task `C-D-AGL-1250K-rho100-followup01` may resume from the exact bound step-2000 checkpoint. Steps 1000–2000 may not be replayed. Continuation is limited to the existing approved plan through global step 60000, with verified archive/fresh-fetch/restore after each chunk. No other trajectory, training, or DFT is authorized.
