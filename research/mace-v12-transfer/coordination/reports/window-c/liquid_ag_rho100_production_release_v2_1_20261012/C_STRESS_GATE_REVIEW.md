# C review and execution release — Ag64 rho100 QA v2.1

Status: **APPROVED** for the single bounded automatic continuation bound in `C_EXECUTION_AUTHORIZATION.json`. C reviewed D commit `e228f4624e4b022023b995f12b7e05b28c92cd94` fetched from GitHub and verified the runner, criteria, configuration, archive backend, tests, harness, locks, source checkpoint, potential, and environment hashes against the submitted values.

## Pressure and stress decision

The residual circular moving-block bootstrap is accepted for this finite pilot trend screen. The implementation estimates integrated autocorrelation from the OLS-detrended residuals, uses the initial-positive sequence through the first non-positive lag, and chooses `ceil(2*tau_int/dt)` samples per block. It fits each bootstrap slope on the original fixed time axis. This accounts for serial correlation better than an iid bootstrap, while remaining a screening method rather than a general stationarity proof.

A two-sided 95% bootstrap slope interval excluding zero is accepted as the operational definition of resolved drift. An interval containing zero means “no drift resolved at this screen's power,” not proof of stationarity. Other locked temperature, energy, structural, diffusion, dual-window, and end-of-run checks remain in force.

No multiplicity adjustment across pressure and six stress components is accepted because this gate fails closed when any component detects drift. Under independent null components, seven nominal 5% tests could produce at least one false drift flag with probability about 30.2%; the components are dependent, so the actual rate is unknown. The practical risk is unnecessary continuation or HOLD, not permission to pass on a detected drift. These are per-series nominal intervals and must not be described as simultaneous 95% coverage.

The 4 ps minimum interval sampled every 0.01 ps, with at least eight complete autocorrelation-sized blocks, is accepted as a minimum data sufficiency rule for this bounded pilot only. For 400 samples, eight blocks requires block length at most 50 samples; greater autocorrelation leaves the gate `PENDING` until more eligible data exist. Missing/invalid input remains `QA_HOLD`; inadequate duration or block count remains `PENDING`; neither can default to PASS. No absolute pressure threshold is introduced. The pressure/stress gate is wired into the joint liquid QA pass condition.

## Execution bindings and limits

The verified D artifacts match the SHA values in the authorization. D reports 17/17 engineering tests, a synthetic GitHub archive upload/fresh-fetch/restore test, and a read-only live preflight of `PREPARED_HOLD`; the latter is expected before C authorization and must be rerun just before starting. The archive test uses synthetic payloads, so each real chunk still has to pass its own remote upload, fresh-fetch, member-hash, and restored-state checks before the runner can continue.

Authorization is limited to task `C-D-AGL-1250K-rho100-followup01`, resuming the step-1000 checkpoint. The runner may continue automatically only under the locked 1600 K heat, 5 ps cool, 1250 K equilibration (maximum 10 ps), and 5 ps candidate sequence, and only while each chunk is at most 1000 steps / 3600 seconds, the global step is at most 60000, the existing state is restored without reinitialization, and the per-chunk archive gates pass. The EAM, density, cell, timestep, friction, RNG, and existing scientific criteria are unchanged. No other trajectory, DFT, or training is authorized.

C does not execute MD. D may read this release from the C branch and start after validating the C commit/blob SHA and its own live preflight.
