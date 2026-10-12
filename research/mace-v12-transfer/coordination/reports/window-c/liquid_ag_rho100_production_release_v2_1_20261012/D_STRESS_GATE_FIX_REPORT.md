# D stress gate repair and C release request

## Scope

This patch repairs only the Ag64 rho100 QA v2.1 pressure/stress gate. It changes the production criteria module, its runner input path, the method/configuration hashes, and synthetic regression evidence. It does not alter the EAM potential, density, cell, integration step, temperature schedule, checkpoint, RNG, or any other scientific criterion. The configuration still has `execution_authorized=false`; no MD, DFT, or training was run.

## Implemented gate

At each fixed-temperature phase boundary, the runner joins `stress.csv` to `thermal.csv` by global step, verifies all six named ASE Voigt components and sample times, and retains only rows after the phase's already frozen transient discard and at the exact target temperature. The stationary interval is all such eligible samples through the current boundary; at least the existing 4 ps QA interval and eight complete autocorrelation-sized blocks are required. The logged cadence must remain 0.01 ps; missing or irregular samples hold the gate.

The stress component order is `xx, yy, zz, yz, xz, xy`, in eV/Å³. Installed ASE documentation states stress is positive in tension and pressure is positive in compression. The gate calculates hydrostatic pressure as `P=-(sigma_xx+sigma_yy+sigma_zz)/3`; it reports pressure in eV/Å³ and GPa. It sets no absolute pressure acceptance value.

For pressure and each of the six components, the report includes the sample count, mean, two-sided 95% moving-block-bootstrap interval and half-width, OLS slope per ps, two-sided 95% slope interval, estimated integrated autocorrelation time, block length, and effective block count. The provisional method uses the frozen circular moving-block bootstrap family, 2,000 deterministic replicates, block length `ceil(2*tau_int/dt)`, and an initial-positive-sequence estimate of `tau_int` on OLS-detrended residuals. A two-sided 95% slope interval excluding zero marks resolved drift; any such drift fails the pressure/stress subgate. This test specification and per-observable multiplicity policy are explicitly submitted to C for confirmation; they are not represented as C-approved until C confirms them.

Missing or malformed columns, non-finite values, duplicate/mismatched samples, or irregular timestamps return `QA_HOLD`. A stationary interval shorter than 4 ps or with fewer than eight effective blocks returns `PENDING`, allowing only the already bounded runner to collect its next internal chunk. The joint liquid assessment cannot pass around this gate.

## Test results

`test_production_runner_v2_1.py` retains its original 12 cases and now passes **17/17**. Added synthetic checks establish:

- Constant pressure and all six constant components: `PASS`, including sign and units.
- Known linear stress/pressure drift: `FAIL_DRIFT`.
- Correlated AR(1) data: block-bootstrap uncertainty exceeds the iid standard-error estimate and uses blocks longer than one sample.
- Insufficient duration/effective blocks: `PENDING`; missing values, a missing component, and NaN: `QA_HOLD`.
- An altered logging cadence: `QA_HOLD`.
- Other liquid metrics passing cannot override a stress drift.
- A fresh archive restore yields byte-equivalent samples and an identical deterministic stress-gate result to the uninterrupted input.

The read-only production preflight reports `PREPARED_HOLD`, confirms the updated runner/criteria/config hashes, and verifies the original step-1000 checkpoint SHA, EAM SHA, ASE 3.29.0 parser, runtime environment, resources, and existing RNG restoration. No step was executed.

## Bound artifacts

- Runner SHA-256: `f26be9819baedc6cad5c3d1ed36ee0e7dba33be2786b09ce8bacf826e3a1c6bd`
- Criteria SHA-256: `dadf5647cd1ffb68fd421e8f2742117c604b324e20443c0e2551e7208ae21c69`
- Configuration SHA-256: `3af1d791a0e8afb601bfd10822a56ad09b2fe0d87580e5e16758057674c86c1c`
- Environment SHA-256: `e71e04b8edc989001a702ffb308b1b897d6b521d0ebfdcb41b43592ae2bff1fb`
- Source checkpoint SHA-256: `c820437a2512c4b3d7dc9e60bb60ac3ced6588aa1f2a629d58cfdabb24f34046`
- Ag EAM SHA-256: `958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447`
- Automated regression: 17 passed, 0 failed; formal MD steps: 0
- Test suite SHA-256: `eecdc40f00abaa72f33243471e89df2477c77e665c096225418fa70da4d3d5f9`

## C review requested

C should review the final runner and criteria SHAs and explicitly confirm or reject:

1. Residual circular moving-block bootstrap for the OLS slope interval, with `tau_int` estimated from detrended residuals and `ceil(2*tau_int/dt)` blocks.
2. The two-sided 95% interval excluding zero as the resolved-drift rule.
3. Applying that rule separately to pressure and six components with no multiplicity adjustment, failing closed if any one resolves drift.
4. Treating the current 4 ps minimum stationary interval and eight complete correlation-sized blocks as sufficient-data requirements for this subgate.

Only after C accepts the implementation and issues a separate signed release bound to the final D commit, runner, criteria, config, environment, checkpoint, and potential hashes may the authorized bounded MD begin. This patch does not authorize or start that run.
