# Ag64 step2000 pre-production validation report

## Scope

Validated the proposed offline step 2000→3000 continuation using the immutable D step2000 archive. No MD was run, no production potential was loaded, and no dynamics integration step was advanced. The candidate template is unsigned and has `execution_authorized: false`.

## Files added

- `tests/test_step2000_preproduction.py` — six offline unit tests, using an ASE stub because ASE is absent in this runtime.
- `ENVIRONMENT_CHECK_REPORT.json` — environment, import, archive, dry-run, and test status.
- `CANDIDATE_CONTINUATION_CONFIG_STEP2000.json` — unsigned step 2000→3000 template with pinned input and code hashes.
- `REPORT_STEP2000_PREPRODUCTION.md` — this report.

The adapter, core, archive backend, potential, and signed step1900→2000 package were not modified in this validation pass.

## Validation results

| Check | Result | Evidence |
| --- | --- | --- |
| Step2000 archive SHA | PASS | Matches `c3a88a759e0bd258e74632aba0a68e4f1a765fcbe8324034d5bf127b5aa69965` |
| Checkpoint state and arrays SHA | PASS | Match the supplied state and arrays digests |
| Stored checkpoint step | PASS | Manifest and state both report step 2000 |
| RNG state parse and PCG64 restore | PASS | `packing`, `velocity`, and `thermostat` states accepted by NumPy PCG64 |
| Candidate config parse and 2000→3000 boundary | PASS in mock harness | `start_step=2000`, `end_step=3000`, `max_added_steps=1000` |
| Mock suite | PASS, 6/6 | Includes checkpoint restore step, RNG consistency, step3001 refusal, wrong potential SHA rejection, wrong checkpoint SHA rejection, and boundary acceptance |
| Actual ASE adapter/core import | BLOCKED | ASE is not installed; required version is 3.29.0 |
| Actual ASE checkpoint restore dry-run | BLOCKED | Could not exercise the real ASE `Atoms`/`FixCom` path in this runtime |

The successful restoration unit tests traverse the production Python restoration code with a small ASE `Atoms`/`FixCom` test double. They are mock evidence, not a real ASE interface result. Archive parsing, checksum checks, state JSON, NPZ arrays, and NumPy PCG64 state restoration use the supplied checkpoint copy. The test suite did not call a production calculator or advance the underlying integrator.

## Changed code and hash pins

No existing production source file was changed. The candidate template records:

- adapter SHA: `4dfc0485f90d67cb56ab9ee8c8721ee2bd749aa13027d15d6c0f99c7d2a92fbb`
- core SHA: `adcc939fcb83b3ab3e4ed6e121e4ce3615bef73e1551c78a09f347d95985e3af`
- archive backend SHA: `0a8d8b71ee938e768f47bc9e819c5953f7f7bb0aa6415f248916b4a300a3928f`
- potential SHA pin: `958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447`

## Current blockers

1. ASE 3.29.0 is absent, so actual adapter/core imports and the requested real-interface dry-run cannot pass here.
2. A C-signed release config does not exist. The generated candidate template intentionally remains unsigned and unauthorized.

## C signature readiness

**Not ready to request C signing yet.** First rerun the import and checkpoint restore dry-run under the approved environment with ASE 3.29.0, then attach its real-interface PASS report. This report does not authorize execution and does not generate an execution-authorized config.
