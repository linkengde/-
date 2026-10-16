# H acceptance report: TSC48 segmented runner

Audit date: 2026-10-16 UTC
Runner source commit: `7e77a8e2ee05e4ec1a9bdb2b0afe325c616cc003`
Runner SHA-256: `0d9b97d2ed138985fc75e292485288ef7755e63399e2e4ecdaaeb30804c47e8f`

## Result

The pinned runner and offline acceptance checks pass. Static review confirms that GPAW 26.7.0 accepts the configured `scf_step` hook and that the runner enters `calc.write(..., mode='all')` on every GPAW rank. The three tests use fake checkpoint bytes. No real GPAW checkpoint was written or resumed, so checkpoint restart behavior, persistent storage, and resource clearance remain unverified.

```text
RUNNER_SHA_VERIFIED=YES
PYTHON_SYNTAX=PASS
OFFLINE_TESTS=3_PASS (fake checkpoint bytes only)
INPUT_CONFIG_VERIFY=PASS
GPAW_SCF_HOOK_STATIC_REVIEW=SUPPORTED
FOUR_RANK_GPAW_MPI_IMPORT=PASS
COLLECTIVE_GPW_WRITE=STATIC_PATTERN_REVIEWED / NOT EXECUTED
DUAL_SLOT_FAKE_FILE_TEST=PASS
REAL_GPAW_CHECKPOINT_WRITE=UNVERIFIED
REAL_GPAW_RESTART_AND_SCF_CONTINUATION=NO
PERSISTENT_STORAGE_READY=NO
RESOURCE_READY=NO
DFT_EXECUTION_READY=NO
DFT_STARTED=NO
```

## Reproducible checks

The files were read from branch `chatgpt-tsc48-segmented-runner-20261016` at the exact source commit above. The local project worktree was clean before audit. The runner, input, and H config hashes match the README and approved lock values:

| File | SHA-256 | Result |
| --- | --- | --- |
| `tsc48_segmented.py` | `0d9b97d2ed138985fc75e292485288ef7755e63399e2e4ecdaaeb30804c47e8f` | match |
| `TSC48_static_perturb_std0p025_seed2026101603_UNLABELED.extxyz` | `da3820c477023e0bb423f87550ed82ae0c2a571f53d294c055f78d7d052d5e8b` | match |
| `calculation_config_preflight.json` | `ca6f365e8b0b8ec847e329d0312ebb51381b21e996050949494f2d3e0669dfe7` | match |

Executed with `/workspace/.venvs/gpaw-mpi/bin/python`:

- `python -m py_compile .../tsc48_segmented.py` — passed.
- `python -m unittest discover -s .../TSC48_SEGMENTED_RUNNER_20261016 -p 'test_offline.py' -v` — 3/3 passed: fake dual-slot fallback/tamper, fail-closed resource check, atomic JSON metadata.
- `python .../tsc48_segmented.py verify --input ... --h-config ... --workdir /tmp/tsc48-offline-verify-no-workdir` — `OFFLINE_INPUT_CONFIG_PASS`; it did not import GPAW or create the work directory.

Runtime inspection after sourcing `/workspace/.setup/runtime.sh` found Python 3.12.14, GPAW 26.7.0, ASE 3.29.0, and the expected `GPAW_SETUP_PATH`. Ti, Si, and C PAW SHA-256 values matched the runner pins. `mpiexec -n 4` launched four processes and each reported its GPAW rank and version. This validates MPI process initialization only; no GPAW calculator or checkpoint write was run.

## GPAW 26.7.0 source review

- `gpaw/new/ase_interface.py` permits `scf_step` and invokes it after each yielded SCF context. The runner's callback name and argument shape are consistent with this installed source.
- `ASECalculator.write()` maps `mode='all'` to `GPWFlags(include_wfs=True)`. `gpaw/new/gpw.py` uses a real ULM writer on rank 0 and dummy writers on other ranks, then completes a communicator barrier. The runner correctly calls `calc.write()` from all ranks rather than rank 0 alone.
- `restart(filename)` constructs a GPAW calculator from the GPW file. A full-mode checkpoint includes PW wavefunction coefficients, so the design is a wavefunction-based warm restart. However, this does not preserve the complete in-memory SCF loop: installed `SCFLoop.iterate()` resets the mixer and starts its iteration counter at 1 on each invocation. It is not a bitwise continuation of mixer history. No real checkpoint has been loaded to validate convergence behavior after such a restart.
- The A/B slot selection, staged `.writing.gpw`, metadata size/SHA validation, and preference for the highest valid sequence are consistent in static inspection. The offline tests prove these rules only for fake byte files, not for GPAW ULM/GPW data or MPI I/O failures.

## Findings requiring correction or evidence

1. **MPI stop-decision synchronization:** each rank independently evaluates elapsed wall time and `STOP_REQUESTED` in `on_step()`. The iteration-count conditions are deterministic, but the wall-clock threshold can be crossed between ranks. If ranks disagree, only some may enter collective `calc.write()` and the job can hang. Before execution, make the stop/checkpoint decision collective (for example, rank-0 decision broadcast or an all-rank reduction) and test the resulting control flow without starting an unauthorized DFT.
2. **Collective write failure handling:** the checkpoint `calc.write()` is outside `root_action()`. A rank-0 filesystem error during the GPAW collective writer is not converted into a guaranteed all-rank failure; MPI hang behavior is untested. The staged file is excluded by slot verification, but that alone does not prove every rank exits safely after an I/O failure.
3. **Per-segment audit record is incomplete:** slot metadata records a sequence, checkpoint hash, and iteration; `segment_status.json` is overwritten and records only the latest status/sequence and local segment iteration count. The runner does not persist cumulative elapsed time, SCF convergence trend, or the exact prior checkpoint SHA used as the resume basis for each segment. Add an append-only segment ledger before treating the requested per-segment audit trail as satisfied.
4. **Resource clearance remains absent:** the delivered template has all readiness booleans false; no site-specific clearance JSON or GPW file was found. `/workspace` has 28,827,324,416 bytes free on an overlay filesystem; `/tmp` has 17,297,489,920 bytes free on tmpfs. `/workspace/scratch` is not proven to be separate durable storage. A dual-slot update can temporarily require the two retained checkpoints plus a third staged checkpoint, in addition to scratch. Actual checkpoint size and scratch peak are unmeasured.
5. **Persistence and external wall clock remain unverified:** no scheduler command or platform wall-clock variable was available; `ulimit -t unlimited` is not an external job deadline. No evidence demonstrates that `/workspace` survives a new session. The four-rank import test does not establish a persistent checkpoint path or write-time reserve.

No history was deleted or edited. No science parameters, input structures, or model files were changed. No DFT, MD, relaxation, or training process was started. `DATA_ROLE` remains `HOLD`.

## Release gate

Do not create a site clearance file or invoke `start`/`resume --execute` on the present evidence. The next review needs the synchronized MPI callback fix and offline validation, a real persistent work volume, measured full-checkpoint/staging and scratch budget, externally enforced wall-clock and checkpoint-write reserve evidence, and an independently authorized real GPAW restart validation. Until then, `RESOURCE_READY=NO` and `DFT_EXECUTION_READY=NO`.
