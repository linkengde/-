# H acceptance report: TSC48 segmented runner v2

Audit date: 2026-10-16 UTC
Reviewed source commit: `779f9b2e16ec33421e52f62a27ff09b72e8b8107`
Runner SHA-256: `c177981e313ed95605c9e10896a2c8b6967c0df0821d99861e2bf36a19c9ac01`

## Required status

```text
V2_STATIC_REVIEW=FAIL
V2_OFFLINE_TESTS=PASS
MPI_DECISION_SYNC=PASS
REAL_GPAW_RESTART_VALIDATED=NO
RESOURCE_READY=NO
DFT_EXECUTION_READY=NO
DFT_STARTED=NO
```

The MPI stop/disk decision fix works in static review and a four-rank GPAW communicator exercise. V2's offline suite passes. Static review fails because the cumulative elapsed counter drops checkpoint-write time at every segment boundary, and its structured SCF samples do not extract the installed GPAW context's energy/convergence fields. The per-segment GPAW text logs are referenced and remain the available source for the full convergence trend. No scientific calculation was run.

## Hashes and offline checks

V2 was fetched from `chatgpt-tsc48-segmented-runner-v2-20261016` at the exact commit above. The local worktree was clean before audit.

| File | SHA-256 | Result |
| --- | --- | --- |
| `tsc48_segmented.py` | `c177981e313ed95605c9e10896a2c8b6967c0df0821d99861e2bf36a19c9ac01` | pinned v2 source |
| `TSC48_static_perturb_std0p025_seed2026101603_UNLABELED.extxyz` | `da3820c477023e0bb423f87550ed82ae0c2a571f53d294c055f78d7d052d5e8b` | approved input match |
| `calculation_config_preflight.json` | `ca6f365e8b0b8ec847e329d0312ebb51381b21e996050949494f2d3e0669dfe7` | approved config match |

Checks with `/workspace/.venvs/gpaw-mpi/bin/python`:

- `python -m py_compile .../tsc48_segmented.py` — passed.
- `python -m unittest discover -s .../TSC48_SEGMENTED_RUNNER_20261016 -p 'test_offline.py' -v` — 5/5 passed: dual-slot fake-file fallback/tamper, fake root-only decision broadcast, append-only ledger linkage, fail-closed resource clearance, and atomic JSON metadata.
- Runner `verify` — `OFFLINE_INPUT_CONFIG_PASS` with the exact input, config, and v2 runner SHA values above. It does not load GPAW.

An additional no-DFT MPI check invoked `shared_decision()` on four real GPAW MPI ranks with deliberately different rank-local wall-clock and free-space inputs. All ranks received identical wall-stop `(True, True, False)` and disk-block `(True, False, True)` decisions. This validates the helper's broadcast behavior, not the complete SCF callback or checkpoint path.

## Static review

- `shared_decision()` uses only rank 0's stop-file, elapsed-time, and free-space values, then broadcasts one fixed-size `int32` array. `collective_stage_clear()` also broadcasts rank 0's staged-file check. The SCF hook consumes the same broadcast decisions on every rank. This removes the V1 rank-local wall-clock branch that could send only some ranks into `calc.write()`.
- `save_checkpoint()` calls `calc.write(..., mode='all')` on all ranks and calls `world.abort(9)` after a Python-visible writer exception. GPAW 26.7.0 maps `mode='all'` to full wavefunction inclusion; its ULM writer writes on rank 0 while other ranks use dummy writers in the collective GPAW write path. The installed GPAW communicator exposes `abort`. No injected failure or actual GPW write was run, so hangs inside MPI/I/O and actual abort behavior remain untested.
- The A/B slots retain the alternate valid checkpoint while the target slot's `.writing.gpw` is staged. Size and SHA metadata select valid slots, and incomplete staged files are not restart candidates. The tests use fake bytes only.
- `ACTIVE_RUN.json` uses exclusive file creation to block a second runner sharing the work directory. A crash or unhandled exception intentionally leaves the lock, requiring independent process inspection before manual recovery. This is fail-closed; lock recovery itself was not exercised.
- The append-only ledger records source and output checkpoint SHA values, segment number, logged-step counter, elapsed values, and log path. However, `save_checkpoint()` stores `prior_seconds + elapsed` before the checkpoint write. `SEGMENT_PAUSED` later records a larger elapsed value after the write, but `latest_verified_counters()` restores the next segment from checkpoint metadata, not that pause event. With the documented 20-iteration segment, each boundary checkpoint's write time is omitted from the next cumulative total.
- `scf_context_sample()` probes `energy`, `e_total`, `e_error`, `density_error`, and `eigenstates_error`. GPAW 26.7.0's `SCFContext` exposes `ham.e_total_extrapolated`, `dens.error`, `eig_error`, and `wfs.eigensolver.error` instead. Its structured samples therefore contain iteration/time but no numeric convergence trend. The runner records the GPAW log path, whose text contains the trend, but does not record a log SHA-256.
- The pinned input, config, and science parameters are unchanged: fixed geometry, PBE/PW500, MP 10×10×4, Fermi–Dirac 0.10 eV, PPCG `niter=2`, Mixer(0.05, 8, 100), original convergence thresholds, and four ranks.

## Resource and restart gates

No `.gpw` file or site-specific clearance record was found. The only clearance JSON is the template, whose readiness fields are all false. Current free space is 28,827,013,120 bytes on `/workspace`'s overlay; `/workspace/scratch` is on the same overlay. `/tmp` is tmpfs with 17,297,383,424 bytes free. No scheduler or external wall-clock limit evidence was available; `ulimit -t unlimited` is not a platform wall-clock guarantee.

The two retained full-wavefunction slots plus a third staging file are the correct transient peak to budget. If the earlier approximate 11.5 GB per checkpoint estimate applies, three files require about 34.5 GB, a shortfall of about 5.67 GB against current free space before scratch. That size remains an estimate; no actual checkpoint write measured it. Therefore storage capacity and cross-session persistence remain unverified and resource clearance cannot be granted.

No history was deleted or changed. No scientific input or model was modified. No DFT, MD, relaxation, or training process was started. `DATA_ROLE` remains `HOLD`.

## Required follow-up before execution

Clarify whether `cumulative_elapsed_seconds` is intended to represent total wall time or SCF-active time. If it is total wall time, carry forward the post-checkpoint `SEGMENT_PAUSED` elapsed value (or equivalent durable checkpoint-write timing). Map SCF samples to the actual GPAW 26.7.0 context fields or explicitly rely on retained and integrity-recorded GPAW logs. Then obtain a real persistent work volume, measured three-file checkpoint plus scratch budget, external wall-clock and write-reserve evidence, and separately authorized real GPAW checkpoint/restart validation. Keep `RESOURCE_READY=NO` and `DFT_EXECUTION_READY=NO` until those gates pass.
