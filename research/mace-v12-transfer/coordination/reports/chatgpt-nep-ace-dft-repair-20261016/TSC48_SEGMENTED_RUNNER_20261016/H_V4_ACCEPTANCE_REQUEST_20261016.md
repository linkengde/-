# H V4 acceptance request — TSC48 segmented GPAW runner

Source branch: `chatgpt-tsc48-segmented-runner-v4-20261016`.
Scientific input and H-approved PBE/PW500/standard-MP 10x10x4 configuration are unchanged.
Latest runner code must be verified directly from this branch; V1/V2/V3 SHA are obsolete.

## Previous findings

H V3: 7 offline tests passed, but cumulative checkpoint accounting still
counted the post-write MPI barrier in `checkpoint_write_elapsed_seconds`;
manifest cumulative elapsed was sampled before subsequent metadata/ledger
operations; SCF sample did not actually consume `ctx.niter`.

## Surgical V4 modifications

1. Timestamp actual `calc.write(...,mode='all')` return **before** MPI barrier.
2. Publish full GPW safely with an initial manifest explicitly marked
   `checkpoint_finalized=false`, after GPW fsync and SHA-256 verification.
3. After publishing preliminary manifest, ledger and status, measure the
   cumulative elapsed time and atomically publish the FINAL manifest with
   `checkpoint_finalized=true` and `cumulative_elapsed_seconds`.
   Verified slot selection now rejects all unfinalized checkpoints.
4. SCF sampler reads `ctx.niter` into `scf_iteration_from_context` while
   preserving separate local hook iteration counter.
5. Still **seven** offline tests; updated to test a purposely delayed
   MPI barrier, differing `ctx.niter`/local count, and unfinalized-slot refusal.

## H V4 offline verification (NO scientific computations)

- Fresh checkout and SHA256 of runner, approved input, approved H config.
- `python -m py_compile tsc48_segmented.py`.
- `python -m unittest discover -s TSC48_SEGMENTED_RUNNER_20261016 -p test_offline.py -v`.
- Runner `verify --input ... --h-config ... --workdir ...`.
- Static inspect `save_checkpoint` timing order and finalized manifest gate.
- Confirm rank-0 broadcast and 4-rank collective write design unchanged.
- Check restart warm-start semantics: GPW wavefunctions, not a bitwise
  preservation of mixer state/SCF iteration sequence.

Return V4_STATIC_REVIEW, V4_OFFLINE_TESTS, CHECKPOINT_TIMING_ACCOUNTING,
SCF_CONTEXT_MAPPING, MPI_DECISION_SYNC and specific unpassed evidence.

**Do not launch DFT/MD/TTM/training or set `resource_ready=true`**.
Actual GPAW GPW write/read, SCF restart, persistent checkpoint storage,
scratch/staging peak and platform wallclock remain separately unverified.
On the previously reported 28.83-GB overlay volume the estimated three-GPW
staging peak (~34.5 GB before scratch) is a blocking risk, not proof of
resource readiness.
