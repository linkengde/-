# H V3 acceptance request — TSC48 segmented GPAW runner

Parent V2 acceptance report:
- branch: audit-h-tsc48-segmented-acceptance-20261016
- commit: 115069e0e015395badbd2d77fa946ad6d1df127b
- result: V2_STATIC_REVIEW=FAIL; V2_OFFLINE_TESTS=PASS; MPI_DECISION_SYNC=PASS

V3 branch:
- branch: chatgpt-tsc48-segmented-runner-v3-20261016
- latest source must be read from this branch; do not use V1/V2 runner SHA.
- scientific input, approved H configuration and scientific parameters are unchanged.

V3 fixes:
1. Checkpoint cumulative elapsed time is sampled after GPW publication and propagated
   back to all ranks, so the next segment inherits elapsed time that includes the
   checkpoint operation. The manifest also records a separate measured
   checkpoint_write_elapsed_seconds covering the actual calc.write() interval.
2. SCF structured samples use the GPAW 26.7 SCFContext interface:
   ctx.niter, ctx.ham.e_total_extrapolated, ctx.dens.error,
   ctx.wfs.eigensolver.error, ctx.eig_error, ctx.wfs.nvalence.
   Units are converted consistently with GPAW convergence criteria.
3. V2 MPI rank-0 decision broadcast, dual-slot fallback and ACTIVE_RUN lock remain.

Required V3 offline acceptance:
- python -m py_compile tsc48_segmented.py
- python -m unittest discover -s TSC48_SEGMENTED_RUNNER_20261016 -p 'test_offline.py' -v
- expected: 7/7 tests PASS
- verify input/config SHA and runner SHA from this branch
- no GPAW calculation is required for this offline step

Required static checks:
- verify no V2 bad SCF field names remain
- verify checkpoint_write_elapsed_seconds measures only the GPW write interval
- verify cumulative_elapsed_seconds is measured after checkpoint publication
- verify all MPI ranks receive the same decision before collective write
- verify fallback slot remains intact until the new GPW has been published
- verify no scientific parameter changed

After V3 offline/static PASS, keep real execution blocked until the already-required
real 4-rank checkpoint write/readback, interruption behavior, persistent storage,
scratch peak and wallclock evidence are separately validated.

Do not start DFT, MD, AIMD, TTM or training in this acceptance step.
