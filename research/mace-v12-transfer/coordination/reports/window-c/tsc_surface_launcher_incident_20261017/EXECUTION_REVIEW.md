# E surface DFT launch incident review

Review date: 2026-10-17 UTC

## Decision

`C_EXECUTION_REVIEW = HOLD_FOR_LAUNCHER_REPAIR`

This is an execution-wrapper defect, not a failed SCF calculation. The first authorized task did not reach the GPAW worker and produced no DFT labels. Preserve the archived failed attempt; do not classify it as `FAILED_SCF`.

## Evidence reviewed

- E branch: `cloud-e-b1-dft`
- E commit: `a8c354fb3cb11efaec033782e716aabdc438425c`
- Task: `TSC100_TI0001_COD9009647_2x2x2_REF`
- Frozen input SHA-256: `b1948a44b89efe17ee7c2891b984167982067feeffa774b7f6fa77d9e9158d86`
- Config SHA-256: `e9a4e27d23c4646b873e77de7acbb79844422c8708c9f24d209827f1849312d2`
- Execution-lock SHA-256: `678d1d287482e82039fc59e0c6961a92069a817f2c856c8745c82cf58bbd5863`
- Worker SHA-256: `95f8661d9f29b607ae96dc421b81a3694972c0becca9dc3862c2ba2e6a9cd408`
- Supervisor SHA-256: `42dbe008db8145daec13d874a12efe104f2000bef4f3a386da5e819beb7553d5`
- Archived `RUN_STATUS.json` SHA-256: `6327b9aac9c4ffe0a1177e5040ae5599ce7222436afd0390e392b7d38e971dae`
- Archived `stderr.log` SHA-256: `dfbe9a9b7e8ea4d864e713b513e7c6966e107f14241c5e89c2b9dc82e006f9f6`

The archived status records exit code 2 after 2.0009 seconds, no checkpoint, no result, and no automatic retry. `stderr.log` contains the GPAW CLI usage error: it rejects worker arguments such as `--task-id`, `--input`, and `--run-dir`. The supervisor launches `gpaw python <worker> ...`; those worker arguments are being parsed by the GPAW command-line layer. The worker therefore did not construct a `GPAW` calculator or enter SCF. There is no `RESULT.json` or GPAW calculation log.

The failed run manifest confirms the execution preflight had passed before the output directory was created: GPAW 26.7.0, ASE 3.29.0, gpaw-data 1.2.1, 4 MPI ranks, 1 OMP thread, 4-core quota, 32 GiB cgroup memory limit, and 30,155,169,792 bytes initially free on `/workspace`. This immediate CLI failure does not measure the actual DFT peak memory or runtime.

The existing four-input scientific review remains unchanged: fixed geometries, PBE/PAW, PW 500 eV, FD 0.10 eV, standard MP `(6,6,1)`, MSR1, default Davidson `niter=2`, SCF `maxiter=160`, force-label candidates only. No geometry, role, or scientific parameter change is needed.

## Required repair before another launch

1. In the supervisor, invoke the worker directly with the same interpreter used by the successful MPI probe:

   `mpiexec -n 4 <sys.executable> <run_surface_singlepoint.py> <worker arguments>`

   Do not route the worker arguments through the `gpaw` CLI. The existing preflight already checks MPI using `mpiexec -n 4 sys.executable -c ...`.

2. Preserve the failed attempt at `results/TSC100_TI0001_COD9009647_2x2x2_REF/`. The supervisor currently requires that exact path to be new, so rerunning with the same path will fail preflight. Add an explicit, validated attempt ID and use a new leaf such as `results/TSC100_TI0001_COD9009647_2x2x2_REF/attempt-0002/`. Check that leaf for nonexistence; never overwrite or move the archived attempt.

3. Update the execution lock's `supervisor_sha256`, package `SHA256SUMS.txt`, and command instructions to bind the repaired supervisor and new lock SHA. Keep the worker, config, input structures, PAW hashes, and all scientific settings unchanged. Keep automatic retry disabled.

4. Before DFT, run a non-scientific MPI CLI smoke test (`mpiexec -n 4 <sys.executable> <worker> --help`) and confirm success with no result/checkpoint/calculation output. Then run the task's fresh resource/hash preflight for `attempt-0002`. A passing smoke test alone is not a release; C must review the new committed supervisor/lock hashes before the calculation command is considered valid.

After C verifies the repaired commit and bindings, E may run only the listed surface tasks sequentially, one at a time, with the existing per-task resource and fresh-fetch archive gates. No other DFT, fitting, MD, or training is authorized by this record.
