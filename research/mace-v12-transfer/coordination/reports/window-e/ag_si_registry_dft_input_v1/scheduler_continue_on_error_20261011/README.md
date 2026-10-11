# Ag–Si registry screening v2: continue-on-error scheduler

This package only prepares and tests the scheduler. No DFT was launched while preparing this policy. The active R1 attempt was stopped when the instruction changed to scheduler-only; its partial SCF/stdout/stderr logs remain in the separate `phase1_6x6x1_20261011_exec03/R1/` directory. R2–R4 were not started.

`SCHEDULER_POLICY.json` binds the unchanged C configuration and authorization SHA values. `run_registry_batch.py` defaults to dry-run; a future launch requires an explicit `--execute` flag and a fresh work package with R1–R4 copies of the authorized `input.py`, structure, and per-registry C release record. The script verifies those files and the runtime before starting. It runs one registry at a time with 4 MPI ranks and OMP=1.

For each registry, stdout, stderr, exit code, SCF status, and run status are written separately. If SCF reaches the frozen 160-iteration limit without convergence, GPAW raises `KohnShamConvergenceError`; the job is marked `FAILED_SCF`, the complete `results/coarse/gpaw.log` and any partial files are preserved, and dispatch continues to the next registry. Other per-registry process failures are recorded and also do not stop later jobs.

The scheduler never changes the input, structure, geometry, k-point mesh, smearing, mixer, cutoff, MPI/OMP allocation, or role. It does not automatically start 8×8×1. After all four jobs, `RANKING_6x6.json` contains `SUCCESS` and `FAILED` lists with failure reasons. Fewer than two successes yields `INSUFFICIENT_FOR_8X8`; two or more yields `READY_FOR_8X8_VALIDATION`. A separate human decision is required for any convergence repair.

Non-DFT checks:

- `python run_registry_batch.py --self-test` exercises continue-on-error and the ranking threshold using synthetic outcomes only.
- `python run_registry_batch.py --work-root <fresh-package>` validates a package, runs a four-rank filesystem/MPI smoke test, prints planned commands, and does not execute GPAW inputs.

Example future explicit invocation (not run during this task):

```sh
python run_registry_batch.py --work-root /path/to/fresh-package --execute
```
