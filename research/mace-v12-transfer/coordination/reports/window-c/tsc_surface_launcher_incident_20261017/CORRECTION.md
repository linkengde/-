# Correction to the E surface DFT launcher review

Date: 2026-10-17 UTC

The first review's proposed fix to bypass the GPAW CLI and invoke `sys.executable` directly was incorrect. The root cause is narrower: the `gpaw python` invocation omitted GPAW's required `--` separator before the worker script and its options.

The correct form is:

```text
mpiexec -n 4 <gpaw executable> python -- <run_surface_singlepoint.py> <worker arguments>
```

Evidence: in the installed GPAW 26.7.0 environment, `gpaw python --help` states that `--` forces remaining arguments to the target program. I ran the corrected four-rank invocation against the existing worker with only `--help`; all four ranks displayed the worker's argparse help and exited successfully. This smoke test imported the GPAW modules but did not construct a calculator, call energy/force evaluation, or start SCF.

The separate retry-path issue remains: the failed attempt already occupies `results/TSC100_TI0001_COD9009647_2x2x2_REF/`. Preserve it and use a fresh attempt-specific leaf such as `results/TSC100_TI0001_COD9009647_2x2x2_REF/attempt-0002/`. Update the supervisor hash in the execution lock and regenerate package checksums. No scientific parameters, structure, input, worker, or role should change.

No DFT was started by this correction test. The existing execution lock is stale once the supervisor is edited; C must review the updated commit and hashes before E starts the calculation.
