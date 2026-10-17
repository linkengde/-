# E repair instruction: surface DFT worker launcher

Source branch head reviewed: `a8c354fb3cb11efaec033782e716aabdc438425c`.

The first task exited in the GPAW command parser before SCF. Preserve its archived directory unchanged. Do not call this an SCF failure and do not alter any input structure or scientific setting.

For the single task `TSC100_TI0001_COD9009647_2x2x2_REF`:

1. Change only the supervisor launch vector by inserting GPAW's required `--` separator: `mpiexec -n 4 <gpaw executable> python -- <run_surface_singlepoint.py> <worker arguments>`. Keep the existing `gpaw python` launcher.
2. Add a validated attempt selector. Leave the existing failed archive untouched and target a new path ending in `attempt-0002`. The output-path preflight must check that exact new path, not the already-existing task parent.
3. Refresh the supervisor SHA in the execution lock, then regenerate the package SHA list and command file. Do not alter the input/config/worker/PAW hashes or DFT parameters. Keep no-retry behavior.
4. Test only the CLI path with four MPI ranks using `gpaw python -- <worker> --help`; verify success and that no `RESULT.json`, GPAW log, or checkpoint is created. Run the updated `--preflight-only` against `attempt-0002`.
5. Commit and push the repair to `cloud-e-b1-dft`, then report the commit SHA, supervisor SHA, lock SHA, package checksum result, and smoke/preflight results to C. Do not start the DFT until C confirms the updated bound hashes.

The calculation scope and parameters remain those in the already-reviewed four-task package. This repair grants no additional task or retry beyond a single new attempt for the failed pre-SCF launch.
