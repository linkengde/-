# A review of the V17 development-validation runner

- The first authorized AgC development-validation launch passed the package preflight but stopped at GPAW import before SCF. GPAW reported: `Please use \"gpaw python\" to run in parallel`. No calculator/SCF result or label was produced.
- Preserved the launcher traceback at `../development_validation_attempts/AgC_v17_dev_validation_01.launcher_attempt_01.log` (SHA256 `c8131e80b0fcb3bfbd7f701d07d76142cc8d553fa94ad59819163d3a5df72297`).
- Updated the guarded runner to execute `execution_preflight.py` with the selected GPAW Python environment, and launch the four-rank calculation as `mpirun ... <same-env>/gpaw python dft_single_point.py`. No input geometry, recipe, role, hash or DFT setting changed.
- Verified the repaired `gpaw python` launch through four MPI ranks (`gpaw.mpi.world.size == 4`) without starting a calculator. Full role/input authorization and geometry preflight passed again before any DFT restart.
- This is a runner correction, not a calculation result. The failed attempt remains documented; the label is still pending SCF and archive verification.
