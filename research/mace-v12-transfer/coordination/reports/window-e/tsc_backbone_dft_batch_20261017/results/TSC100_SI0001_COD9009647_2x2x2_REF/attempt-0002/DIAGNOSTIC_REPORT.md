# Si-terminated 100-atom reference SCF diagnostic

## Outcome

`TSC100_SI0001_COD9009647_2x2x2_REF`, attempt `attempt-0002`, stopped safely before the five-hour walltime limit after 37 SCF iterations. The SCF did not converge. No energy, force, or stress labels were emitted or accepted.

The supervisor requested a safe stop with about 600 seconds reserved for saving the latest complete checkpoint. GPAW exited with code 75 after the worker preserved that checkpoint. The stop reason is `SAFE_WALLTIME_STOP_NO_CONVERGENCE`; this was not an out-of-memory stop. The cgroup memory peak observed by the supervisor was 21,756,133,376 bytes under a 34,359,738,368-byte limit, and no cgroup OOM event occurred.

## Calculation and parameter record

The run used GPAW 26.7.0, ASE 3.29.0, four MPI ranks, and one OpenMP thread per rank. It retained the fixed input geometry and physical settings: PBE/PAW, 500 eV plane-wave cutoff, 6×6×1 Monkhorst–Pack mesh, 0.10 eV Fermi–Dirac smearing, non-spin-polarized calculation, and 1e-5 energy, density, and eigenstate convergence thresholds.

The only numerical mixing adjustment from attempt-0001 was MSR1 beta 0.05 → 0.02, chosen to damp the oscillatory SCF behavior recorded in the earlier attempt. The configuration diff is in `CONFIG_DIFF.json`. The full iteration log is in `gpaw.log`.

## Checkpoint and integrity

The GPAW 26.7.0 full-wavefunction checkpoint write/read path passed an MPI-4 preflight before the target run. During the run, full-mode single-precision checkpoints were SHA-256 verified after writing at steps 1, 11, 21, and 31. At safe stop, the worker wrote and verified the latest checkpoint at step 37, replacing the older checkpoint only after the new file had been hashed.

- Local checkpoint: `checkpoints/latest.gpw`
- SCF boundary: step 37
- Size: 3,113,760,276 bytes
- SHA-256: `2edbcb7b3803105c5157defcefaab6d63ddae179c73be7fa6463a466413c5b9a`
- Checkpoint mode: `all`, single precision, wavefunctions included

The checkpoint binary remains local and is not included in the Git archive. `CHECKPOINT_STATUS.json` records its hash. `SHA256SUMS.txt` covers the archived small files, including this report, and excludes the multi-gigabyte checkpoint itself. The last checkpoint was not reread while SCF was active; the write/read interface was validated by the preflight.

## Reproduction evidence

`RUN_MANIFEST.json`, `RUN_STATUS.json`, `SCF_PROGRESS.json`, `SCF_STOPPED.json`, `SAFE_STOP_REQUEST.json`, `CHECKPOINT_STATUS.json`, `CONFIG_DIFF.json`, `TARGET_PREFLIGHT.json`, and the GPAW iteration log preserve the run configuration, stop reason, progress, and checkpoint evidence. The attempt-0002 small-file archive was pushed in commit `f98053792d1f5c1241939cc3ce8a3afe0e285b0b` and verified after a fresh fetch.

No automatic resume or second parameter set was started. The result remains an incomplete SCF diagnostic, not a DFT training label.
