# R2 SCF recovery package v2

This package prepares one hash-bound R2 diagnostic point. It does not authorize or start DFT by itself. `R2_RECOVERY_CONFIG.json` remains `execution_authorized: false`, and `runner.py` defaults to a read-only preflight.

The recovery uses only the Phase2 collector lock at `registry_phase2_result_collector_20261012/PHASE2_PARAMETER_LOCK.json`, SHA-256 `d97de3fe5a26ab6529d5e6e8ceef15adc1f3a2f548f03b59f9f092361c2c0c78`. The runner checks the exact path and digest and rejects the older top-level lock. It also rechecks the frozen R2 input and structure SHA-256 values.

The only numerical recovery changes from the Phase2 lock are `PPCG(niter=4)` and the finite SCF limit from 160 to 200. GPAW 26.7.0's new ASE interface is required because its `scf_step` hook is called after each completed SCF iteration; its parameter mapping for PPCG is `{"name":"ppcg","niter":4}`. PBE, PAW hashes, 500 eV, 8×8×1, 0.10 eV Fermi–Dirac smearing, the original 1e-5 convergence values, Mixer 0.05/8/100, symmetry, Poisson settings, geometry, MPI=4 and OMP=1 remain pinned.

The runner writes a full-wavefunction GPW at every 40th completed SCF iteration. Each file is written to a temporary name, read back through GPAW's `read_gpw`, checked for atom/cell/PBC identity and wavefunctions, hashed, and paired with an atomic JSON sidecar binding that SHA-256 to the callback iteration. Checkpoint failure is broadcast after collective barriers; all ranks stop, prior verified checkpoints remain, and the runner does not retry. The sidecar iteration is diagnostic metadata: GPAW's GPW restart state contains density and wavefunctions, while the SCF iteration number is not an internal resume counter.

For orderly stopping, the runner writes a checkpoint and stops at an SCF boundary at 12,600 seconds, leaving a 1,800-second margin. `launch_r2.sh` sends TERM at 14,399 seconds and uses a one-second kill window to enforce the 14,400-second outer limit if an SCF iteration is stuck. In that fallback case the last verified checkpoint and captured logs remain; labels are rejected unless the run status and GPAW's `Converged in N steps` marker both confirm convergence.

After C issues a new authorization bound to the exact runner, collector, config, input, structure and lock digests, run with:

```bash
launch_r2.sh /absolute/path/REGISTRY_R2_RECOVERY_AUTH.json AUTH_SHA256 /workspace/e-r2-recovery-run
```

The work root must be new and outside the Git repository. Without a signed C authorization and matching SHA values, `--run` exits before creating it. A failed SCF or runtime error emits `FAILED_SCF` with no energy, force or stress labels. No automatic retry or continuation is implemented.

The repeatable checkpoint readback test uses one hydrogen atom and zero SCF iterations in a single process; the separate four-rank failure-injection test verifies collective stop behavior. An optional four-rank synthetic hydrogen readback completed once, but a final repeat spent over a minute in LCAO initialization before reaching the checkpoint hook and was stopped. That repeat is not counted as a pass. No R2 DFT was started in this preparation round.

The detailed C `REJECTED` review text was not present in the fetched C refs or E working tree. This package implements the explicit acceptance checklist supplied with this task and does not infer additional review comments.
