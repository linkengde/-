# R2 Phase2 SCF recovery report

## Outcome

**Recovery status: NOT_RECOVERED.** No diagnostic run reached SCF convergence, so total energy, force-consistent free energy, atomic forces, and stress remain unavailable. Partial SCF energies are not reported as labels.

## Baseline and bindings

- Phase2 source commit: `0dd8628e4a91e88f31aba3fd7bbce1c8a40a2aa2`
- Original R2 status: `FAILED_SCF`, 160 steps, `KohnShamConvergenceError` at the locked 160-iteration limit.
- Original input SHA-256: `7c1327f403c94240e6181da990fe7b078f420cc32f92a68bee3282792b0bcc66`
- R2 structure SHA-256: `cb2e16376887f7cdeb05c7ddf4aacfacacbea0a000837a3888a0f3ea91785f70`
- Copied `input.py` and `structure.xyz` hashes match those original SHA values exactly. The original files were not modified.
- Fixed geometry, R2 registry, PBE/PAW, 8×8×1 k-points, 500 eV cutoff and 0.10 eV smearing were kept fixed. No relaxation, d0 scan or role assignment was performed.

## Environment preflight

- GPAW 26.7.0; ASE 3.29.0; gpaw-data 1.2.1; OpenMPI 5.0.7.
- MPI launch was verified through `mpiexec -n 4 gpaw python`; the diagnostic ran with 4 MPI ranks and `OMP_NUM_THREADS=1`. A plain `mpiexec python` probe correctly failed because GPAW MPI was disabled in that mode; the supported `gpaw python` launcher passed. One diagnostic-script startup then failed to import the sibling `input.py` before any SCF iteration; the module path was fixed, and the successful preflight and SCF attempts are recorded separately.
- Ag, Ti, Si and C PAW SHA-256 checks matched the Phase2 parameter lock.
- CPU quota: 4 cores; cgroup memory limit: 32 GiB. No other active GPAW or MPI calculation was present before the diagnostic.

## Diagnostic attempts

| Attempt | maxiter | Mixer (beta, nmaxold, weight) | Steps observed | Outcome |
|---|---:|---|---:|---|
| 1 | 240 | 0.10, 8, 100 | 13 | Stopped after persistent large oscillations; partial log retained. |
| 2 | 300 | 0.02, 8, 100 | 34 | Density residual and energy remained unstable; partial log retained. |
| 3 | 300 | 0.05, 3, 100 | 123 | Stopped while still oscillatory and unconverged; partial log retained. |

Only maxiter and mixer settings varied. All runs used independent output directories under this R2 diagnostic copy. No parameters from these tests replace the original Phase2 result.

## Final data status

- Total energy: `NOT_AVAILABLE`
- Force-consistent free energy: `NOT_AVAILABLE`
- Atomic forces: `NOT_AVAILABLE`
- Stress: `NOT_AVAILABLE`
- Final diagnostic SCF: not converged after 123 observed steps (300 maxiter setting); no final result JSON was emitted by GPAW.
- Machine-readable status: `R2_PHASE2_SCF_RECOVERY_RESULT.json` (SHA-256 `411395ccd1e478cd1f1456625c0ee4393c88bd58528fcd40b95addf9b31d090e`).

Full GPAW logs, launcher output, attempt status files, diagnostic runners and SHA lists are retained alongside this report. The original Phase2 R2 artifacts and R1 result remain unchanged.
