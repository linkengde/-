# Ti3SiC2 (10-11) reference restart API diagnostic

Date: 2026-10-23 UTC
Task: `TSC69_GEOMHOLD_PYR_10m11_COD9009647_TITI_REF`

## State before restart

The fixed-geometry PBE/PW500, standard Monkhorst-Pack 6x6x1 run reached SCF
iteration 20 without convergence. The energy and residuals remained strongly
oscillatory; no energy, force, or stress label is valid. A versioned checkpoint
was committed at iteration 20:

- File: `results/TSC69_GEOMHOLD_PYR_10m11_COD9009647_TITI_REF/checkpoint_iter_0020.gpw`
- Size: 4,583,533,942 bytes
- SHA-256: `6a2d52c2fc3eabf7e60aa39680a5828914e0568422137c8230573e2a59f9098d`
- Input SHA-256: `4b84cfc8f76827c7f81d9401f60e0c161355d61adc57cfe9ba867a78c780f473`
- Geometry SHA-256: `fd3351b1ba798ea9e0567013740fba8651d801fa8b62199b89da2c4d5ef37669`
- Configuration SHA-256: `15ddc8f596281298656a5a21df67b2915f0aebd5ba3bb18316fd2935ee8729c2`

An independent four-rank, read-only reload verified 69 atoms and iteration 20.
Rehashing the file after a failed restart attempt returned the same SHA.

## Resume failure and correction

The first resume attempt exited before starting an SCF iteration. GPAW 26.7.0's
new ASE frontend rejects `ASECalculator.set(maxiter=...)`; the runner had used
that unsupported call. No checkpoint or input was changed, and there was no
additional DFT step.

The runner now sets the remaining cap through the restored GPAW SCF loop:
`calc.dft.scf_loop.maxiter = locked_maxiter - checkpoint_iteration`. A
four-rank, read-only API probe loaded the exact checkpoint, confirmed its stored
limit was 160, assigned the remaining 140 iterations, and exited without
starting SCF. Therefore the cumulative ceiling remains 160; scientific
parameters and geometry remain unchanged. The runner and execution-lock hashes
were updated together.

## Current decision

The reference is still **not converged** and supplies no DFT labels. Its latest
checkpoint remains valid and locally available for an unchanged-parameter
resume. Disk free space at the checkpoint was 11,995,377,664 bytes, enough for
the next atomic checkpoint replacement under the runner's 4 GiB reserve plus
1 GiB transient-space guard. No other active GPAW calculation was present.

The response structure is blocked until the reference has a validated,
converged result. No TRAIN/DEV/TEST role is assigned: this is a same-COD-source
facet diagnostic, not an independent-source blind family.

The 4.58 GB GPW checkpoint is not included in normal Git objects. The repository
has no existing LFS-tracked checkpoint path, and this file exceeds common
per-object Git LFS limits. The checkpoint remains in the current workspace;
the committed manifests identify it by size and SHA but do not claim it was
remotely archived.
