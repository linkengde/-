# Ag-only fixed-geometry DFT result

**Status: CONVERGED.** This is the Ag₂₁ subset of the optimized 135-atom interface, with its original coordinates, cell and PBC. No coordinate or cell relaxation was performed.

## Binding and method

- Full source structure: `source_optimized_structure.extxyz`, SHA-256 `eda2586c63c6c3c793ffce99417e7c04569941b0a2a916fda1b5402b7ba52d9f`.
- Subset input: `Ag_only/input_subset.extxyz`, SHA-256 `e75ed95aa31ff7340adeb6a8c0644e117a9c8b22bf13515262b86e8872afc35d`.
- Atom IDs map to source IDs 115–135 in order; composition Ag₂₁, 21 atoms, neutral charge, 357 PAW valence electrons.
- Cell (Å): `[[13.38617865561341,0,0],[6.693089327806705,11.592770775358236,0],[0,0,37.9451596]]`; PBC `[true,true,true]`.
- GPAW 26.7.0, ASE 3.29.0, gpaw-data 1.2.1; PBE/PAW, plane-wave cutoff 500 eV, Γ 1×1×1, Fermi–Dirac width 0.10 eV, non-spin-polarized, neutral.
- SCF thresholds: energy, density and eigenstates each `1e-5`; maximum iterations 333. GPAW log records PPCG (`niter=2`) and default MSR1 mixer. Execution used 4 MPI ranks and OMP=1.
- Ag PBE setup SHA-256 matches the E baseline: `bff00ac79a7fa6e8c3a71c8b34d3a9c488e6bf1c00b3e9c4c91c231551145e97`.

## Results

GPAW reports **SCF converged in 13 steps**, with an SCF-loop duration of 461.217 s.

| Quantity | Result |
|---|---:|
| Total energy, extrapolated (`force_consistent=False`) | −44.714091632522 eV |
| Force-consistent free energy, FD 0.10 eV | −45.116114998658 eV |
| Force vectors | 21, finite; eV/Å |
| Maximum force norm | 0.372812057566 eV/Å |
| Stress (eV/Å³) | `[[0.003163175397, −5.83231e-8, 1.21492e-8], [−5.83231e-8, 0.003163272400, 4.62393e-9], [1.21492e-8, 4.62393e-9, 0.000146470451]]` |
| Runtime occupation check | sum 357.0, min `3.6039e-24`, max 2.0; 13 fractional bands; Fermi level −3.600775 eV |

Ag₂₁ has an odd neutral valence count. Under the retained non-spin-polarized setup, fractional Fermi–Dirac occupations are expected and the runtime occupation sum matches 357. This result does not test whether a spin-polarized Ag₂₁ state would be lower in energy; that is a limitation of this same-method component comparison.

## Collection and integrity

The first run converged and logged energy, free energy and forces but stopped before stress and checkpoint collection. Its log is preserved in `Ag_only_attempt01_collector_failure/`. A second run with the identical structure and numerical settings wrote `Ag_only_retry02/RESULT.json`, 21 forces and the checkpoint. The output was independently reread: checkpoint energy, free energy, forces, stress, coordinates, cell and PBC match the recorded result; input, force and checkpoint SHA checks pass.

The launcher ultimately returned a nonzero status after rank 0 had written the verified result: non-root ranks tried to read the rank-0-only text log and misclassified their empty log as SCF failure. This was postprocessing only; the common MPI calculation had already completed and the checkpoint is readable. The collector source has been corrected so convergence is based on GPAW's successful completed calculation and public iteration count, not rank-local log contents. No third SCF run was performed.

- Checkpoint SHA-256: `e4fcf89daa7bfca09c233fc30d492876c7c1e29a232d21c2f6d7787c4edb7a75`.
- Force CSV SHA-256: `1c1f3bdc0cb0ff9177bdf723bb52ddd2b40e8c80b5cc9d3fb5f4521f0c88dbb1`.
- GPAW log SHA-256: `f81af61fa36eb112909337a4a4fbecd643b2ece7abf8eff2094297143d633b66`.
- Result JSON SHA-256: `6079447f0bfd5e1a2c59ce7cd7f03ccb0a57b41bec241b4d67e8d55053fee6a0`.

