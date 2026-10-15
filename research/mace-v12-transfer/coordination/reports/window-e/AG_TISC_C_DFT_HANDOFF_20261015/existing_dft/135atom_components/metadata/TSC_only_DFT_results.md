# Ti₃SiC₂-only fixed-geometry DFT result

**Status: CONVERGED; recovered from the completed GPAW checkpoint.** This is the Ti₅₇Si₁₉C₃₈ subset of the optimized 135-atom interface, with its original positions, cell and PBC. No coordinate or cell relaxation was performed.

## Binding and method

- Full source structure: `source_optimized_structure.extxyz`, SHA-256 `eda2586c63c6c3c793ffce99417e7c04569941b0a2a916fda1b5402b7ba52d9f`.
- Subset input: `TSC_only/input_subset.extxyz`, SHA-256 `260ae9b736cb95a641b4d1eae9ffbec5711f2ffe895ed1f0acc7e25ff637fdf1`.
- Atom IDs map to source IDs 1–114 in order; composition Ti₅₇Si₁₉C₃₈, 114 atoms, neutral charge, 912 PAW valence electrons.
- Cell (Å): `[[13.38617865561341,0,0],[6.693089327806705,11.592770775358236,0],[0,0,37.9451596]]`; PBC `[true,true,true]`.
- GPAW 26.7.0, ASE 3.29.0, gpaw-data 1.2.1; PBE/PAW, plane-wave cutoff 500 eV, Γ 1×1×1, Fermi–Dirac width 0.10 eV, non-spin-polarized, neutral.
- SCF thresholds: energy, density and eigenstates each `1e-5`; maximum iterations 333. GPAW log records PPCG (`niter=2`) and its default MSR1 mixer. Execution used 4 MPI ranks and OMP=1.
- Ag, Ti, Si and C setup hashes match the E baseline: `bff00ac79a7fa6e8c3a71c8b34d3a9c488e6bf1c00b3e9c4c91c231551145e97`, `bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f`, `2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0`, `fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983`.

## Results

GPAW reports **SCF converged in 44 steps**, with an SCF-loop duration of 2682.108 s. Energies are not directly comparable to the full or Ag-only systems because their compositions differ.

| Quantity | Result |
|---|---:|
| Total energy, extrapolated (`force_consistent=False`) | −872.877542503438 eV |
| Force-consistent free energy, FD 0.10 eV | −874.741340407191 eV |
| Force vectors | 114, finite; eV/Å |
| Maximum force norm | 0.975716679118 eV/Å |
| Stress (eV/Å³) | `[[0.002733689774, 5.24870e-7, −1.05687e-7], [5.24870e-7, 0.002733223502, −1.07035e-7], [−1.05687e-7, −1.07035e-7, −0.002161869100]]` |
| Occupations | finite; sum 912.0, min `1.0131e-10`, max 2.0; 95 fractional bands; Fermi level −1.398515 eV |

The large forces are on the unchanged interface-derived frozen geometry after removing Ag. They are diagnostic forces, not a claim that this isolated slab is relaxed.

## Recovery and integrity

The first runner completed the SCF and wrote `TSC_only/partial.gpw`, then its collector failed by reading a nonexistent `ASECalculator.scf` attribute. No SCF rerun was made. `recover_tsc_checkpoint.py` verified the convergence marker in the GPAW log, all 114 positions, cell, PBC, element order, finite forces/stress and the neutral electron/occupation sum before writing `RECOVERED_RESULT.json` and `recovered_forces.csv`. The filename `partial.gpw` is historical; its checkpoint contains the converged result.

- Checkpoint SHA-256: `fed4ea635d1563f2d1f963f16f91c482544233bbd0d9945835b3b610d54afdbb`.
- Recovered force CSV SHA-256: `299d3cad8e984bbb6a0c084322e763104c0da658d442c0e0862929b96799c1d5`.
- GPAW log SHA-256: `a9c8f815372c511702f5e9447ebc51a68972a796f0341f990c72352ed714fc54`.

