# H P2 32-atom Γ-centered 1×8×6 fixed-geometry audit

P2 was run after the user directed continuation of the same 32-atom k-point ladder. Input SHA-256 `575686f5c9749c23c2258e7491633cd265367bad00432d0b9078feb89fe0f6b8` and all geometry bytes match P0 and P1. Only the Γ-centered mesh changed. Settings: PBE, PW 500 eV, Fermi-Dirac 0.10 eV, spin-unpolarized, four MPI ranks, GPAW 26.7.0, ASE 3.29.0, identical four-element PAW hashes. No structure relaxation or other calculation ran.

## P2 SCF result

- SCF converged: **YES**, in **51** iterations (`Converged in 51 steps`).
- Native energy: **-207.158494776088 eV**.
- Force-consistent free energy: **-207.385349300085 eV**.
- Maximum force: **3.305340665 eV/Å**.
- Wall time: **3657.905 s**; peak sampled process-tree RSS **3541258240 bytes** (shared pages may be counted more than once).

## P1 (1×6×4) ↔ P2 (1×8×6) gate

| Metric | P2 − P1 | Gate | Result |
|---|---:|---:|---|
| Native total energy, absolute difference per atom | 0.218012 meV/atom | ≤ 2 meV/atom | **PASS** |
| Force-consistent free energy, absolute difference per atom | 0.281323 meV/atom | ≤ 2 meV/atom | **PASS** |
| All-atom force-vector RMS | 0.012524772 eV/Å | ≤ 0.01 eV/Å | **FAIL** |
| Maximum per-atom force-vector difference | 0.022233895 eV/Å at atom 10 | report | — |

**Adjacent-mesh gate: FAIL.** The two energy criteria pass; the force-vector RMS criterion fails. Keep both results as independently SCF-converged labels, but do not claim full k-point convergence for this geometry. As context, P0 (1×4×2) to P2 (1×8×6) still differs by 9.226814 meV/atom native energy, 8.668152 meV/atom free energy, and 0.039042789 eV/Å force RMS.

This task is complete. No subsequent mesh, other structure, MD, or training job was started.
