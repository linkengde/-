# H P3 32-atom Γ-centered 1×10×8 fixed-geometry audit

This is one user-directed extension beyond the frozen queue's 1×8×6 rung. The exact frozen input SHA-256 is `575686f5c9749c23c2258e7491633cd265367bad00432d0b9078feb89fe0f6b8`; structure, atom order, cell, coordinates, and PBC were unchanged. Method remained PBE, PW 500 eV, Fermi-Dirac 0.10 eV, spin-unpolarized, GPAW 26.7.0, ASE 3.29.0, four MPI ranks, and the same Ag/Ti/Si/C PAW hashes. No relaxation or other scientific task ran.

## SCF result

- SCF converged: **YES**, in **53** iterations (`Converged in 53 steps`).
- Native energy: **-207.199171421666 eV**.
- Force-consistent free energy: **-207.413242146206 eV**.
- Maximum force: **3.304526569 eV/Å**.
- Wall time: **6438.995 s** (107.32 min).
- Peak sampled process-tree RSS: **4361879552 bytes** (shared pages may be counted more than once).

## Adjacent mesh gate: 1×8×6 → 1×10×8

| Metric | P3 − P2 | Gate | Result |
|---|---:|---:|---|
| Native total energy, absolute difference per atom | 1.271145 meV/atom | ≤ 2 meV/atom | **PASS** |
| Force-consistent free energy, absolute difference per atom | 0.871651 meV/atom | ≤ 2 meV/atom | **PASS** |
| All-atom force-vector RMS | 0.010746257 eV/Å | ≤ 0.01 eV/Å | **FAIL** |
| Maximum per-atom force-vector difference | 0.019446618 eV/Å at atom 13 (Ti) | report | — |

**Adjacent-mesh gate: FAIL.** Energy criteria pass; force-vector RMS remains above the threshold. This result is an independently SCF-converged label, but this adjacent pair does not pass the full numerical-resolution gate. The previous P1→P2 force RMS also failed (0.012524772 eV/Å). No next mesh or other job was started.
