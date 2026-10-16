# H P1 32-atom Γ-centered 1×6×4 fixed-geometry audit

Instruction: `NEXT_SINGLE_1x6x4_CLOUD_INSTRUCTION.md` at review commit `3774b4c23abcb7351507349232c2960d04c0b44f`. P0 reference commit: `ce501377a33a36b059621adefbf2c18ec9b6a2f6`; P0 result SHA-256: `0a6fe5669af67200beea1af442d078408e277c850b5e1464b8e9313b062999df`.

The exact frozen structure bytes, atom order, cell, coordinates, and PBC were reused (SHA-256 `575686f5c9749c23c2258e7491633cd265367bad00432d0b9078feb89fe0f6b8`). Only the mesh changed: Γ-centered 1×4×2 → Γ-centered 1×6×4. Method remained PBE, PW 500 eV, Fermi-Dirac 0.10 eV, spin unpolarized, GPAW 26.7.0, ASE 3.29.0, four MPI ranks, and identical Ag/Ti/Si/C PAW hashes. No relaxation or other scientific task was started.

## P1 SCF result

- SCF converged: **YES**, in **49** iterations (`Converged in 49 steps`).
- Native energy: **-207.165471144659 eV**.
- Force-consistent free energy: **-207.376346949665 eV**.
- Maximum force: **3.297094871 eV/Å**.
- Wall time: **2090.158 s**; peak sampled process-tree RSS **3025080320 bytes** (shared pages may be counted more than once).

## P0 ↔ P1 comparison

| Metric | P1 − P0 | Gate | Result |
|---|---:|---:|---|
| Native total energy, absolute difference per atom | 9.444825 meV/atom | ≤ 2 meV/atom | **FAIL** |
| Force-consistent free energy, absolute difference per atom | 8.386828 meV/atom | ≤ 2 meV/atom | **FAIL** |
| All-atom force-vector RMS | 0.037974370 eV/Å | ≤ 0.01 eV/Å | **FAIL** |
| Maximum per-atom force-vector difference | 0.077263348 eV/Å at atom 20 | report | — |

**K-point convergence gate: FAIL.** All three criteria fail. Retain both converged calculations as reproducible labels, but do not mark this structure as k-point converged or use this pair as a validated production reference. Stop here; 1×8×6 and other jobs remain unstarted.

Historical 1×4×2 centering remains unverified; it does not affect this new P0/P1 pair because both runs explicitly use Γ centering.
