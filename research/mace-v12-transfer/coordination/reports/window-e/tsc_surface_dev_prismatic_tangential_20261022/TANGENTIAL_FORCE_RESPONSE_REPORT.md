# Tangential Ti-displacement response diagnostic

Status: **DFT and frozen-potential run-0 checks passed**. One 120-atom fixed-geometry PBE/PAW point was constructed from the existing COD 9009647 (11-20) Ti3SiC2 slab. This remains a same-parent diagnostic; it does **not** create an independent DEV/TEST family or change any data role.

## Structure and calculation

The pair of terminal Ti source IDs 4 and 113 was displaced by +0.02 and -0.02 Angstrom along the first in-plane Cartesian direction, respectively. Atom count/composition are Ti60Si20C40; cell, PBC, atom sequence, and all other coordinates match the frozen source reference. Input SHA-256: `f3ff64c30213d41bdac54d667d1a99a0e7dcd5acad2f8c93a1f0336cb9450b8d`; geometry SHA-256: `bfef1abd56d22c94d2cc01a679fe9b9a44dcf10fc67b89e99f289ff9a7e253b8`.

The fixed-geometry GPAW calculation converged after 52 total SCF iterations (40 steps in the first window and 12 after checkpoint resume). GPAW 26.7.0, ASE 3.29.0, PBE, PW 500 eV, standard Monkhorst-Pack 7x2x1, Fermi-Dirac 0.10 eV, Davidson niter=2, MSR1 beta=0.05/nmaxold=10/weight=20; 4 MPI ranks, OMP=1. No atomic or cell relaxation.

- Total energy: **-939.261473902527 eV**.
- Force-consistent free energy: **-940.098020842865 eV**.
- Maximum force magnitude: **1.727108 eV/Angstrom** (frozen, unrelaxed slab).
- Stress (eV/Angstrom^3): `[[0.00498956283573553,-4.18327418219964e-05,-0.00010100236452401131],[-4.18327418219964e-05,-0.0021317735613340027,2.1140731755320147e-05],[-0.00010100236452401131,2.1140731755320147e-05,-0.00023509936022669473]]`.

The 40-step full-wavefunction checkpoint was hashed and read back successfully before resumption; no partial energies or forces were used.

## TSCBOP frozen-coordinate comparison

A single LAMMPS `run 0` was performed for the displaced geometry. LAMMPS 22 Jul 2025 Update 4; `units metal`, `pair_style tersoff`, element mapping Ti/Si/C. The source TSCBOP SHA-256 remains `729273a6488d75756ed7927432c77c7347a1b115cecbaefae4981bafcec80fd7`. The simulation advanced zero time steps and did not minimize coordinates. Source IDs and elements mapped exactly; the maximum wrapped-coordinate difference was 1.422e-14 Angstrom.

Raw TSCBOP-vs-DFT force-vector RMSE at this displaced geometry is **1.947398 eV/Angstrom** (component RMSE 1.124331). This large absolute mismatch remains; a paired response can look much better because the baseline error cancels.

For the paired displacement response, defined as `(F_TSCBOP(displaced)-F_TSCBOP(reference)) - (F_DFT(displaced)-F_DFT(reference))`, the all-atom vector RMSE is **0.023213 eV/Angstrom**, P95 0.041815, maximum 0.129738 at source atom 119. Within 3.3 Angstrom of either displaced Ti (18 atoms), vector RMSE is **0.059368 eV/Angstrom**. Element vector RMSE: Ti 0.020940, Si 0.007061, C 0.030559 eV/Angstrom.

For context, the previously completed paired normal (z) Ti displacement in this same parent family had all-atom response vector RMSE 0.011171 eV/Angstrom and local (<3.3 Angstrom) RMSE 0.026783; the present tangential response is larger. These are two directional probes of one source family, not independent validation statistics.

## Interpretation and limits

The unmodified TSCBOP gives a much larger absolute-force mismatch than its differential response error. It captures much of this one small tangential response mode after subtracting the same-parent reference, but the local maximum and larger error than the normal probe show the response is not exact. No broad transferability claim follows. No TRAIN/DEV/TEST assignment changed. No potential fitting, additional DFT, MD, or TTM was performed.

Detailed machine-readable metrics are in `analysis/TSCBOP_RUN0_TANGENTIAL_20261023/SUMMARY.json`; per-atom data are in `TANGENTIAL_FORCE_COMPARISON.csv` and `TANGENTIAL_RESPONSE_DIFFERENCE.csv`.
