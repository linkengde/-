# Stage C-D1 final DFT input audit

## Input identity and geometry

Source package: C commit `1fec365ecedb8555ae2c857509d1018953639b86`, directory `research/mace-v12-transfer/coordination/reports/window-c/ti3sic2_ag_interface_fit_dft_inputs_20261014/`. The source package's complete `SHA256SUMS.txt` was rechecked successfully before this audit. `STRUCTURE_RESULTS.csv` records per-structure geometry/config hashes and measured metadata.

| Group | Count | Identity | Ordered labels / PBC | Minimum periodic separation |
|---|---:|---|---|---|
| A — Ti₃SiC₂ backbone | 6 | 48 atoms, Ti24Si8C16 | Ordered species and unique `local_new_id` checked; TTT | 2.074173–2.102018 Å (matches E7 geometry audit) |
| B — Ag/Ti₃SiC₂ interface | 9 | 217 atoms, Ag84Ti57Si38C38 | Ordered species and unique `atom_id` checked; TTT | 1.8 Å minimum across proposed gap scan; no pair below 1.5 Å |

For the interface scan, the parent structure SHA is `945c5e214c1ceaeba4125ff785f4092e14858c75a985061a7841da3e5ce6a9a2`. The 9 inputs differ only by rigidly translating the Ag layer to plane gaps 1.8–3.4 Å; no ions, cell, lateral registry, atom IDs, or PBC were changed within the input package. Cell matrices and `pbc` are recorded per row. No input contains DFT labels.

## Electronic-structure parameters

| Setting | A — bulk backbone | B — interface slab |
|---|---|---|
| Software | GPAW 26.7.0; ASE 3.29.0; gpaw-data 1.2.1 | Same |
| Method | PBE / PAW / plane waves | Same |
| Cutoff | 500 eV | 500 eV |
| Occupations | Fermi–Dirac, 0.10 eV | Fermi–Dirac, 0.10 eV |
| SCF thresholds | Energy, density, eigenstates each 1e-5; maxiter 160 | Same |
| Mixer | Pulay β=0.05, nmaxold=8, weight=100 | Same |
| k mesh | 10×10×4; GPAW default Monkhorst–Pack integer mesh, no added offset | 2×2×1 Gamma-centered, carried forward from prior geometry proposal; **not convergence-validated** |
| Parallel | 4 MPI ranks, 1 OpenMP thread/rank | Same |

PAW file hashes were measured in the active setup directory and match the configs: Ag `bff00ac79a7fa6e8c3a71c8b34d3a9c488e6bf1c00b3e9c4c91c231551145e97`; Ti `bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f`; Si `2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0`; C `fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983`.

## Audit decision

- Structure identities, order, IDs, cells, PBC, input hashes, cutoff, smearing, SCF thresholds and PAW setup identities: **PASS**.
- Group A method matches the existing E6/B fixed-cell bulk protocol; calculation resource release: **HOLD**.
- Group B method fields are explicit, but 2×2×1 convergence is unproven; calculation release: **HOLD** pending a same-geometry denser in-plane comparison and a stated numerical tolerance, plus adequate resource headroom.
- Geometry is fixed; no relaxation, atom motion, or cell optimization is included.

The report does not claim that merely having an input file makes its DFT label numerically adequate.
