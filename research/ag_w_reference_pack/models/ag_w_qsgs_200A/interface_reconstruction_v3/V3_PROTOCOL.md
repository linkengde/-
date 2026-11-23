# Ag–W local interface reconstruction v3: preregistered protocol

Date: 2026-11-23
Status: criteria fixed before the v3 interface relaxation runs.

## Scope

V3 starts with small planar Ag/W slabs, not the 482,471-atom QSGS model. It tests whether the Dai EFS potential can represent the two literature interfaces before any three-dimensional rebuild. The raw QSGS geometry and v1/v2 audits remain untouched.

The primary references are He et al., *J. Phys. Chem. C* (2018), Table 4: Ag(100)/W(100), `Wsep=2.80` and `gamma_int=0.72 J/m²`; Ag(111)/W(110), `Wsep=2.45` and `gamma_int=0.86 J/m²`. The paper reports seven atomic layers as sufficient for its surface-energy convergence tests, 15 Å vacuum, and about 3.44% in-plane mismatch. It does not publish the complete interface supercell vectors, exact strain partition, atomic coordinates, or a reusable DFT input deck. Therefore, v3 values compared to Table 4 remain a screening comparison, not an exact reproduction.

## Criteria fixed before relaxation

### A. Coincidence-cell selection

1. Keep the paper's stated orientation relationships and independently search integer repeats along the two orthogonal in-plane vector pairs.
2. Use the validated Dai-EFS equilibrium lattice constants: Ag 4.0900219316 Å and W 3.1599967891 Å.
3. Assign the common vector length as the geometric mean of the two natural repeat lengths. This is a symmetric log-strain convention, not a claim of elastic-energy optimality.
4. A candidate cell must have no more than **1.00% absolute in-plane strain in either phase or direction** and area no larger than **1,800 Å²**. Choose the smallest-area cell satisfying both; never relax these limits after seeing energies. The limits are project screening rules, not universal physical constants.
5. If no cell meets the predeclared limits, stop that orientation and report it; do not force a higher-strain fit.

### B. Local atomic relaxation and interface thermodynamics

1. Use seven atomic planes per phase and 15 Å vacuum, consistent with the source paper's slab-thickness check. Preserve every atom after constructing the slab. Enumerate every distinct sequence offset through one oriented bulk repeat: 2 Ag × 2 W termination combinations for (100)/(100), and 3 Ag × 2 W combinations for (111)/(110).
2. Hold the central three planes of each phase fixed; relax all other atoms with the unconstrained Dai EFS potential. No temporary springs or soft-core forces are allowed in final relaxation.
3. For Ag(100)/W(100), enumerate top `(0,0)`, half of the first matched Ag surface vector `(1/2,0)`, and BCC-type `(1/2,1/2)` translations. For Ag(111)/W(110), sample the four coarse half-vector translations `(0,0)`, `(1/2,0)`, `(0,1/2)`, `(1/2,1/2)` of the matched Ag oriented surface cell; these are translation samples, not claims that the paper's exact atomic registries were reconstructed. Scan initial interface gaps 2.0–3.6 Å at 0.2 Å increments plus 6.0 Å separated reference. Keep all minima and do not select only the value closest to DFT.
4. A relaxed point passes the mechanical gate only when LAMMPS reports force-tolerance termination and maximum mobile-atom force is at most **0.03 eV/Å**, the force tolerance stated in the DFT paper. No atom may be lost. Report minimum same-species and cross-species distances; flag a same-species distance below 85% of the corresponding Dai-EFS equilibrium nearest-neighbor spacing as local collapse.
5. Compute surface references with the same strained in-plane cell, seven planes, vacuum, and core constraints. Compute bulk references with the same in-plane strain and the out-of-plane repeat relaxed. Check `Wsep = gamma_Ag_surface + gamma_W_surface - gamma_int` closes within **0.01 J/m²**.
6. For a provisional potential-screening pass, both `Wsep` and `gamma_int` must be within **±0.30 J/m²** of the published values for each orientation, with positive `gamma_int`. This is a project screening tolerance, not an uncertainty quoted by He et al.; it will not be loosened after results are known. Because the exact DFT cells are unavailable, satisfying it would still not establish an apples-to-apples reproduction.

### C. Gate before any 3D expansion

For the **quantitative-validation track**, do not expand to QSGS until both orientations pass the local mechanical and thermodynamic gates. Any later QSGS rebuild must preserve the original material box, add/delete no atoms, keep Ag/W mass ratio within 0.1% of 1:1 and each phase's number density within 1% of its corresponding target, and retain both 3D-connected phases (largest component at least 98%, x/y winding rank 2, and spanning the 200 Å material thickness). Check static forces and local collapse after full relaxation. Do not start 300 K MD or TTM until these structural checks pass.

After the v3 screen, the project scope was clarified to permit **qualitative graduate-level exploratory use** despite disagreement with DFT. This does not change the numerical gates or convert the quantitative screen to a pass. It permits development of a separately labeled exploratory 3D candidate after it passes its own static geometry, density, atom-count, and connectivity checks. Any such results must carry the known potential limitations; quantitative interface-property claims remain outside scope.

## Frozen construction choices

- `Ag(100)/W(100)`: `[010]Ag || [011]W`, as reported by He et al.
- `Ag(111)/W(110)`: `[-2,1,1]Ag || [1,-1,0]W`, as reported by He et al.
- Match search is recorded by `search_low_strain_cells.py`; all selected repeats, common lengths, area, and phase strains are saved to `cell_match_report.json`.
- The first run uses one Ag/W boundary and two outer free surfaces in a z-periodic cell with 15 Å vacuum. It is compared to the paper's one-interface slab convention. The surface and bulk references are generated from the same EFS-strained structures.

## Known limits

- The 1%/1,800 Å², 0.03 eV/Å, 85%, 0.01 J/m², and ±0.30 J/m² gates above were written before the v3 relaxation outputs. They must be cited as project criteria, not literature standards (except the DFT paper's 0.03 eV/Å force tolerance and its seven-layer/15 Å setup).
- The match search balances logarithmic strain equally between phases because verified orientation-specific elastic constants are not yet available. It does not claim an elastic-energy optimum.
- The EFS potential already has an unresolved B2 Table 7 energy discrepancy and prior Ag(111)/W(110) interface disagreement. A v3 pass would be a limited interface screen only.
