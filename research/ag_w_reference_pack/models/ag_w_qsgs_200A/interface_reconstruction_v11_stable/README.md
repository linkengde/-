# Ag–W V11 Stable Reference Model

**Purpose:** produce a runnable exploratory Ag–W bicontinuous reference model for TTM-MD and Gaussian laser heating. V11 is allowed to make small, explicitly measured changes to the simulation cell and, only as a last resort, a small mass-balanced atom-count adjustment.

**Status at branch creation:** execution plan committed; no V11 structure has yet been generated and no LAMMPS/TTM run has yet been executed. Do not interpret this directory as a passing model.

## Frozen evidence

- V10 parent commit: `90ddbf8bbbf4cdf29677771dd5996f52e8c10000`.
- V9 full-structure state: 482,601 atoms; Ag 304,046; W 178,555; Ag/W mass ratio 0.999127; 200 × 200 × 300 Å; 46,117 original fixed-bottom IDs.
- V10 5 Å mapping, best phase (0,0) Å: 1,181.24 Å³ remained unassigned within 20 Å.
- V10 10 Å mapping: 6,533.53 Å³ remained unassigned within 20 Å.
- V10 40 Å diagnostic could assign nearly all volume, but its 34.87 Å volume-weighted mean transfer distance is not acceptable as the production map.
- V10's highest-overflow 5 Å cells are concentrated at the first/last z layers (centres z=42.5 and 197.5 Å). This motivates testing z guard cells before modifying atoms.
- V10 found 7,619 atoms from the static Newton group outside the old 40–200 Å mapping slab. The mapping domain should cover the actual material envelope rather than silently dropping these atoms.

## V11 decision order

1. **Fix the mapping domain first, not the material.** Test a 5 Å grid with one guard layer below and above the old slab: z=35–205 Å, shape 40 × 40 × 34, spacing 5 Å. Include every original static-Newton atom currently inside this padded region; report and reject any excluded atom. Keep the same 20 Å direct-transfer limit and the same Ag/W atomic-volume constants.
2. Scan all four 5 Å x/y phases (0,0), (0,2.5), (2.5,0), (2.5,2.5) Å with the exact production max-flow solver. Run an actual LAMMPS mapping-only check for the best candidate. Do not count guard-cell capacity as material volume; empty guard cells remain zero-fraction cells.
3. Recompute Ag-only Gaussian source power on the complete mapped Ag set. For this exploratory nominal-plane source, keep the original z=200 Å reference and explicitly enable the documented above-plane depth clamp for protruding surface atoms; record this approximation. Do not silently move the surface plane to 205 Å.
4. If the padded grid still fails, test small in-plane cell dilations of 0.25%, 0.50%, and 0.75% with coordinates scaled about the x/y box centre. Preserve z and the free-surface vacuum. Recompute grid spacing from the new periodic box and preserve 40 × 40 x/y cells; use anisotropic spacing if needed and supported by the solver. Re-run the full four-phase mapping scan after each candidate.
5. Only if all preceding candidates fail, permit a targeted mass-balanced adjustment of no more than 0.1% of atoms. Select candidates from the actual residual donor cells, not randomly. Preserve both Ag and W percolation and keep whole-model Ag/W mass ratio within 0.1% of unity. Record every removed/added ID and its reason. Do not change the EAM atomic-volume constants or relax the per-cell capacity bound.
6. For the selected geometry, perform a spring-free energy minimization with the original bottom fixed IDs. Then reinitialize all mobile atoms from one common 300 K Maxwell distribution (one thermostat for all mobile atoms, not separate Ag/W thermostats), remove initial COM momentum, run NVT for 2 ps and extend to at most 5 ps if the potential-energy slope or phase temperatures have not plateaued. Follow with 0.1 ps NVE.
7. Only after mapping and cold-baseline gates pass, run source/no-source paired TTM smoke tests and a low-energy exploratory Gaussian pulse. Recompute the source integral on the final geometry; never reuse V7/V8/V9 source power.

## Required acceptance

- Whole-model Ag/W mass ratio: 0.999–1.001.
- No severe same-phase or cross-phase overlaps; report nearest-neighbour distributions and force hot spots.
- Both phases retain 3D/xy percolation appropriate to the original QSGS target.
- Every mapped Newton atom is either mapped or explicitly identified as outside the declared physical mapping domain; no silent clipping.
- Per-cell Ag+W mapped fraction ≤ 1 within the solver tolerance.
- Direct transfer ≤ 20 Å for the accepted production mapping; report actual volume-weighted and maximum transfer distances.
- Source power agrees between atom quadrature and grid integration; W direct source remains zero.
- NVT tail is plateau-like; Ag/W temperatures are reported with both standard and COM-removed conventions. Do not force phase temperatures equal with separate thermostats.
- NVE has no runaway energy drift or atom loss.
- TTM source/no-source pair conserves energy to the established numerical tolerance. Label the provisional Ag/W interface conductance and pulse envelope as exploratory.

## Reproducibility

Keep V5–V10 untouched. Store each geometry variant, data/restart SHA-256, LAMMPS inputs/logs, phase-scan JSON, mapping-only JSON, cold metrics, TTM accounting, and a complete SHA256SUMS file in this directory. Each candidate must have an explicit PASS/FAIL gate. No step may be marked PASS unless it was actually executed.
