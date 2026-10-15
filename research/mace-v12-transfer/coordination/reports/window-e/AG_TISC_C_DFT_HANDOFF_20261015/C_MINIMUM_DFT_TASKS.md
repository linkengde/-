# C window: minimum supplemental DFT task list

## Decision gate before spending compute

The audit says `CROSS_POTENTIAL_DATA_READY=NO`, `MORSE_FIT_FEASIBLE=NO`, `ADDITIONAL_DFT_REQUIRED=YES`. First resolve the recorded Ag–C model-form issue: the archived central-radial/Morse lower-bound test leaves 0.906–0.909 eV/Å vector residual on two stored Ag–C diagnostic inputs against a 0.05 eV/Å gate. C should explicitly retain/change the pair form or the gate before authorizing a fit-oriented acquisition. This is a decision gate, not a new audit or a claim that a new potential is fitted.

## What already exists

| Existing record | Verified content | Reuse boundary |
|---|---|---|
| 135-atom relaxed full interface | 1 converged PBE/PAW record on Ag21Ti57Si19C38; native energy −939.821159610571 eV, force-consistent free energy −941.236116683098 eV, 135 forces, stress; final structure SHA `eda2586c63c6c3c793ffce99417e7c04569941b0a2a916fda1b5402b7ba52d9f`; minimum Ag–Si 2.52319455 Å, Ag–Ti 3.8069 Å, Ag–C 5.0949 Å. | One structural family/registry anchor. Counts as one TRAIN grid point only if C's frozen target grid actually includes it. Not an independent validation family. |
| Matched Ag-only and TSC-only controls | Both converged on the same final cell/PBC and retained atom coordinates/IDs as the 135-atom full result. Full per-atom forces, logs and energies included. | Reusable for rigidly translated/separated versions only after exact component cell/PBC/coordinate/ID invariance is confirmed. Not a general clean-surface reference. Ag21 is neutral odd-electron, non-spin-polarized with FD fractional occupations; preserve that caveat. |
| Original 2.000 Å full interface single point | A separate converged full-system fixed-geometry result on the original (unrelaxed) 135-atom input; max force 7.2162 eV/Å. | No matched Ag-only/TSC-only outputs at that geometry. Do not use as an interaction-energy/cross-force label unless C adds matching component controls and accepts the strongly repulsive geometry. It may be considered as a compressed-side point only after review. |
| V19 reference bundle | Complete verified archive of 16 converged numerical records: six local Ag–Ti/Si/C acquisitions and ten k-point/smearing controls on one 52-atom geometry. | Six acquisitions are from reused training-family local references and lack matched phase-only DFT; the ten controls are one geometry, not ten structures. Do not count these as independent validation or repeat them. |
| Existing failed Γ-centered 2×2×1 attempt | Failure manifest/logs only. | No usable energy/force/stress label. Do not count as completed data. |

DFT interaction-energy values from the matched 135-atom final geometry are already in the audit inventory. Native-energy difference (−22.229525474612 eV) and force-consistent free-energy difference (−21.378661277249 eV) are distinct and must remain separate. Never compare unrelated absolute energies across compositions.

## Prioritized minimum acquisition

| Priority | C action | Minimum job set | Output / stop gate |
|---|---|---:|---|
| P0 | Decide the accepted Ag–C functional form and whether the 0.05 eV/Å vector gate remains binding, in light of the archived lower-bound evidence. Freeze target termination(s), registries and normal gap values. | No DFT | If unresolved, keep acquisition/fit on HOLD. |
| P1 | Build and freeze exact, hashed structures for a TRAIN coverage floor: two non-equivalent lateral registries × three normal separations (repulsive/compressed, near-minimum, expanded/attractive). Keep the TSC slab fixed and describe the rigid/allowed transformations. Ensure the selected geometry family activates Ag–Ti, Ag–Si and Ag–C at useful distances; if not, add a justified pair-relevant termination/gap series. | 6 full-interface geometries | New structures are NOT included in this handoff; C must generate/review/freeze them before authorization. Exact gap values are intentionally not guessed here. |
| P2 | Run fixed-geometry full-interface DFT for the TRAIN grid. Reuse the existing final 135-atom full/component point only if it occupies an actual frozen grid slot. | 5 new full SCFs if the existing point occupies one slot; otherwise 6 new full SCFs | Keep native energy, free energy, forces, stress if requested, convergence, atom IDs, cell/PBC and hashes. No optimization. |
| P3 | Freeze an independently constructed parent surface/registry family for held-out validation, at the same three gap classes. It must exercise the same pair environments and cannot share a parent structure with TRAIN. C controls sealed roles. | 3 full-interface SCFs plus up to 2 phase-only SCFs (Ag-only and TSC-only) | Existing component references may be reused only if exact phase geometry/cell/PBC/IDs are invariant under the permitted rigid transforms. Otherwise run matched neutral phase controls. No role changes without C's registry. |
| P4 | Verify coverage and numerical rank using only static model-side calculations after data arrive; inspect pair-specific distance/force coverage and train/held-out family separation before any fitting request. | No DFT unless coverage is deficient | If Ag–C is not meaningfully sampled by the P1 set, add a separate pair-relevant termination with a three-gap curve and matched components; the baseline job count is then insufficient. |

### Baseline count and conditions

The audit's lower-bound estimate is **at least 10 new SCF single points** under the narrow assumptions that (a) the existing relaxed 135-atom point is on the six-point TRAIN grid, (b) its Ag-only/TSC-only references exactly transfer over the TRAIN rigid transformations, and (c) the independent validation family needs two new component references. That is 5 new TRAIN full + 3 held-out full + 2 held-out component jobs. If the 135-atom point does not occupy the frozen grid, add one full-interface job (at least 11 new SCFs). If TRAIN phase geometries change beyond verified rigid transforms, add the missing component jobs. If pair coverage requires a new Ag–C termination/gap series, add its full and component calculations; exact count depends on frozen structures.

This is a coverage floor, not a guarantee that all Morse parameters are identifiable or that the force threshold can be met. Existing V19 acquisitions and controls are already available and must not be duplicated.

## Initial common DFT method proposed by the archived request

- GPAW 26.7.0, ASE 3.29.0; PBE/PAW.
- Plane-wave cutoff 500 eV; Γ-centered 1×1×1; Fermi–Dirac width 0.10 eV.
- SCF thresholds: energy, density, eigenstates each 1e-5 as archived.
- Fixed atomic positions and cell; no relaxation, MD, or thermal trajectory.
- Record native total energy and force-consistent free energy separately, all three force components per atom, stress if requested, SCF status/steps, occupations/electron count, atom IDs/order, composition, cell/PBC, exact input/result/configuration SHA-256.
- Phase-only references must be neutral and match cell/PBC/retained coordinates and PAW/numerical settings. Flag odd-electron/spin/occupation concerns instead of concealing them.
- C may issue a corrected locked method before any execution authorization. Exact frozen inputs and an authorization bound to their SHA values are still missing.

## Ready-to-execute status

`EXISTING_REFERENCES_PRESENT=YES`  
`NEW_TRAIN_AND_VALIDATION_INPUTS_FROZEN=NO`  
`DFT_AUTHORIZATION_PRESENT=NO`  
`DFT_EXECUTION_STARTED=NO`
