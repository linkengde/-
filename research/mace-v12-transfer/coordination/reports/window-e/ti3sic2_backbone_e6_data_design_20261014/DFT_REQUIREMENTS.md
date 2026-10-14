# DFT requirements and existing-run0 boundary

## Existing labels to reuse

Reuse the actual converged PBE/PAW force labels already audited in E2/E4 and consumed by E5. Do not rerun them for this design. The five 48-atom bulk cases are one source family. The 114-atom TSC-only slab can remain a diagnostic family but cannot stand in for independent clean surfaces. The DFT setup and units must remain method-matched; record total energy and force-consistent free energy separately.

## New static DFT labels required for stated coverage

| Priority | Proposed data | Minimum count | Why existing data cannot substitute | Status |
|---|---|---:|---|---|
| 1 | 48-atom parent volume `−2%/+2%` and shear `−0.5%/+0.5%` frozen-geometry points; collect forces, energy, and stress. | 4 | Existing 48-atom force archives have no stress and do not sample strain curvature. The proposal inputs are not DFT labels. | Inputs proposed in E3; DFT not run in this stage. |
| 2 | One symmetry-audited Si site, `±0.02 Å`, same parent and frozen settings. | 2 | No Si-displacement DFT labels exist; Ti/C movement cannot supply Si's local force response. | Missing. |
| 3 | Ti- and Si-terminated clean surface reference and local-response points, with frozen cell/PBC/vacuum and a declared surface convention. | 4 | Existing 114-atom structure is derived by deleting Ag from a relaxed interface; it is frozen, highly forced, and has no independent clean-surface equilibrium reference. | Missing. |
| 4 | Independent DEV and TEST families, each preferably a reference plus deterministic response point. | 4 recommended; 2 weak lower bound | Current roots are one bulk parent plus one derived slab; k-mesh changes do not create independent geometry. | Missing. |
| Conditional | A physically relevant Si–C-active structure family, only if target-domain screening shows the pair is inside the Tersoff taper. | No fixed count; likely multiple DFT-labelled states plus held-out family | Current minimum distance exceeds cutoff, so E5 provides no Si–C target response. A potential-only `run 0` cannot create target forces. | Activation geometry not identified; Si–C should stay frozen meanwhile. |
| Conditional | Establish a trusted PBE equilibrium cell/internal-coordinate reference if the audited CIF-derived structure is not replaceable by an accepted reference. | 1 variable-cell/internal-coordinate workflow or trusted audited equivalent | Current fixed CIF-derived bulk point has residual force and is not verified as PBE equilibrium. | Source/reference gap remains. |

The four volume/shear proposals and two Si points use one bulk parent family and must all remain in one role. Every surface or perturbation family must likewise remain intact across role assignment.

## What can be done with LAMMPS `run 0` alone

- For existing DFT-labelled geometries, new tied-block potential variants can add force-response columns without a DFT rerun. E5 already performed the requested eight-block ±1% analysis on the current six structures.
- A `run 0` on a candidate surface/defect geometry can check whether Si–C distances activate the unchanged taper. This is a screening step, not a DFT label and not a scientific validation.
- `run 0` cannot provide DFT force targets, stress, energy differences, PBE equilibrium response, or independent DEV/TEST evidence. The new geometry families in the table require DFT if they are to appear in a fit or validation objective.

No DFT, LAMMPS run, fitting, MD or TTM was launched for E6.
