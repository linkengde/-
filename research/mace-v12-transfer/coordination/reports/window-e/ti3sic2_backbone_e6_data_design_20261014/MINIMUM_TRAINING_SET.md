# Minimum training-set and family split proposal

## Current archive: reusable, but not a role assignment

- Reuse the converged 48-atom bulk baseline and the existing Ti/C `±0.02 Å` z-displacement labels after the E2 provenance/units audit. These five geometries share one parent family; all k-grid variants are numerical checks of the same geometry, not extra independent configurations.
- The 114-atom TSC-only slab is derived from a 135-atom Ag interface, frozen, and has a reported maximum backbone force near `0.976 eV/Å`. Keep it as a diagnostic comparison. It is not a clean-surface reference and does not supply an independent test family.
- No current records are assigned to TRAIN, DEV or TEST. Preserve the empty role freeze until C approves a split manifest.

## Proposed family-level roles

| Role | Proposed source families | Keep together | Purpose |
|---|---|---|---|
| TRAIN | The audited 48-atom bulk parent family, including its baseline, Ti/C/Si displacement descendants and volume/shear variants; plus the chosen clean Ti- and Si-terminated surface training families and all their local-response descendants. | All variants from each root structure stay in TRAIN. Do not split ± displacement pairs, volume/shear pairs, k-grid controls, or surface relaxations from their root. | Fit bulk force response, elastic stiffness and the two target surface terminations. |
| DEV | A separately rooted, independently constructed surface/defect family with a reference and a small deterministic response configuration. | Keep all reference and perturbation descendants in DEV. It must not be cut from, mapped from, or relaxed from a TRAIN root. | Select model size, regularization and weights. |
| TEST | A second distinct root family, preferably an Ag-adjacent solid-skeleton/interface environment with an independently archived DFT reference and a response point; otherwise an independently generated clean-surface/defect family matching the target domain. | Keep all descendants in TEST; do not inspect it during model or weight selection. | One-time final check of transfer to a new family and target environment. |

A different supercell, k-grid, atom ordering, or small displacement of the same parent is not an independent family. If a DEV/TEST family cannot be obtained without shared lineage, report independent validation unavailable and keep the fit exploratory/blocked. Proposed roles above are for a future C-reviewed manifest only; no current role was changed.

## Minimum coverage matrix

For any fit claim covering the stated solid-backbone goals, cover bulk forces, strain response, Si local response, and surfaces. An algebraically identifiable force subset without these targets is not a material-quality fit.

| Coverage | Minimum useful observations | Parameter/data purpose | Limitation if omitted |
|---|---|---|---|
| Bulk equilibrium and local Ti/C response | Existing audited baseline + Ti/C displacement family; no duplicate k-grid entries. | Anchor local Ti–C/Ti–Si force directions. | One parent family cannot establish transfer or robust bounds. |
| Hydrostatic and shear response | Both signs of volume and shear proposals (4 static DFT points) with energy and stress. | Constrain curvature/elastic response and distinguish amplitude/decay changes. | Forces from the unstrained parent do not constrain elastic tensor. |
| Si local response | Symmetry-audited Si site, `±0.02 Å` (2 static DFT points). | Add a missing element/site response and improve Ti–Si directions. | Ti/C displacement labels cannot replace Si force response. |
| Surface response | Ti- and Si-terminated clean slabs, each reference plus a deterministic local response (4 static DFT points). | Constrain under-coordinated backbone forces for the target surface/interface role. | The existing interface-derived slab is not an equilibrium clean surface. |
| DEV and TEST | Two entirely new source families; recommended reference + local-response pair in each (4 static DFT points). | Tune on DEV and reserve TEST for one final transfer check. | One point per family (2 total) is a lower bound but cannot test within-family response; no family independence means no generalization claim. |
| Si–C block, conditional | Only if the intended domain has Si–C pairs within taper: at least two distinct active distance/coordination states in a training family and a physically distinct held-out active family; recompute the Jacobian before selecting Si–C A/B. | Activate otherwise-zero Si–C columns and test distance dependence. | Artificially compressed contacts or one active snapshot cannot justify fitting both Si–C A and B. |

The recommended known-minimum supplement is 14 new static DFT points: 4 volume/shear + 2 Si displacement + 4 clean-surface reference/response + 2 DEV + 2 TEST. This is consistent with E3's recommended plan. A 12-point lower bound reduces DEV and TEST to one point each and weakens validation. If the surface structures do not naturally activate the Si–C taper and Si–C is needed for the intended domain, add a separately approved active-contact DFT family; its count is conditional and not included in 14. A trusted PBE bulk equilibrium lattice/internal-coordinate reference is also conditional: the audited CIF-derived baseline is not established as a PBE equilibrium structure.
