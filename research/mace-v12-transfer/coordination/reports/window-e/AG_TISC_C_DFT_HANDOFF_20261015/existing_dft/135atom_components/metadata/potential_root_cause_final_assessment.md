# 135-atom interface force-error attribution: final assessment

## Evidence status

| Check | Status | Basis |
|---|---|---|
| `FORCE_MAPPING_VALID` | PASS | Full DFT, component DFT and LAMMPS rows map to source IDs, element order and frozen coordinates. The archived full-force columns match the GPAW force file exactly; LAMMPS force decomposition closes within `8.8e-13 eV/Å`. Units are eV/Å throughout. |
| `TSCBOP_IMPLEMENTATION_VALID` | PASS (implementation only) | Prior audit verified the TSCBOP hash and Ti/Si/C/NULL Ag mapping. The present large TSC-only force residual is evidence about accuracy, not an implementation error. |
| `CROSS_FORCE_IMPLEMENTATION_VALID` | PASS (implementation only) | Prior `run 0` audit verified EAM/TSCBOP/Morse routing and the Morse difference force. The remaining effective-cross residual is a model mismatch, not a force-mapping defect. |
| `TSC_INTRAPHASE_ERROR_ISOLATED` | PASS for this frozen geometry | Direct comparison of 114-atom TSC-only DFT and TSCBOP-only LAMMPS gives 2.5624 eV/Å vector RMSE; Ti 2.6242 and C 2.4433 eV/Å. |
| `CROSS_ERROR_ISOLATED` | PASS as an effective-force comparison | Full-minus-phase-only DFT forces can now be compared with LAMMPS Morse forces. All-atom effective-cross RMSE is 1.3513 eV/Å; interface subset RMSE is 1.5178 eV/Å. This is not a unique DFT pair-potential decomposition. |
| `READY_FOR_POTENTIAL_REFIT` | **NO** | Only one interface geometry and one registry have component-resolved references. This is inadequate coverage for a transferable fit; no fitting was done. |

## Root cause on this structure

The most important finding is that Ti/C force errors are mainly carried by the intraphase skeleton term: their TSCBOP-only errors are 2.624/2.443 eV/Å, while the Morse-vs-effective-DFT-cross errors are 0.332/0.095 eV/Å. That rules out changing only Ag cross terms as a sufficient repair for this structure.

The Si termination and Ag atoms show the opposite emphasis. Si has a 2.742 eV/Å effective-cross discrepancy, and Ag has a 2.150 eV/Å effective-cross discrepancy, despite the Ag-only EAM discrepancy being 0.0565 eV/Å. Cross-term work is therefore justified for Ag/Si, but it must follow separate backbone review. The maximum effective-cross error is 2.954 eV/Å at Si atom 76, 2.523 Å from Ag. Full hybrid-model force RMSE remains 2.2547 eV/Å.

These numbers describe the one optimized, frozen 135-atom geometry under Γ-only PBE/PAW and the current non-spin-polarized electronic setup. They do not establish transferability or predict the errors on other interfaces.

## Additional DFT configurations needed before a fit

For a useful cross-term validation set, add full-system fixed-geometry DFT forces at:

1. At least two additional, non-equivalent lateral Ag registries over the same Si-terminated surface.
2. Several normal separations around the current relaxed gap, including one compressed point and points on both sides of equilibrium. A practical first grid is three separations per registry, spanning roughly ±0.15–0.30 Å from the local minimum, with the exact distances frozen before calculation.
3. Two or more deterministic small displacements of interface atoms (about 0.02–0.05 Å) to test local force response rather than one stationary geometry.
4. One separated-slab configuration with the Ag/substrate distance beyond the Morse cutoff and verified no interaction with periodic images, to test the zero-cross limit.

For attribution at each new gap/registry, also calculate neutral Ag-only and TSC-only references at the matching frozen cell/coordinates. For a direct full-potential fit, full-system force labels across those configurations are the essential data; phase-only references are required when continuing to diagnose cross versus intraphase contributions. The TSC-only discrepancy also warrants a small strained/surface Ti₃SiC₂ control set before considering a backbone refit. Any new calculations or fitting remain separate work and were not started here.

## Answers

- TSCBOP-only force RMSE: **2.5624 eV/Å** (Ti 2.6242, Si 2.6072, C 2.4433).
- EAM-only Ag force RMSE: **0.0565 eV/Å**.
- Morse effective-cross force RMSE: **1.3513 eV/Å** over all atoms; **1.5178 eV/Å** on the ≤8 Å interface subset.
- Ti/C error is mainly in the skeleton term on this geometry. Si has substantial skeleton and cross discrepancies; Ag is dominated by cross discrepancy.
- The current result supplies a force-attribution reference, but it is **not enough** to fit a transferable cross potential.
- Needed next: alternate registries, a gap grid, deterministic near-equilibrium perturbations, a separated reference, and additional strained/surface TSC-only controls. No refit or new DFT was started after these two component single points.

