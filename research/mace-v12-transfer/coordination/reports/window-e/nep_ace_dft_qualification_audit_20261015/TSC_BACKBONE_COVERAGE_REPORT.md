# Ti₃SiC₂ backbone coverage

## Confirmed DFT-labelled coverage

Five unique 48-atom bulk geometries are confirmed with converged PBE energy and forces: one baseline parent plus Ti and C displacements of ±0.02 Å. The raw archive contains seven k-mesh calculation variants of the baseline geometry; those variants are one structure, not seven training samples. The 48-atom cell is C16Si8Ti24. The Ti/C displacements provide a small local-force response set. Explicit stress/virial is absent for these five bulk controls.

The E2 archive also contains a 114-atom Ti₃SiC₂-only subset and a 135-atom interface-derived slab/interface parent. Their forces are useful for that frozen local geometry and decomposition, but do not form independent clean-surface or bulk families. The legacy surface-distance family has 13 geometries, of which one is classified as a complete confirmed DFT label and 12 lack enough source/convergence evidence for strict use. They are not a validated independent surface set.

The V19 repository contains proposed Ti₃SiC₂ volume/shear and additional structure points, but the inventory marks the proposals without converged DFT labels as `NON_DFT_DATA`. No confirmed Si displacement pair was found.

## Missing physics

No complete Si displacement response, Ti₃SiC₂ volume or shear stress response, independent clean surface/termination family, or decorrelated high-temperature solid-vibration family is available with verified DFT labels. No curved/thin-neck surface environments representative of a bicontinuous framework are covered.

## Decision

`TSC_BACKBONE_DATA_READY = INSUFFICIENT_DATA`. Current labels can support a format/pipeline smoke test and a narrow bulk-force baseline. They cannot establish elastic response, Si response, free-surface transfer, thermal robustness or independent generalization. Additional DFT is required before trusting a specialized potential to preserve the solid backbone during Ag melt/freeze simulations.
