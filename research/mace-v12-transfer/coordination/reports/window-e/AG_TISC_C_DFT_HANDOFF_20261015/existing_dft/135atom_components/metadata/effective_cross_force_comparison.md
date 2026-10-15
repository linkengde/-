# Effective Ag–Ti₃SiC₂ cross-force comparison

## Definition and scope

All forces are compared on the same frozen optimized 135-atom structure, with source atom IDs 1–135, the same periodic cell and PBC, and units eV/Å. The DFT full forces are from the converged GPAW geometry result. The component-only forces are from the two converged fixed-geometry PBE/PAW single points in this package. LAMMPS component forces are the already audited A/B/C `run 0` results using the unchanged EAM, TSCBOP and calibration-v0 Morse files.

For each source atom, the diagnostic effective DFT cross force is defined by

```text
Ag:       F_cross_DFT = F_full_DFT − F_Ag-only_DFT
Ti/Si/C:  F_cross_DFT = F_full_DFT − F_Ti3SiC2-only_DFT
```

The corresponding LAMMPS cross force is the existing `F_cross_Morse` column, independently checked against `F_full − F_Ag-only − F_TSCBOP` (maximum closure error `8.8e-13 eV/Å`). These DFT differences include electronic redistribution and polarization after removing the other phase; they are effective force differences, not a rigorous two-body decomposition of the DFT Hamiltonian.

Vector RMSE is `sqrt(mean_i(sum_xyz((F_LAMMPS−F_DFT)^2)))`, taken over atoms in the named group. P95 and maximum are percentiles/maxima of each atom's vector-error norm. The interface subset is defined geometrically as nearest atom in the opposite phase ≤8 Å, matching the current Morse cutoff; 107 atoms qualify, 28 are outside.

## Force errors

| Comparison | Element/region | N | Vector RMSE (eV/Å) | P95 | Max | Max atom |
|---|---|---:|---:|---:|---:|---:|
| TSCBOP-only vs DFT TSC-only | Ti | 57 | 2.624203 | 4.096982 | 4.172368 | 68 |
| TSCBOP-only vs DFT TSC-only | Si | 19 | 2.607155 | 3.273784 | 3.273835 | 106 |
| TSCBOP-only vs DFT TSC-only | C | 38 | 2.443264 | 3.109278 | 3.235305 | 72 |
| TSCBOP-only vs DFT TSC-only | all Ti/Si/C | 114 | **2.562410** | 3.906345 | 4.172368 | 68 |
| EAM-only vs DFT Ag-only | Ag | 21 | **0.056485** | 0.085432 | 0.085457 | 130 |
| Morse vs effective DFT cross | Ag | 21 | 2.149582 | 2.377281 | 2.600507 | 129 |
| Morse vs effective DFT cross | Ti | 57 | 0.331636 | 0.731828 | 0.849349 | 2 |
| Morse vs effective DFT cross | Si | 19 | 2.742213 | 2.953676 | 2.953679 | 76 |
| Morse vs effective DFT cross | C | 38 | 0.094648 | 0.159219 | 0.159271 | 5 |
| Morse vs effective DFT cross | interface atoms ≤8 Å | 107 | **1.517760** | 2.897280 | 2.953679 | 76 |
| Morse vs effective DFT cross | atoms >8 Å | 28 | 0.034959 | 0.059507 | 0.059812 | 102 |
| Morse vs effective DFT cross | all 135 atoms | 135 | **1.351321** | 2.861368 | 2.953679 | 76 |

The original full hybrid potential remains at **2.254673 eV/Å** vector RMSE against the full DFT forces (reported previously as 2.2547 eV/Å); its 95th-percentile error is 3.482461 eV/Å and maximum is 3.519362 eV/Å at Ti atom 38. It is not valid to add the component RMSE values because vector errors can reinforce or cancel.

## Spatial interpretation

The largest TSCBOP-only discrepancies are Ti atoms in the substrate layer near z≈6.79 Å. For example, Ti 68 has 4.172 eV/Å phase-only error while its Morse cross-force error is below the TSCBOP discrepancy scale. C atoms also show large phase-only errors (C RMSE 2.443 eV/Å) but only 0.095 eV/Å effective-cross RMSE. This localizes the dominant Ti/C mismatch to the skeleton-force contribution on this geometry.

The largest effective-cross discrepancies are at the Si-terminated interface: Si 76, 64 and 112 each have about 2.954 eV/Å cross error at nearest Ag distance 2.523 Å and z≈8.88 Å. Ag atoms also have 2.150 eV/Å cross RMSE although the EAM-only Ag error is just 0.0565 eV/Å. In contrast, the Morse-vs-effective-cross residual outside its 8 Å range is small here (0.0350 eV/Å RMSE).

Per-atom forces, component errors, positions and distances are in `DFT_LAMMPS_component_force_errors.csv`; the top-error locations are also tabulated in `_SPATIAL_ERROR_TABLE.md`.

## Interpretation limits

- A single relaxed interface geometry is enough to identify where the present model's force error sits for this frozen configuration; it is not enough to validate a potential over other registries, gaps, temperatures or strains.
- The phase-only DFT structures are neutral slabs at the unchanged full cell and frozen coordinates. Their forces include surface and strain response of each isolated phase.
- The Ag₂₁ component is neutral with an odd valence count and uses the same non-spin-polarized FD-smearing setup as the full run. A spin-polarized Ag-cluster check was not part of this task.
- No absolute total energies of the different stoichiometries were subtracted to claim a binding energy or energy accuracy.

