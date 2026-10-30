# Ag70 Case A/Case B matched comparison and interface-parameter screen

## Matched-time comparison

The checkpoints use the same Ag70 QSGS geometry hash, 100×100×50 native grid, 1× refinement, Gaussian source power (1.551358946 mW integrated over Ag), 30 Å lateral radius, and 20 Å depth e-fold. Case A is at 14.510 ps and Case B at 14.506 ps (4 fs apart). Their accumulated source energies differ by only 0.028%; this is a close time-and-energy comparison, not a cross-composition comparison.

| Metric | Case A: thermally inactive scaffold | Case B: scaffold storage and conduction |
|---|---:|---:|
| Ag peak temperature | 2120.9 K | 1472.8 K |
| Ag mean temperature | 1100.4 K | 878.5 K |
| Ti₃SiC₂ peak temperature | 300 K | 1431.5 K |
| Ag liquid fraction, normalized by Ag volume | 4.397% | 0.3327% |
| Largest connected cluster (`f_l ≥ 0.5`) | 15,367 voxels; 440.91 nm³ | 1,161 voxels; 33.31 nm³ |
| Total source energy | 2.25102×10⁻¹⁴ J | 2.25040×10⁻¹⁴ J |
| Relative energy error | 7.15×10⁻¹² | 6.94×10⁻¹² |

At these nearly matched conditions, the Case A largest cluster is 13.2 times the Case B cluster volume. Case B stores 6.990×10⁻¹⁵ J in the Ti₃SiC₂ phase by heat-off, about 31.1% of system enthalpy. This is the modeled reason for the smaller, cooler Ag melt region. It is an effect of the thermal assumptions, not a visualization-scale difference.

Case A is a thermal-localization limiting scenario: its scaffold has zero heat capacity and conductivity. It should not be described as a real insulating Ti₃SiC₂ skeleton. Case B is closer to a two-phase composite response but still assumes perfect Ag/Ti₃SiC₂ contact and insulated external boundaries.

## What the existing property record allows

The repository's `MATERIAL_SOURCE_AUDIT.md` says the exact primary records behind the present Ti₃SiC₂ Cp/k anchors are not archived, and the Ag liquid Cp/k anchors are also incomplete or tentative. Case B stopped at 14.506 ps because the next accepted step would leave the Ti₃SiC₂ candidate table range of 300–1473.15 K. Do not extend heating by extrapolating that table or by increasing source power. The accepted heat-off state remains inside the table.

## Interface-resistance scale check

`compare_caseA_caseB_matched_14p51ps.py` recreates the matched comparison in `caseA_caseB_matched_14p51ps.csv`; the 14.510 ps Case A checkpoint is archived beside the Case B checkpoints for this comparison. The same-scale 3D image is `caseA_caseB_melt_pool_same_scale_14p51ps.png`.

`interface_contact_resistance_sweep.py` and `interface_conductance_resistance_sweep.csv` compare the areal interface resistance `1/G_K` with the sum of the two cell-centred bulk half-cell resistances,

`R_bulk = Δx/(2 k_Ag) + Δx/(2 k_Ti3SiC2)`.

At the current ~0.30 nm voxel pitch and the candidate conductivity anchors, these resistances are comparable near `G_K ≈ 1.8–2.3×10¹¹ W m⁻² K⁻¹`, depending on axis and temperature anchor. A logarithmic test sweep from 10⁸ to 10¹³ W m⁻² K⁻¹ spans interface-dominated to near-perfect-contact face transfer. This sweep is a numerical bracketing range only; it is **not** a measured or recommended physical range. The existing finite-interface kernel is a fixed-capacity conduction prototype without Ag phase change and has not been coupled to the Case B enthalpy solver. Its results must not be presented as melt-pool predictions.

## Next decision gate

Before using a finite `G_K` in a physical Case B melt calculation, locate an interface conductance range for the actual Ag/Ti₃SiC₂ preparation or set a clearly labelled, non-calibrated sensitivity bracket. Then add the interface law to the enthalpy model with matched energy accounting and verify `G_K→0` and `G_K→∞` before any longer heating run. Independently, replace the present candidate Cp/k inputs with source-backed high-temperature data or stop before the table ceiling, as done here. No long transient was run for this screen.
