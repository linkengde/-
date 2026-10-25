# QSGS heat-localization screen

## Scope

This screen addresses the 10 ps result in which the temperature response covered the small periodic QSGS box before Ag melted. It compares the existing volume-fraction electron-conduction links with links weighted by exact Ag contact area on each projected face, then repeats the same QSGS unit cell 3×3 in periodic x/y to test the thermal-box-size effect. A 5 Å calculation checks whether the 10 Å result is resolution dependent.

These are continuum two-temperature/enthalpy calculations, not TTM-MD. Case A is used: the Ti₃SiC₂ skeleton is geometrical only, there is no bottom sink, x/y are periodic, and z is insulated. The tiled case repeats the same saved QSGS unit cell; it does not create an independent larger random microstructure.

## Method

The heat source and material parameters were held fixed at the existing definitions: $I_0=5\ \mathrm{eV/(ps\,Å^2)}$, $R=30\ \mathrm{Å}$, $l_{skin}=20\ \mathrm{Å}$, spatial profile $e^{-r^2/R^2}e^{-d/l_{skin}}$, 0.3125 fs outer step, and 0.1 ps heating. Three-dimensional fields were saved at 0.001, 0.01, 0.05, and 0.1 ps. The source was not renormalized after Ag masking; the grid-integrated power is reported for each grid.

For the geometry-corrected cases, each coarse face conductance is multiplied by the area fraction of that face where the original binary QSGS mask has Ag on both sides. This gives a symmetric conservative electron-diffusion operator and suppresses links across faces with no source-geometry Ag contact. It is a finite-volume geometry correction, not an interface resistance or a microscopic transport law. The existing Ag volume fraction still supplies each cell's electronic energy capacity.

The expanded box repeats the 10 Å QSGS grid three times in each periodic lateral direction. The source remains centered over the middle copy, leaving the surrounding copies as lateral thermal context. The central copy is separately measured in addition to the full 3×3 box.

## Results at 0.1 ps

Peak and weighted metrics exclude projected Ag fractions below $10^{-12}$. Those fractions are floating-point remnants with negligible total Ag volume (below $3.6\times10^{-17}$ of the Ag weight in every case). The unfiltered numerical maximum in the geometry-corrected runs reached 5227 K in a cell with Ag fraction $1.11\times10^{-16}$; it is not a physical hotspot. The archived pilot fields were written just before this cleanup; their analysis applies the cutoff. The solver now clips such round-off-only phase slivers to zero when loading projected fields, preventing the diagnostic artifact in future runs.

| Case | Source power (mW) | Ag-weighted (T_e) mean (K) | Valid (T_e) max (K) | Ag volume above half-rise threshold | Hotspot (2r_{rms}) (Å) | Max Ag lattice (T_l) (K) | Ag liquid fraction |
|---|---:|---:|---:|---:|---:|---:|---:|
| 10 Å, original box, volume-fraction links | 0.863925 | 724.75 | 1063.80 | 100.0% | 246.8 | 300.62 | 0 |
| 10 Å, original box, exact QSGS face contacts | 0.863925 | 724.67 | 1130.56 | 55.0% | 196.8 | 300.69 | 0 |
| 10 Å, 3×3 repeat box, volume-fraction links | 0.863925 | 367.65 | 850.63 | 0.141% of full box; 1.27% of central tile | 66.9 | 300.51 | 0 |
| 10 Å, 3×3 repeat box, exact QSGS face contacts | 0.863925 | 367.47 | 929.96 | 0.103% of full box; 0.925% of central tile | 62.2 | 300.60 | 0 |
| 5 Å, original box, exact QSGS face contacts | 0.880833 | 730.19 | 1242.32 | 6.85% | 106.7 | 300.81 | 0 |

The 5 Å source integral is 1.96% above the 10 Å integral because the unchanged peak source is sampled on a finer cell-centered grid. No power renormalization was applied. The energy-balance residual magnitude is below $5\times10^{-31}$ J in all five cases. The maximum lattice temperature stays near 301 K, far below the Ag melting point; no liquid Ag forms in this 0.1 ps screen.

## Interpretation

The original 30 nm periodic box has a strong finite-domain effect: with legacy links, its whole Ag phase exceeds the run-relative half-rise threshold at 0.1 ps, while the mean electron temperature in the central 30 nm tile of the 3×3 case is about 425 K. The geometry correction changes the original-box response but does not eliminate broad heating there. At 3×3 size, corrected and uncorrected central-tile half-rise fractions are both about 1%, so the periodic-box effect is more prominent than this first-order face correction in this screen.

The 5 Å corrected result differs materially from the 10 Å corrected result: the half-rise fraction falls from 55% to 6.85%, and the estimated hotspot diameter also changes. This means the current localization metric is not spatially converged. The roughly 2% source-power difference also contributes to the resolution comparison. The 10 Å result must not be used to claim a converged local melt pool, and the 5 Å field is still only a short transient, not a melt-pool or solidification result.

The next decision should be to retain the exact-face operator and round-off cleanup, improve source cell averaging or report the grid-integrated power with each resolution, then perform a shorter time-step and spatial convergence check on 5 Å before increasing heating duration. Do not interpret this screen as evidence for melting or extend the old 10 Å run to 20–100 ps under unchanged conditions.

## Reproduction

From the model directory, regenerate face-contact fields and rerun a case with:

```bash
python phase_field_qsgs/qsgs_face_contacts.py --spacing 10A
python phase_field_qsgs/qsgs_face_contacts.py --spacing 5A
python phase_field_qsgs/qsgs_ttm_enthalpy.py --spacing 10A --boundary A --heat-ps 0.1 --cool-ps 0 --dt-fs 0.3125 --sample-ps 0.01 --snapshot-ps 0.001,0.01,0.05,0.1 --face-contacts phase_field_qsgs/qsgs_ag_face_contacts_10A.npz --output phase_field_qsgs/qsgs_localization_base_exact_contacts.csv
python phase_field_qsgs/qsgs_ttm_enthalpy.py --spacing 10A --boundary A --heat-ps 0.1 --cool-ps 0 --dt-fs 0.3125 --sample-ps 0.01 --snapshot-ps 0.001,0.01,0.05,0.1 --tile-xy 3 --face-contacts phase_field_qsgs/qsgs_ag_face_contacts_10A.npz --output phase_field_qsgs/qsgs_localization_3x3_exact_contacts.csv
python phase_field_qsgs/qsgs_ttm_enthalpy.py --spacing 5A --boundary A --heat-ps 0.1 --cool-ps 0 --dt-fs 0.3125 --sample-ps 0.01 --snapshot-ps 0.001,0.01,0.05,0.1 --face-contacts phase_field_qsgs/qsgs_ag_face_contacts_5A.npz --output phase_field_qsgs/qsgs_localization_5A_exact_contacts.csv
python phase_field_qsgs/analyze_heat_localization_screens.py
python phase_field_qsgs/plot_heat_localization_comparison.py
```

The four 10 Å comparison runs and the 5 Å check are archived as CSV histories, JSON summaries, and 3-D field snapshots. Reanalyzed, round-off-filtered metrics are in `qsgs_localization_screen_metrics.json`; the comparison image is `qsgs_heat_localization_screen_0p1ps.png`.
