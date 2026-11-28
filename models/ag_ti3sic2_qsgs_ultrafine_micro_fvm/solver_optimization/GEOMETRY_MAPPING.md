# QGS mapping and connectivity checks

The original seed 1105 150³ Ag fraction was conservatively volume-overlap remapped onto 50³ and 100³ grids. The input has mean Ag fraction 0.70. The report checks mixed-cell fraction, thresholded phase connectivity using 6/18/26 neighbors, XYZ spanning, and the source-rate volume integral. Connectivity here is computed after thresholding a mixed control volume at phase fraction 0.5; it does not prove a physical Ag channel at subcell scale.

| Grid | Mixed-cell fraction | Ag components (6/18/26) | Ag spans XYZ | Local source integral error |
|---:|---:|---:|---|---:|
|50³|28.532%|2 / 2 / 2|yes|4.501%|
|100³|14.4072%|3 / 1 / 1|yes|1.8704%|
|150³ native|0%|2 / 2 / 2|yes|0%|

Ti₃SiC₂ also spans all axes at each tested resolution. Its thresholded component counts change with neighborhood choice, reinforcing that a binary connected-component count is a diagnostic of the mapped voxel geometry, not a resolved transport topology. The source-integral error comes from remapping the volumetric heating mask and is not a thermal energy residual; it should be accounted for when comparing resolutions.

At the 100³ safe-run heated endpoint, Ag-only bulk liquid fraction was 3.33011×10⁻⁴ and the actual Ag-fraction-weighted liquid volume was 2.33108×10⁻¹⁹ m³. With Ag fraction threshold 0.5, the thresholded geometric liquid volume decreased from 3.27×10⁻¹⁹ m³ at liquid threshold 0.25 to 2.45×10⁻¹⁹ m³ at 0.5 and 1.27×10⁻¹⁹ m³ at 0.75. The component stayed at one under those thresholds. These geometric threshold volumes differ from the weighted Ag liquid volume by design; do not substitute one for the other.

Full connectivity and threshold data are in `results/QSGS_MAPPING_CONNECTIVITY.csv` and `results/checkpoints/qgs_n100_numba_ew_safe_t4_full/liquid_threshold_sensitivity_heated.csv`.
