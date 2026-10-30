# Ag70 1× Case B ti3sic2 candidate property limit reached

## Run outcome

- Status: `AG70_CASE_B_HEAT_SEGMENTED_RUN_COMPLETED`
- Source on during run: `True`
- Stop reason: `ti3sic2_candidate_property_limit_reached`
- End time: `14.506 ps`; completed segments: 132.
- Final valid checkpoint: `checkpoints/caseB_ag70_1x_heat_final_14p506ps.npz`
- Final Ag peak/mean temperature: 1472.8031 / 878.5024 K.
- Final Ti₃SiC₂ peak/mean temperature: 1431.5269 / 880.1530 K.
- Final Ag liquid fraction (within Ag): 0.0033270131; latest CSV sample has 2 connected clusters, with the largest containing 1,161 voxels (33.3113 nm³), centroid (145.368, 150.577, 150.252) Å. The previous summary field showed 0 nm³ and has been corrected from the segment CSV.
- Final system enthalpy: 2.250401287426e-14 J; source ledger: 2.250401287410e-14 J; residual: 1.5610e-25 J.

## Sampled segment history

`heat_metrics.csv` records phase-resolved temperature, Ag liquid fraction, connected-cluster metrics, phase enthalpy, source energy, and energy residual at completed segment boundaries. Connectivity persistence is sampled, not a proof of hydrodynamic or mechanical stability.

## Interpretation limits

The calculation uses first-pass candidate properties, ideal Ag/Ti₃SiC₂ contact, and an insulated outer boundary. It includes no finite interface resistance, external heat sink, melt flow, vaporization, thermal expansion, or mechanical skeleton response. A source-off Case B result describes internal redistribution in the composite domain, not heat rejection into an electrode or the environment.
