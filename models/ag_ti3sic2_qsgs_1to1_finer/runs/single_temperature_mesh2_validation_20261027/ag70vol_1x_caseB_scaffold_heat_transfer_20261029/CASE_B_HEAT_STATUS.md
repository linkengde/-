# Ag70 1× Case B duration completed

## Run outcome

- Status: `AG70_CASE_B_HEAT_SEGMENTED_RUN_COMPLETED`
- Source on during run: `True`
- Stop reason: `duration_completed`
- End time: `2 ps`; completed segments: 20.
- Final valid checkpoint: `checkpoints/caseB_ag70_1x_heat_final_2ps.npz`
- Final Ag peak/mean temperature: 769.1222 / 385.2151 K.
- Final Ti₃SiC₂ peak/mean temperature: 743.7482 / 386.3340 K.
- Final Ag liquid volume fraction: 0; largest sampled connected cluster: 0 nm³.
- Final system enthalpy: 3.102717892475e-15 J; source ledger: 3.102717892477e-15 J; residual: -2.0964e-27 J.

## Sampled segment history

`heat_metrics.csv` records phase-resolved temperature, Ag liquid fraction, connected-cluster metrics, phase enthalpy, source energy, and energy residual at completed segment boundaries. Connectivity persistence is sampled, not a proof of hydrodynamic or mechanical stability.

## Interpretation limits

The calculation uses first-pass candidate properties, ideal Ag/Ti₃SiC₂ contact, and an insulated outer boundary. It includes no finite interface resistance, external heat sink, melt flow, vaporization, thermal expansion, or mechanical skeleton response. A source-off Case B result describes internal redistribution in the composite domain, not heat rejection into an electrode or the environment.
