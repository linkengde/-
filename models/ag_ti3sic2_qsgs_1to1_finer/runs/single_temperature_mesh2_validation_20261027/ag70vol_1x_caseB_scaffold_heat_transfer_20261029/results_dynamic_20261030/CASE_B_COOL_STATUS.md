# Ag70 1× Case B all ag solid

## Run outcome

- Status: `AG70_CASE_B_COOL_SEGMENTED_RUN_COMPLETED`
- Source on during run: `False`
- Stop reason: `all_ag_solid`
- End time: `14.666 ps`; completed segments: 16.
- Final valid checkpoint: `checkpoints/caseB_ag70_1x_cool_final_14p666ps.npz`
- Final Ag peak/mean temperature: 1209.7612 / 879.6238 K.
- Final Ti₃SiC₂ peak/mean temperature: 1156.3067 / 881.1116 K.
- Final Ag liquid volume fraction: 0; largest sampled connected cluster: 0 nm³.
- Final system enthalpy: 2.250401287426e-14 J; source ledger: 2.250401287410e-14 J; residual: 1.5610e-25 J.

## Sampled segment history

`cool_metrics.csv` records phase-resolved temperature, Ag liquid fraction, connected-cluster metrics, phase enthalpy, source energy, and energy residual at completed segment boundaries. Connectivity persistence is sampled, not a proof of hydrodynamic or mechanical stability.

## Interpretation limits

The calculation uses first-pass candidate properties, ideal Ag/Ti₃SiC₂ contact, and an insulated outer boundary. It includes no finite interface resistance, external heat sink, melt flow, vaporization, thermal expansion, or mechanical skeleton response. A source-off Case B result describes internal redistribution in the composite domain, not heat rejection into an electrode or the environment.
