# Ag70 1× Case B status through 2 ps

## Scope

This is a fresh 300 K run of the Ag70, 100×100×50 native QSGS geometry with the single-temperature enthalpy solver's existing Case B: Ag and Ti₃SiC₂ both store and conduct heat, the interface is ideal contact, only Ag receives the original Gaussian source, x/y are periodic, and z is insulated. No Case A checkpoint was converted or reused, and no 2×, TTM, or MD job was started here.

The geometry SHA-256 is `03ef47b1dd534b3215cf03f42c8cf7e887dcdfcc185b33bab7453cd838c1faac`; its grid fractions are 69.9996 vol% Ag and 30.0004 vol% Ti₃SiC₂. Integrated source power is 1.551358946238 mW. The 300 K explicit stability step is 0.0404211 fs.

## Preflight

`CASE_B_MODEL_VALIDATION.md` and `caseB_preflight.json` record the passed checks: 300 K static source-off evolution, uniform 500 K no-flux equilibrium, nonuniform Ag-to-scaffold heat transfer, energy conservation, explicit time-step rejection, Ag-only source support, and rejection above the Ti₃SiC₂ table ceiling of 1473.15 K. Four independent regression tests pass in `test_ag70_1x_caseB.py`.

The nonuniform conduction test warms the scaffold from 300 K to 332.040 K while the Ag weighted mean cools from 600 K to 595.654 K. System energy changes by only 2.04×10⁻¹⁶ relative.

## Heating result

Heating began from a fresh 300 K Case B state. It is saved in 0.05 ps segments; complete restart checkpoints exist at each segment boundary. At 2 ps:

| Metric | Case B result |
| --- | ---: |
| Ag peak temperature | 769.122 K |
| Ti₃SiC₂ peak temperature | 743.748 K |
| Ag weighted mean temperature | 385.215 K |
| Ti₃SiC₂ weighted mean temperature | 386.334 K |
| Ag liquid volume fraction | 0 |
| Ag sensible enthalpy | 2.1420745×10⁻¹⁵ J |
| Ag latent enthalpy | 0 J |
| Ti₃SiC₂ stored enthalpy | 9.6064342×10⁻¹⁶ J |
| Total enthalpy / cumulative source input | 3.1027179×10⁻¹⁵ J / 3.1027179×10⁻¹⁵ J |
| Energy residual | −2.10×10⁻²⁷ J |

No liquid has appeared by 2 ps, so there is no Case B pool to cool and no heat-off equilibrium calculation is warranted yet. The maximum scaffold temperature is safely below its current property ceiling. Case B heating must continue from `checkpoints/caseB_ag70_1x_heat_final_2ps.npz` if the melt-onset window is to be found.

## Matched 1 ps comparison with Case A

`caseA_caseB_matched_1ps.csv` and `.json` restore and compare the original Case A and new Case B checkpoints at the same 1 ps time. They have matching Ag70 geometry hashes, initial temperature, source definition, and cumulative source energy. Case B's Ag peak is 193.18 K lower (21.32%) than Case A. Ag mean temperature is 19.11 K lower. Both have zero liquid at this time. Case B has 4.8385346×10⁻¹⁶ J stored in the scaffold, 31.2% of the total input energy.

This is the effect of the Case A/Case B thermal assumptions in this geometry and source setup. It is not evidence that either composition is more weld-resistant in an experiment.

## Runtime and next step

The completed 0–2 ps path consumed about 49.9 minutes of wall time, averaging roughly 25 minutes per simulated ps on this machine. This measured cost should guide the next staged extension. Continue from the 2 ps checkpoint in bounded segments, retain every completed segment checkpoint, and stop if Ag reaches 2300 K, Ti₃SiC₂ reaches 1473.15 K, or the energy ledger fails. Do not increase source power to force melting.

Candidate material properties remain incompletely source-verified. Case B still assumes zero interface resistance, an insulated outer boundary, a common temperature field, and no melt flow, vaporization, expansion, or mechanical skeleton response.
