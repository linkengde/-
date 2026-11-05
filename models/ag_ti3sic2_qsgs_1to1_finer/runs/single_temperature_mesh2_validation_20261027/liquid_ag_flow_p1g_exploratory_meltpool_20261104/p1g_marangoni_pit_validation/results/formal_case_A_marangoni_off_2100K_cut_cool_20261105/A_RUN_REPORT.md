# Formal P1G Case A: Marangoni off

- Run: `formal_case_A_marangoni_off_2100K_cut_cool_20261105`
- Solver: isolated P1G CLSVOF/enthalpy prototype; no production solver changes.
- Initial state: 300 K, common Ag70 QSGS geometry, frozen 64³ grid, Δx = 0.301857 nm.
- Heat source: normalized Ag + Ti₃SiC₂ Gaussian/depth-exponential source, fixed at 8.931899 mW. Marangoni: off. No evaporation or ejection.
- Heat-off policy: actual Ti₃SiC₂ maximum reaches 2100 K, or full-power input and bottom sink differ by ≤0.01 mW for 10 consecutive steps. Hard property ceiling 2300 K. The run stopped by the temperature condition.

## Results

- First liquid Ag: 0.132666250 ps; 0.0113234 nm³.
- Maximum total liquid Ag: 121.852732 nm³ at 0.624017255 ps.
- Maximum connected liquid cluster: 110.090850 nm³, 4035 cells, at 0.624017255 ps; bounding dimensions X/Y/Z = 10.2631545 / 10.565012 / 3.32043233 nm.
- Heat-off: HEAT_OFF_AT_TARGET t_ps=0.624017255 Tmax_skeleton_K=2117.73401 Tmax_Ag_K=2093.287 liquid_Ag_nm3=121.852732 cut_K=2100
- Complete solidification: DONE_SOLIDIFIED t_ps=0.861280016 Tmax=1189.41287 K Tmax_Ag=1056.42561 K Tmax_skeleton=1062.78488 K
- Maximum surface pit depth in the VOF column-height audit: 0 nm; maximum rim rise: 0 nm.
- Peak velocity: 0 m/s. With Marangoni disabled, this control remained at zero in the saved history.
- Ag volume drift maximum: 0%.
- Maximum relative energy residual (residual / stored enthalpy): 2.87204e-08.
- Stored time-series energy fields at the final step: H = 4.43951510072e-15 J, cumulative input = 5.7327953547e-15 J, cumulative ideal-sink removal = 1.29328025937e-15 J.
- Time-resolved fields: `thermal_flow_fields.h5`, 19 frames, 64×64×64, finite temperature field; raw snapshots retained as `.bin.gz`. Time series and audits: `p1g_history.csv`, `P1G_AG_MASS_BALANCE.csv`, `P1G_ENERGY_BALANCE.csv`, `P1G_PIT_GEOMETRY.csv`, `p1g_peak_pool_dimensions.csv`.

This is the no-Marangoni control only. It does not establish the Marangoni effect; Case B is required for that comparison. The ideal bottom sink and the exploratory microstructure/VOF coupling remain model assumptions, not an external heat-transfer or atomistic prediction.
