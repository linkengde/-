# P1P final decision: heating, evaporation, and permanent-pit assessment

Date: 2026-11-07
Branch: `codex/ag70-p1p-evaporation`
Base commit: `6aecbac46cf7dda59bea72767983ac86b85d2142`

## Decision

**No resolvable permanent crater was produced in the completed A or B runs. C is incomplete, so its post-cooling crater status is unknown.** This is an exploratory result for the present continuum VOF/enthalpy model; it is not a conclusion about real Ag–Ti₃SiC₂ contacts.

The final solidified A and B surface fields contain no positive center depression at the current `0.301857 nm` grid spacing. Their center-averaged displacement is slightly upward: A `−0.023 nm`, B `−0.075 nm` (the CSV convention is positive downward). The largest possible measured downward displacement in the center zone is `0 nm`; surface heights are voxel-quantized. Thus neither completed case meets the requested resolved-crater criterion.

## Full A/B/C results

All campaigns began from the same verified 300 K initial field (SHA256 `abd706826e90a174c6eb37322b60637b54b59f1bab9fab18b56ab06801b61e73`), then continued from case-specific checkpoints after the short-run phase. The continuation checkpoint hashes are listed in the comparison CSV. All used a 4 mW maximum composite heat source and 3× Marangoni. A used zero arc pressure and evaporation off; B used 501 MPa and evaporation off; C used 501 MPa and evaporation on with `p_Ag,v = 0 Pa`. The 501 MPa and 3× settings are exploratory inputs, not calibrated arc conditions.

| Case | Peak liquid Ag | Largest 6-neighbor cluster | Hold/cooling result | Final center surface | Decision |
|---|---:|---:|---|---:|---|
| A: no pressure, no evaporation | 93.399 nm³ | 82.985 nm³; 7 components | Completed 5.020 ps hold; fully solidified at 12.86318 ps | −0.023 nm | No resolved pit |
| B: 501 MPa, no evaporation | 85.878 nm³ | 75.287 nm³; 11 components | Hold stopped after 2.640 ps by pre-step property guard; then fully solidified at 9.52330 ps | −0.075 nm | No resolved pit; shortened hold |
| C: 501 MPa, evaporation | 85.878 nm³ before stop | 75.287 nm³; 11 components | Stopped at 9.30331 ps by P1P property/composition failure; no cooling or solidified field | Not available | Incomplete; pit unknown |

At the peak/hold-start snapshot, the center 0–3 nm radial Ag volume had fallen by `0.347 nm³` in B/C and the 3–5 nm annulus gained `0.296 nm³`; outer annuli gained the remainder to within the global mass drift. This is evidence of some redistribution in this numerical model, but the final B surface still has no downward center displacement. A's corresponding center loss was `0.061 nm³`.

Connectivity shows 7 thresholded components in A and 11 in B/C at the peak, but the per-component table shows the secondary components total only about `0.3 nm³` and each spans just 1–5 cells. They are unresolved voxel-scale remnants, not multiple substantial melt pools. The resolved field is dominated by one main pool.

By B's final solidified snapshot, center 0–3 nm Ag volume was down `0.723 nm³`; the 3–5, 5–8, and >8 nm zones gained `0.618`, `0.062`, and `0.015 nm³`. Total Ag volume drift was `−0.000688%`. A's center loss was `0.158 nm³`, largely offset by `+0.141 nm³` in the 3–5 nm annulus; total drift was `−0.0000024%`. These radial changes do not establish a crater: the surface height analysis shows no downward center depression, and the measured heights are quantized at one grid interval.

The completed energy ledgers close closely: A input `5.0520×10⁻¹⁴ J`, bottom sink `4.4841×10⁻¹⁴ J`, residual `−8.64×10⁻²¹ J`; B input `3.7154×10⁻¹⁴ J`, bottom sink `3.1587×10⁻¹⁴ J`, residual `−6.49×10⁻²¹ J`. C has no final ledger because it stopped before a completed final thermal step/cooling phase.

## Stops and physical interpretation

- **B:** the guard logged `Tmax_skeleton_K=5000 K`, but 5000 K is the enthalpy-inversion upper sentinel, not a measured field temperature. The last valid field/history had a global maximum of `1769.8 K` and skeleton maximum `1706.8 K`. The pre-step guard shut off heating and pressure conservatively; this run completed cooling but did not complete its planned 5 ps hold. The 5000 K value must not be reported as a physical temperature.
- **C:** the original run stopped before the thermal update at `9.30330989 ps` after 7 interface cells failed validation against the `1234–2436 K` property table; no extrapolation was used. An isolated instrumented replay from the verified `9.00332084 ps` restart returned `P1P_ERR_RANGE` (status `3`) in 5 cells at `9.28035320 ps`. Their reconstructed interface temperatures were `4925–5320 K`, while stored material fractions were extremely small (`f ≤ 1.17×10⁻⁹`, `Ag ≤ 5.06×10⁻¹¹`, `Ti₃SiC₂ ≤ 1.12×10⁻⁹`). With almost all cell volume classified as gas, `Hadv` inverted through the model's gas heat capacity alone implies roughly `4293–4911 K` before the interface-gradient extrapolation. This supports a gas-dominated enthalpy / vanishing-Ag-interface artifact as the immediate mechanism, but does **not** identify how the anomalous enthalpy reached those cells. The replay stopped about `0.023 ps` earlier and found 5 cells rather than the original 7, so it is diagnostic evidence, not an exact reproduction. There is no final cooled field and no permanent-pit conclusion for C.
- C's accumulated evaporated mass was `1.0911×10⁻³⁰ kg`, about `6.09×10⁻⁶` Ag atomic-mass equivalents. This is a continuum flux integral, not a literal fractional atom, and is negligible at this scale. C and B were indistinguishable in pool volume and velocity through the last valid state to displayed precision.
- B's maximum effective imposed pressure was `500.4 MPa` with integrated normal force magnitude `1.945×10⁻⁸ N`. This confirms the exploratory load entered the solver; it does not prove the value represents a real arc.

## Interpretation

The completed evidence supports: **A and B each had one resolved dominant liquid pool plus disconnected one-to-five-cell remnants; the remnants are too small to identify as physically resolved secondary pools. The B pressure case showed modest center-to-annulus Ag redistribution, but no resolvable permanent crater after cooling.** C cannot be compared after solidification because its run stopped before cooling. Evaporation at the tested conditions is far too small to explain a pit.

Do not claim that pressure or Marangoni has no physical role, or that real nanoscale contacts cannot form pits. The present grid, voxelized surface representation, constitutive assumptions, and incomplete C run limit that inference. The immediate C stop is confirmed as a property-table range rejection in an instrumented replay, but the anomalous reconstructed interface temperature and tiny cell fractions still need root-cause analysis before any further evaporation-enabled full case. No new long calculation was started here.

## Files

- `P1P_ABC_FULL_FINAL_COMPARISON.csv`: case settings, stop/completion status, peak pool, mass/energy totals, evaporation, and pressure load.
- `P1P_ABC_FULL_RADIAL_MASS.csv`: initial-to-final radial Ag redistribution for completed cases.
- `P1P_ABC_FULL_PEAK_RADIAL_MASS.csv`: initial-to-peak/hold-start radial redistribution.
- `P1P_ABC_FULL_SURFACE_SHAPE.csv`: final surface height displacement for completed cases.
- `P1P_ABC_FULL_PEAK_CLUSTERS.csv`, `P1P_ABC_FULL_LIQUID_CLUSTERS.csv`: 6-neighbor connectivity at thresholds 0.1, 0.5, and 0.9.
- `P1P_ABC_FULL_PEAK_COMPONENTS.csv`: per-component volume and bounding coordinates.
- `P1P_C_FAILURE_CELL_DIAGNOSTIC.csv`: per-cell status and local state from the isolated instrumented C replay.
- `P1P_ABC_FULL_HOLD_START_3D.png` and `P1P_ABC_FULL_HOLD_START_ZOOM.png`: phase-colored A/B/C visualizations.
- `P1P_ABC_LIQUID_COMPONENTS.png` and `P1P_ABC_SECONDARY_POOLS.png`: pool components shown separately.
- `P1P_COMBINED_DIAGNOSIS.md`: interface-temperature, replay, archival energy/mass, and full-run diagnostics.
