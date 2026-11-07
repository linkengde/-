# P1P combined diagnosis: interface temperature, replay, energy, and pit mechanism

Date: 2026-11-07
Branch: `codex/ag70-p1p-evaporation`
Base commit: `6aecbac46cf7dda59bea72767983ac86b85d2142`
Production solver, source geometry, temperature reconstruction, and force multipliers were not modified for this post-processing.

## Executive result

The production interface-temperature algorithm remains the cross-interface centered-difference extrapolation. Replaying the exact same checkpoint is deterministic; the previous discrepancy came from comparing a surface-source calculation before the thermal update with a checkpoint saved after the update. The archived 3.0128 ps case and full A/B/C cases show liquid Ag in several connected components, but no completed case produced a resolvable permanent center crater. A and B cooled to fully solidified states; C stopped before cooling due to a P1P property/composition failure and therefore has no final crater assessment.

The evaporation module's integrated mass remains negligible at this scale. No evaporation coefficient or pressure was increased.

## Interface temperature reconstruction

The production method was retained. For 403 fit-eligible interface cells at the archived state, comparison to liquid-side one-sided linear fits gave an area-weighted absolute temperature difference of `10.389 K` and a maximum local difference of `114.870 K`. Recomputing Hertz–Knudsen flux on each method's reclassified active area changed total mass rate by `−6.7%`; using a common area gave `−10.8%`. The absolute rates were still extremely small.

Manufactured checks included planar, oblique-plane, spherical, and rough-interface fields. The one-sided method reproduced ideal planar/oblique linear profiles well, but was not uniformly more accurate on curved or rough interfaces. The synthetic fields do not establish which method is more accurate on an actual rough Ag/Ti₃SiC₂ interface. This supports keeping the production centered-difference algorithm for these runs.

Detailed data: `P1P_INTERFACE_TEMPERATURE_COMPARISON.csv`, `p1p_replay_interface_temperature_diagnostic.csv`, and `p1p_manufactured_interface_temperature.csv`.

## Checkpoint replay reproducibility

In the solver, `p1p_apply_surface_sources(dt)` computes and records the evaporation source before `thermal_step`; the thermal update and periodic checkpoint follow. Therefore an evaporation row and checkpoint with the same printed time can represent opposite sides of one update. At the archived label `3.0128048117 ps`, the recorded source row had interface area `2.82245×10⁻¹⁷ m²`, max flux `2.58734×10⁻³ kg m⁻² s⁻¹`, and max interface temperature `1479.762 K`. Re-evaluating the next source call from the checkpoint used the updated state: area `2.84478×10⁻¹⁷ m²`, max flux `2.24205×10⁻³ kg m⁻² s⁻¹`, and max interface temperature `1465.127 K`.

Two replays from the same checkpoint, code, and input produced byte-identical temperature, summary, and flux CSVs. The apparent discrepancy was event ordering, not nondeterministic replay. See `P1P_REPLAY_REPRODUCIBILITY.csv` and `P1P_REPLAY_REPEAT_HASHES.csv`.

## Archived 3.0128048 ps case

This was 300 K initial heating with 4 mW composite source, 3× tangential Marangoni, 501 MPa exploratory arc pressure, P1P evaporation enabled, and vapor partial pressure set to 0 Pa. It did not include a full cooling interval and cannot establish a permanent crater.

At 3.0128 ps: global peak temperature `1474.353 K`; liquid Ag volume `26.803 nm³`; 6-neighbor extraction gave 9 liquid components and a largest component of `20.960 nm³`. The center 5 nm mean height change was about `0.240 nm upward` (downward is positive), smaller than the `0.3019 nm` grid spacing. Center 0–3 nm Ag volume changed by `−1.94×10⁻²⁵ kg` equivalent, while the 3–5 nm annulus gained `1.15×10⁻²⁵ kg`; this is too small to claim meaningful radial transfer.

The energy ledger was: input `1.21232968×10⁻¹⁴ J`, bottom sink `6.63573176×10⁻¹⁵ J`, stored enthalpy `5.48756183×10⁻¹⁵ J`, residual `−3.22477×10⁻²¹ J`. Ag volume drift was `−0.0002146%`. Cumulative evaporation was `1.06864×10⁻³² kg`, about `5.97×10⁻⁸` Ag atomic-mass equivalents; that continuum integral is not a literal fraction of an atom.

See `P1P_ENERGY_LEDGER.csv`, `P1P_MASS_LEDGER.csv`, `P1P_POOL_CONNECTIVITY_HISTORY.csv`, and the `P1P_CHECKPOINT_*` files.

## A/B/C common-state comparison

All three campaigns started from the same verified 300 K initial field (SHA256 `abd706826e90a174c6eb37322b60637b54b59f1bab9fab18b56ab06801b61e73`), then continued from their own case-specific restart checkpoints after the short-run phase. The restart hashes are recorded separately in the comparison CSV because the prior force settings make those evolved checkpoints different. Conditions were maximum 4 mW composite heat source and 3× Marangoni. A used pressure off / evaporation off; B used 501 MPa / evaporation off; C used 501 MPa / evaporation on with `p_Ag,v=0 Pa`. The 501 MPa and 3× settings are exploratory assumptions, not calibrated physical conditions.

### Full-run table

| Case | Peak liquid Ag | Largest 6-neighbor cluster | Hold and final state | Result |
|---|---:|---:|---|---|
| A | `93.399 nm³` at `7.623 ps` | `82.985 nm³`, 7 components | 5.020 ps hold completed; solidified at `12.86318 ps` | No resolved pit |
| B | `85.878 nm³` at `6.663 ps` | `75.287 nm³`, 11 components | Hold stopped at `2.640 ps` by pre-step property guard; cooled/solidified at `9.52330 ps` | Shortened hold; no resolved pit |
| C | `85.878 nm³` through last valid peak at `6.663 ps` | `75.287 nm³`, 11 components | Stopped at `9.30331 ps`; no cooling/final field | Incomplete; pit unknown |

At the peak, the 6-face periodic threshold (`f × Ag_fraction × liquid_fraction ≥ 0.5`) gives 7 components for A and 11 for B/C. The component list shows one dominant cluster (`82.985 nm³` in A, `75.287 nm³` in B/C); all secondary components combined are only about `0.3 nm³` in each case, and consist of one-to-five-cell remnants. Thus the voxel criterion technically finds disconnected liquid cells, but it does **not** resolve multiple substantial pools separated by the skeleton. The physically readable result at this grid is one main pool plus tiny isolated threshold remnants, not several well-resolved secondary pools. B's maximum effective arc pressure was `500.4 MPa` and integrated normal force magnitude `1.945×10⁻⁸ N`. A's balanced hold source reached `3.94761 mW`. B/C were near `3.917 mW` at the last balance event (`9.163 ps`), where the logged source–sink difference was `−9.72×10⁻⁶ mW`.

At hold start, radial Ag redistribution relative to the common initial state was:

| Case | ΔAg, r=0–3 nm | ΔAg, 3–5 nm | ΔAg, 5–8 nm | ΔAg, >8 nm |
|---|---:|---:|---:|---:|
| A | −0.061 nm³ | +0.057 nm³ | +0.004 nm³ | +0.0002 nm³ |
| B | −0.347 nm³ | +0.296 nm³ | +0.020 nm³ | +0.004 nm³ |
| C | −0.347 nm³ | +0.296 nm³ | +0.020 nm³ | +0.004 nm³ |

B's completed final snapshot had center 0–3 nm Ag volume change `−0.723 nm³`; outer bands gained `+0.618`, `+0.062`, and `+0.015 nm³`. Its total Ag volume drift was `−0.000688%`. A's final center change was `−0.158 nm³`, with `+0.141 nm³` in the 3–5 nm band and total drift `−0.0000024%`. These are modest model-space transfers; they did not create a downward center surface displacement.

### Surface and evaporation decision

At the final solidified fields, the 0–3 nm mean downward displacement (positive means down) was `−0.023 nm` for A and `−0.075 nm` for B, i.e. slightly upward. No center column was lower than the initial surface in the measured center zone; the grid spacing is `0.301857 nm`. C has no solidified surface file. Consequently, A and B fail the resolvable permanent-crater criterion; C is unassessed.

Cumulative C evaporation reached `1.0911×10⁻³⁰ kg`, about `6.09×10⁻⁶` Ag atomic-mass equivalents. B and C's liquid volume and velocities matched through the last valid state to displayed precision. The current evaporation calculation therefore contributed no measurable morphology change in this run.

### Run-stop details

- B's `5000 K` guard report is the enthalpy inversion upper sentinel, not a physical cell temperature. The preceding valid stored state had global maximum `1769.8 K` and skeleton maximum `1706.8 K`. The guard conservatively shut off heat and pressure; cooling completed, but the planned 5 ps hold did not.
- C stopped before the thermal update at `9.30330989 ps`: the original run counted 7 interface cells rejected by the P1P property/composition validation for the `1234–2436 K` table; extrapolation was refused. An isolated instrumented replay from the verified `9.00332084 ps` restart returned `P1P_ERR_RANGE` (status `3`) for 5 cells at `9.28035320 ps`, with reconstructed interface temperatures from `4925` to `5320 K`. In those cells, `f ≤ 1.17×10⁻⁹`, `Ag ≤ 5.06×10⁻¹¹`, and `Ti₃SiC₂ ≤ 1.12×10⁻⁹`. Since the solver's gas volumetric heat capacity is `1206 J m⁻³ K⁻¹`, the gas-dominated estimate `300 + Hadv/1206` is already `4293–4911 K`; the free-surface gradient extrapolation raises this to the logged interface temperatures. All five cells pass the current minimum checks (`f > 10⁻¹²`, Ag coverage `>10⁻¹²`) even though their PLIC areas are just `4.45×10⁻⁸–3.86×10⁻⁷` of a nominal cell face. This supports the interpretation that tiny VOF/Ag traces in gas-dominated enthalpy cells are being treated as exposed Ag interface. The anomalous `Hadv` source/transport remains unproven. The replay stops about `0.023 ps` earlier than the original and finds 5 rather than 7 cells, so it localizes the immediate failure class but does not reproduce the exact event. No cooling or solidification was completed for C. See `P1P_C_FAILURE_CELL_DIAGNOSTIC.csv` and the preserved instrumented source/log in `../p1p_cell_failure_diag_20261107/`.
- The launcher recorded status 0 for C even though the solver logged a fatal STOP. The CSV reports this as `INCOMPLETE_STOP`, based on solver output, rather than treating the launcher status as success.

## Conclusion and next diagnostic

The available evidence supports **one resolved dominant liquid pool plus disconnected one-to-five-cell threshold remnants in A and B; those remnants are too small to identify as physically resolved secondary pools**. Pressure increased the model's radial Ag redistribution, but A and B retained no resolvable downward center depression after solidification. Evaporation was negligible; C failed before cooling. These data do not show that real contacts cannot form a pit or that the forces have no physical effect.

The isolated C replay confirmed the immediate status as a property-table range error and recorded cell coordinates, reconstructed interface temperatures, fractions, and advected enthalpy. Gas-dominated inversion of the logged `Hadv` explains the temperature scale, and the current permissive nonzero-fraction checks let these nearly empty interface cells reach the Ag property lookup. The original run's 7-cell event was not exactly reproduced: the replay stopped earlier with 5 cells. Before another evaporation-enabled full run, reconcile that difference and trace why enthalpy remains in cells whose VOF and material fractions are nearly zero. Validate a minimum-interface criterion and phase-consistent temperature reconstruction before changing production code. Also repair the pre-step enthalpy/property guard diagnostic so it distinguishes an inversion sentinel from an actual material temperature. No further heating run was launched as part of this report.

## Generated files

- `P1P_ABC_FULL_FINAL_COMPARISON.csv`
- `P1P_ABC_FULL_RADIAL_MASS.csv`
- `P1P_ABC_FULL_PEAK_RADIAL_MASS.csv`
- `P1P_ABC_FULL_SURFACE_SHAPE.csv`
- `P1P_ABC_FULL_PEAK_CLUSTERS.csv`
- `P1P_ABC_FULL_PEAK_COMPONENTS.csv`: per-component volumes and spatial bounds; reveals that secondary components are only isolated few-cell remnants.
- `P1P_ABC_FULL_LIQUID_CLUSTERS.csv`
- `P1P_ABC_FULL_HOLD_START_3D.png`
- `P1P_ABC_FULL_HOLD_START_ZOOM.png`
- `P1P_ABC_LIQUID_COMPONENTS.png` and `P1P_ABC_SECONDARY_POOLS.png`
- `P1P_C_FAILURE_CELL_DIAGNOSTIC.csv`: per-cell status and local state from an instrumented replay of the C failure.
- `P1P_FINAL_PIT_DECISION.md`
- `analyze_abc_full.py`, `plot_abc_full_hold_start.py`
