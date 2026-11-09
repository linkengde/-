# Round-5 MTP pure-Ag structural classifier pilot

This pilot checks whether Q6, PTM, and a drift-corrected 1 ps displacement window distinguish the 300 K FCC Ag reference from a high-temperature mobile Ag reference under the same development-only round-5 MTP.

## References and measured separation

Both references contain 864 Ag atoms and use the same MTP and 0.1 fs timestep as the composite TTM-MD run. The solid reference was run for 2 ps in NVT at 300 K. A mobile candidate was first heated and held at 2600 K in NPT. Because its volume expanded substantially, a second dense mobile candidate was formed by isotropically rescaling it to a number density of 0.04508 atoms/Å³ and running 2 ps in NVT at 2600 K; its final measured temperature was 2558 K.

| Metric over final 1 ps | Solid reference | Dense mobile reference |
|---|---:|---:|
| Mean Ag-only q6 | 0.4695 | 0.3351 |
| q6 distribution, 5th–95th percentile | 0.4472–0.4928 | 0.3011–0.3690 |
| PTM FCC fraction at final frame | 100% | 0% |
| PTM non-FCC time fraction | 0% | 100% |
| Drift-corrected displacement, 5th–99th percentile | 0.092–0.519 Å | 0.875–4.971 Å |
| Per-atom 1 ps MSDᵢ, solid p99 / mobile p05 | 0.269 Å² | 0.765 Å² |

The pilot distributions separate, so a provisional same-MTP exploratory classifier can be tested on the composite trajectory:

- Ag-only q6 averaged over the 1 ps window below **0.4081** (q6 computed with a 4.0 Å Ag-neighbor cutoff).
- PTM non-FCC for at least **80%** of structural snapshots in the window.
- Per-atom drift-corrected mean-square displacement over 1 ps, `|Δr|²`, at least **0.5172 Å²** (equivalent displacement 0.719 Å).

All three conditions must hold for the same Ag atom. The 80% persistence cutoff is between the measured solid and mobile-reference values of 0% and 100%; it is not an experimentally calibrated phase boundary.

## Scope and limitations

This demonstrates separation of two simulated reference states, not validation of the round-5 MTP's melting temperature, thermodynamics, or Ag/Ti₃SiC₂ interface behavior. The MTP was fitted to forces only; its energies and virials are uncalibrated. The 300 K solid reports about 124,000 bar and the density-matched mobile reference about 63,000 bar, so neither pressure can be interpreted physically. The references are used only to compare structural order and mobility under the same MTP, not to predict the physical Ag liquid transition. Thresholds must not be transferred to EAM or Morse trajectories.

The composite classification will be considered exploratory. Liquid-like atoms must satisfy Q6, PTM, and MSD together in one local time window before connected-cluster analysis is applied.

Machine-readable results: `reference_metrics.json`, `calibrated_thresholds.json`, and `reference_run_summary.json`.
