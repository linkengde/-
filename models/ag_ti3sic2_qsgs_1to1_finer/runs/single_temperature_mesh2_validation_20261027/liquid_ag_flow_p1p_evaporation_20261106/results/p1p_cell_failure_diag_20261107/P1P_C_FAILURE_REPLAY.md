# C-case evaporation guard: instrumented replay

This is a short diagnostic replay from the verified C-case restart at `9.00332084 ps`. It does not alter the production solver or original run outputs. The only source change is additional per-cell logging in `p1p_m1m_evaporation_debug.c`; the production source remains byte-identical to the campaign source SHA below.

## Result

The replay stopped before the thermal update at `9.28035319592 ps`. Five interface cells returned status `3` (`P1P_ERR_RANGE`) from the evaporation property lookup. Their reconstructed interface temperatures ranged from `4925.17` to `5319.75 K`, outside the table range `1234–2436 K`. The cell fractions were extremely small: sample fraction `5.96×10⁻¹¹` to `1.17×10⁻⁹`, Ag fraction `1.11×10⁻¹²` to `5.05×10⁻¹¹`, and Ti₃SiC₂ fraction `7.997×10⁻¹¹` to `1.119×10⁻⁹`; pore fraction was zero.

This identifies the immediate stop as a table-range error produced by anomalously high reconstructed interface temperature. It does not identify the upstream cause. The original C run stopped at `9.30330989 ps` with 7 rejected cells, whereas this instrumented replay stopped about `0.023 ps` earlier with 5. Treat this as diagnostic localization, not exact reproduction. No cooling or solidification was run for C.

The solver's mixture inversion uses `air_fraction = max(0, 1 - Ag_fraction - Ti3SiC2_fraction)` and air volumetric heat capacity `rho_air Cp_air = 1.2 × 1005 = 1206 J m⁻³ K⁻¹`. In these cells the solid fractions are below `1.2×10⁻⁹`, so the gas term dominates. The simple gas-dominated estimate `300 K + Hadv/(1206 J m⁻³ K⁻¹)` is `4293–4911 K`; the actual interface reconstruction is `4925–5320 K` after adding its local gradient extrapolation. Thus the rejected temperatures are consistent with gas-dominated cell enthalpy inversion followed by classifying a vanishingly small VOF/Ag trace as an exposed Ag interface. This is a strong mechanism hypothesis from the logged state, not yet a proven root cause: the source and transport of the anomalous `Hadv` still need tracing.

The current evaporation loop accepts any interface cell with `f > 1×10⁻¹²`, then accepts any Ag coverage above `1×10⁻¹²`. The five replay cells pass those checks: `f` is `5.96×10⁻¹¹–1.17×10⁻⁹`, Ag coverage is `0.0137–0.0432`, and PLIC area is only `4.45×10⁻⁸–3.86×10⁻⁷` of a nominal cell face. They therefore pass the current “nonzero interface” test despite containing vanishingly small total material volume. This explains how gas-dominated cells reach the Ag evaporation property lookup. No threshold or production code was changed; a defensible minimum-interface criterion and a phase-consistent temperature reconstruction need separate validation.

## Preserved artifacts and hashes

- Restart field SHA256: `685f96295a5618a8bb76b40ac65c3e8f4e547bcc147d647a6ecc53ed723fb07a`
- Production solver source SHA256: `098a64e3d7123114929d1cd09f27ea7ea1d96d2eac6b9c4091f7b7afbcd22bdd`
- Instrumented replay source SHA256: `9be6910661b00d5d573d2125de9f02ac0ededf4e686c6ee10dfddd163bc99b5c`
- Instrumented replay executable SHA256: `a4fea1bfd7dcc20a3752b5468c6643f6a5a3b554d4024990e8437b0c89ce9877`
- `debug.stderr`: complete replay log, SHA256 `cc22573e7efdcfd1859d99731d12562a45dea408e2c7e82dde58396bffad3b18`
- `P1P_C_FAILURE_CELL_DIAGNOSTIC.csv` in the sibling combined-diagnosis directory: parsed cell values and gas-dominated inversion estimates.

The original 7-cell C stop remains preserved in the original run log. The replay is not evidence that C completed and does not change the permanent-pit conclusion: C has no cooled final surface.
