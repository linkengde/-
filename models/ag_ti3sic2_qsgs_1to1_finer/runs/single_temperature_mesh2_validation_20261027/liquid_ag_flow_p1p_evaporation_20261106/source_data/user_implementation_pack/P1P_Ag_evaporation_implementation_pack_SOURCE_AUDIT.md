# Uploaded ZIP audit: Ag_evaporation_P1P_P1Q_source_pack(1).zip

## Integrity and contents
- ZIP CRC test: passed.
- 8 members; total archive size about 5.3 kB.
- SHA-256: `fc03754a54995705ef8405ca79bbcadb6c58ed270ed87df68d908569acd02cf3`
- The uploaded ZIP is byte-identical to the earlier generated ZIP.
- Contents are 6 CSV files and 2 Markdown files.
- No journal PDF, Supporting Information, raw experimental data file, NIST JSON, or implementation source code is present.

## Findings
1. The Wang et al. (2024) Ag vapor-pressure equation and the values at 1373 and 1773 K match the publisher-indexed article text.
2. The five Table 8 alpha values match the article. They are effective coefficients derived as actual maximum evaporation rate divided by theoretical maximum rate; they are measured/inferred only at 1373–1773 K.
3. Table 5 rates are indeed in g cm^-2 min^-1, explicitly stated in Eq. (1) of the paper. The old CSV had a needless "verify unit" warning; the replacement CSV supplies the SI conversion. Conversion is 1 g cm^-2 min^-1 = 1/6 kg m^-2 s^-1.
4. Table 7 is a logistic fit against total furnace/system pressure in an argon-controlled vacuum-distillation experiment. That pressure is not the local Ag vapor partial pressure required by the HKL back-pressure term.
5. The vaporization enthalpy table is mathematically consistent with differentiating the Wang vapor-pressure fit via Clausius-Clapeyron, but it is derived under an ideal-vapor approximation and is not a directly measured Hvap(T) table.
6. The NIST Antoine equation is correctly converted from bar to Pa. It differs from the Wang fit materially near the ends of the overlapping range; the two fits should not be averaged. At 1823 K the NIST value is about 35.5% below Wang; at 2000–2200 K it is about 5.1% above; at 2300 K it is about 4.7% below; at 2425 K it is about 19.4% below.
7. The old pack did not clearly separate molar flux from mass flux. The replacement handoff gives both correct forms and SI units. Recomputing Table 8 theoretical maximum rates from the published vapor-pressure fit and HK formula gives values about 1.2% above the rounded Table 8 values across all five temperatures; this small systematic mismatch is recorded in `theoretical_flux_formula_check.csv`. Preserve Table 8 alpha values as published rather than silently recalibrating them.
8. The previous simple recoil-pressure expression is not a universal closure. Wang et al. (2020) derives recoil from the momentum flux through a Knudsen layer; the Anisimov factor 0.54 is tied to a sonic assumption.

## Main implementation blockers
- No physical local Ag vapor partial pressure is provided for the actual arc.
- No measured alpha_ev(T) exists in the supplied Ag paper above 1773 K.
- The second-step enthalpy failure exists even with evaporation disabled; it must be fixed in the baseline before P1P can be accepted.
- At a ~0.302 nm grid and ps times, continuum evaporation may correspond to a fraction of an atom; atom-resolved loss should be validated with MD rather than amplifying alpha.

## Recommended status
- Source/data audit: passed with corrections and limitations recorded.
- P1P numerical implementation: not yet accepted.
- P1Q recoil pressure: remain disabled.
