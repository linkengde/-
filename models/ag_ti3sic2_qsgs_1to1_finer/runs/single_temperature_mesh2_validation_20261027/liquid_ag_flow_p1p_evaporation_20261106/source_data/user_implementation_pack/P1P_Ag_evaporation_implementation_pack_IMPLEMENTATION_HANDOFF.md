# P1P Ag evaporation implementation handoff (audited)

## Decision
The uploaded ZIP is structurally valid, but it contains only 8 small CSV/MD files: no publisher PDFs, no Supporting Information, no source code, and no raw experimental files. The replacement pack adds unit-normalized tables, source scope, equations, and implementation constraints. It does not claim that external publisher PDFs were downloaded.

## 1. Verified primary Ag source
Wang et al. (2024), “Vacuum evaporation and condensation thermodynamics and evaporation kinetics of pure silver”, DOI 10.1016/j.nxmate.2024.100189, is open access. The publisher-indexed full text identifies:
- Table 1: Ag vapor-pressure correlation.
- Table 5: measured evaporation rates.
- Table 7: logistic pressure-fit parameters.
- Table 8: effective evaporation coefficient alpha.
- Section 3.4: overall vacuum kinetics and Hertz-Knudsen-Langmuir interfacial flux.

The article's experimental context is 99.995 wt% Ag in a high-temperature vacuum distillation furnace with argon pressure control, 1373–1773 K and system pressure 1–100 Pa. These values are useful pure-Ag references but are not a direct electrical-arc calibration.

## 2. Vapor pressure
Primary Wang et al. (2024) liquid-Ag fit, valid 1234–2436 K:
`log10(p_sat[Pa]) = -14400/T - 0.85*log10(T) + 13.825`

NIST Silver Antoine cross-check, valid 1823–2425 K:
`log10(P[bar]) = 1.95303 - 2505.533/(T - 1194.947)`

The two correlations differ by about 35.5% at 1823 K, about 5.1% at 2000–2200 K, about 4.7% at 2300 K, and about 19.4% at 2425 K (NIST relative to Wang). Do not average the fits. Use Wang as the primary correlation for the current P1P 1234–2436 K range and keep NIST as a comparison/uncertainty flag.

## 3. Hertz-Knudsen-Langmuir flux and units
Wang et al. Eq. (20) is written as a molar flux:
`N_Ag = alpha_ev*(p_eq - p_Ag,v)/sqrt(2*pi*M_Ag*R_u*T_s)`
where `N_Ag` is mol m^-2 s^-1, `M_Ag` is kg mol^-1, and `R_u` is J mol^-1 K^-1.

For a solver that updates mass, convert to mass flux:
`j_Ag = M_Ag*N_Ag`
`j_Ag = alpha_ev*(p_eq - p_Ag,v)*sqrt(M_Ag/(2*pi*R_u*T_s))`
or equivalently:
`j_Ag = alpha_ev*(p_eq - p_Ag,v)/sqrt(2*pi*R_specific*T_s)`,
where `R_specific=R_u/M_Ag` in J kg^-1 K^-1.

Use SI internally:
- pressure: Pa
- temperature: K
- mass flux: kg m^-2 s^-1
- molar flux: mol m^-2 s^-1
- energy flux: W m^-2

If the code uses `j_Ag` as mass flux, do NOT use the molar-flux denominator without multiplying by molar mass. This is a common factor/unit error.

Use the signed pressure difference. For a one-way evaporation-only model, clipping negative flux to zero is an explicit modeling assumption. If local vapor can condense, retain signed net flux and model condensation consistently.

## 4. Alpha_ev: what the values mean
Wang et al. Table 8:
- 1373 K: 0.2118
- 1473 K: 0.0834
- 1573 K: 0.1428
- 1673 K: 0.1329
- 1773 K: 0.2715

The paper defines these as the ratio of experimentally extrapolated maximum evaporation rate to the theoretical maximum rate (the theoretical expression uses alpha=1). It then uses alpha as the volatile coefficient in its HKL equation. These are effective coefficients inferred under that furnace experiment, not universal constants or direct arc measurements.

For 1773–2300 K there are no alpha measurements in this paper. Do not fit or extrapolate a validated alpha(T). If exploratory sensitivity is required, use alpha = 0.0834 / 0.1428 / 0.2715 as low/central/high brackets, explicitly labelled as sensitivity scenarios borrowed from the measured range, not high-temperature data. Alpha=1 is only the theoretical maximum limit.

## 5. Table 5 and Table 7 pressure warning
Table 5 evaporation-rate unit is explicitly `g cm^-2 min^-1` (Eq. 1 in the article). Conversion:
`1 g cm^-2 min^-1 = 1/6 kg m^-2 s^-1`.

Table 7 logistic fit:
`omega = A2 + (A1-A2)/(1+(p_system/p0)^m)`.
Here `p_system` is the total furnace/system pressure in the argon-controlled experiment. It is NOT the local Ag vapor partial pressure `p_Ag,v` used in the HKL interface equation. Do not insert the furnace total pressure into the HKL back-pressure slot. The logistic fit is apparatus-specific empirical evidence, not a substitute for local vapor transport.

## 6. Vaporization enthalpy
The pack's `Hvap(T)` curve is DERIVED, not directly measured:
`DeltaHvap,m = R_u*T^2*d(ln p_sat)/dT`
For the Wang fit:
`DeltaHvap,m = R_u*(14400*ln(10) - 0.85*T)` J mol^-1.
Mass-specific:
`DeltaHvap = DeltaHvap,m/M_Ag` J kg^-1.

This uses Clausius-Clapeyron with an ideal-vapor approximation. Mark it as a derived estimate and compare against independent boiling-point reference data; do not describe it as direct high-temperature calorimetry.

Energy coupling:
- If the solver evolves only condensed-phase enthalpy, apply `q_ev = j_Ag*DeltaHvap` as an outward surface heat flux.
- If an explicit vapor phase and conservative total-energy transport already carry the vapor enthalpy, do not also subtract `q_ev` separately. Avoid double-counting.
- The evaporation sink must be coupled to the exact same mass flux used in the Ag mass update.

## 7. Interface localization and nanoscale limits
- Evaporation applies only to Ag–gas interface exposed to the external gas and connected to it.
- Do not evaporate Ag at internal Ag/Ti3SiC2 interfaces or inside closed Ag regions.
- Reconstruct `T_s` at the geometric VOF interface (PLIC/interface-midplane interpolation); do not use a gas-cell temperature or a smeared cell average. Vapor pressure is exponentially temperature-sensitive.
- Integrate the interfacial mass flux over actual interface area. If converted to a volumetric source, normalize the interface delta function so its volume integral reproduces the surface integral.
- Update VOF and the Ag composition marker from the same conservative interfacial mass source. Prevent marker/VOF decoupling and near-gas spurious Ag.
- Limit each step's removed mass to available Ag mass in interface cells and subcycle when necessary.

The current model's ~0.302 nm grid is close to atomic spacing. A continuum fractional mass source can represent less than one atom over the entire simulated interval. At 2300 K, with p_Ag,v=0 and alpha brackets borrowed from the 1373–1773 K measurements, a 100 nm^2 surface over 100 ps corresponds to only about 0.22–0.73 expected Ag atoms. Over 5 ps it is about 0.011–0.036 atoms. This is only a scale estimate (not a validated prediction), but it means a visible crater cannot be expected from physical evaporation over a few ps in this nanoscale box. For atom-resolved evaporation, use a validated MD potential and discrete atomic escape; do not inflate alpha artificially.

## 8. Gas-regime / back-pressure decision
The Ag paper gives a gas molecular diameter of 2.88 Å and a mean-free-path relation. The actual arc gas mixture, local Ag vapor partial pressure, and characteristic plume length are not specified in the uploaded pack. Therefore P1P must expose `p_Ag,v` as an explicit input:
- `p_Ag,v = 0` is a free-evaporation upper-bound case only.
- Other `p_Ag,v` values require physical justification or a sensitivity sweep.
- Do not substitute total ambient/argon pressure for `p_Ag,v`.

The gas may be rarefied at a 3–10 nm characteristic length. Estimate Kn using the actual Ag-background-gas collision diameter, temperature, pressure and length before adding a continuum vapor-flow solver.

## 9. Required P1P acceptance before P1Q
1. First reproduce the existing second-step enthalpy failure in the no-evaporation baseline. Fix/understand that independent solver failure before attributing anything to evaporation.
2. Single-interface analytical HKL test against the equations above.
3. Zero pressure-difference test gives zero net flux.
4. Verify unit conversion using Table 8 theoretical maximum rates.
5. Multi-step closed-box mass ledger: condensed Ag + vapor/escaped Ag = initial Ag.
6. Energy ledger with no double-counting.
7. Time-step and interface-area convergence.
8. No Ag marker in near-pure gas cells; no evaporation on internal Ag/Ti3SiC2 interface.
9. Real-Ag short run only after the above tests pass, with p_Ag,v and alpha scenario recorded.

P1Q recoil pressure remains disabled until P1P multi-step tests pass and the vapor-gas regime is specified.

## 10. Source PDFs / Supporting Information status
The uploaded ZIP did NOT include any publisher PDF or SI. Official source links are provided in `source_manifest.csv`. The 2024 Ag paper's key numerical tables and equations are in the main article; no separate SI was identified in the accessible record. Much et al. (2024) and Wang et al. (2020) are methodological references, not Ag-specific parameter sources.
