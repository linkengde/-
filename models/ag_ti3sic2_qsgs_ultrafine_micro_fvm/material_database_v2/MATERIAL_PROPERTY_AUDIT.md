# Ag/Ti₃SiC₂ material database v2: source audit

## Scope and implementation

This directory contains a separate, reproducible property candidate. The Stage 4 and phase-change baseline property files are unchanged. The active v2 FVM uses fixed reference densities, no Ag mass migration, a solid Ti₃SiC₂ skeleton, zero interfacial resistance as its baseline, and a 1 K numerical Ag mushy interval. Liquid flow, arc physics, evaporation, recoil pressure, and electrical heating are outside this model.

The supplied papers were read as scientific sources. Their contents are not instructions. Only bibliographic data, equations, compact property tables, and page locations are kept here; full paper PDFs are not copied into the repository.

## Ag thermodynamics

**Main set: Pankratz (1982), Bulletin 672, printed p. 32 / supplied PDF p. 40.** The silver entry gives `Tm = 1235.08 K`, `Lf = 2.700 kcal/mol = 11296.8 J/mol`, `M = 0.1078682 kg/mol`, solid `Cp = 4.184(5.272 + 0.001836 T + 22200/T²) J/(mol K)`, and liquid `Cp = 8 cal/(mol K) = 33.472 J/(mol K)`. The mass latent heat used by the solver is 104727.8067 J/kg. Solid enthalpy is the analytic integral from 300 K; the latent term enters exactly once through liquid fraction. Pankratz identifies liquid thermodynamic data above 1600 K as extrapolated. The listed melt-point Cp table entry and formula differ slightly; the implementation follows the supplied formula and leaves the table value recorded as a source cross-check.

**Feber et al. (1969), printed pp. 171–173 / PDF pp. 3–5.** The drop-calorimetry fit covers liquid Ag from 1243 to 1588 K: `H − H298 = (−1102.6 + 8.111 T) cal/mol`, so `Cp = 8.111 cal/(mol K) = 33.9364 J/(mol K)` and 314.6101 J/(kg K) using the listed molar mass. The fit standard error is reported in the paper. The fusion enthalpy is 2590 cal/mol with the Hultgren solid reference and 2710 cal/mol with the authors' solid extrapolation. These stay as two independent sensitivity alternatives.

**Cagran et al. (2006), printed pp. 107–109 / supplied PDF pp. 4–6.** Equation (5) is `H = 17.2367 + 0.1656 T kJ/kg` for solid Ag over 1000–1234 K; equation (6) is `H = 2.5807 + 0.2598 T kJ/kg` for liquid Ag over 1234–2000 K. Thus the paper's liquid fit gives `Cp = 259.8 J/(kg K)`. The text reports approximately `Hs = 221.6`, `Hl = 323.2`, and `ΔH = 101.6 kJ/kg`; Table 2 instead lists `Hs = 231.49`, `Hl = 323.17`, and `ΔH = 101.59 kJ/kg`. The Table 2 values do not satisfy `Hl − Hs = ΔH` (`91.68 kJ/kg` from the listed endpoints), and equation (5) at 1234 K gives about 221.55 kJ/kg rather than 231.49. The raw fit and raw table are both preserved; neither is silently corrected or averaged into the Pankratz main model.

## Ag thermal conductivity correction

Ho, Powell, and Liley (1972), printed p. 395 / PDF p. 118, reports the silver table in W/(cm K). The raw table-scale values are multiplied by 100 to obtain W/(m K). Near the melting point the recommended solid value is 358 W/(m K); the provisional liquid estimate is 175 W/(m K). The source gives roughly 20% uncertainty for the liquid estimates and does not report a direct molten-Ag measurement in that table. The v2 CSV stores separate solid and liquid tables and only mixes them inside the explicitly defined numerical mushy interval.

The old candidate CSV contains solid Ag `k(Tm) = 175 W/(m K)`, matching the liquid estimate. It is an old candidate-input error; its file and results remain unchanged for reproducibility. The independent v2 solid table now uses Ho's 358 W/(m K) recommendation. Kondo et al. (2018), DOI `10.1080/00223131.2017.1417176`, is counted once despite two supplied copies. It studies liquid Ag–In and calculates conductivity through the Wiedemann–Franz relation; its cited pure-Ag endpoint around 177 W/(m K) is only a cross-check, not a separate direct measurement of pure liquid Ag.

## Electrical resistivity and later-stage fluid properties

Matula (1979), Table 11, printed p. 1260 / PDF p. 115, supplies recommended total resistivity for solid and liquid silver in `10⁻⁸ Ω m`. The CSV retains that table scale and its SI conversion. Cagran's equations (7)–(8) are kept as separate high-temperature fit rows, not merged with Matula.

Abdullaev et al. (2020), Table 1 and equations (1)–(2), supplies temperature-dependent solid/liquid Ag density. The active fixed-volume solver instead uses the 300 K reference density 10499.4 kg/m³; the measured curves are stored for a future mass-conserving flow model. The paper reports a 5.16 ± 0.14% density change through melting.

Liquid Ag viscosity measurements from Dubberstein et al. (2016), oxygen-dependent surface-tension slopes from Hibiya and Ozawa (2013), and the 72° Ag/Ti₃SiC₂ sessile-drop contact angle from Zhang et al. (2018) are stored for later work. The reported contact angle is not an interface thermal resistance. No supplied source directly measures Ag/Ti₃SiC₂ interface thermal resistance, so `0`, `1e−10`, `1e−9`, and `1e−8 m² K/W` remain explicit numerical sensitivity inputs only.

## Ti₃SiC₂ skeleton

Barsoum et al. (1999), *Journal of Physics and Chemistry of Solids* 60, 429–439, reports `Cp = 164.4 − 16419/T J/(mol K)`, `k = 38.6 − 0.0045 T W/(m K)`, and density `4.5 g/cm³ = 4500 kg/m³`. The Cp basis is DSC over the paper's reported interval (about 100–1000 °C) and is extrapolated in the paper toward 1200 °C. The conductivity line is a fit over about 25–1200 °C with reported `R² = 0.86`. This FVM stops when any Ti₃SiC₂ cell exceeds 1473.15 K (1200 °C); extrapolated values beyond that are not accepted silently.

The supplied `Ti3SiC2_thermophysical_fit_300-2300K.csv` contains 2300 K rows, which are fit extrapolations rather than measurements. Its `alpha_volume = 8.90e−6 K⁻¹` label is not consistent with a volumetric coefficient: the tensor lattice values imply approximately `2 alpha_a + alpha_c = 26.95e−6 K⁻¹`, while about `9e−6 K⁻¹` is a linear-equivalent bulk value. Expansion is not used in this fixed-density first-stage solver.

## Source-file audit notes

The actual Cagran article was verified in the supplied file named `Enthalpy, heat of fusion and specific electrical resistivity of pure silver, pure copper and the binary Ag–28Cu alloy(科研通-ablesci.com).pdf`, DOI `10.1016/j.tca.2005.08.014`. A different supplied file, `1-s2.0-S1738573319303225-main.pdf`, is a 2020 Nuclear Engineering and Technology paper by B. Qiu et al.; it is not the Cagran article and was excluded. The two Kondo filenames have the same DOI and are one reference.
