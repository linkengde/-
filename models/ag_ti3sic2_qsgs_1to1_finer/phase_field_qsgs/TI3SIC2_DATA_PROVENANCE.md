# Ti3SiC2 thermophysical correlation data (300–2000 K)

## What this dataset is

`ti3sic2_barsoum_correlation_curve_300_2000K.csv` is a reproducible temperature grid evaluated from published empirical correlations. It is **not** a digitization of the original experimental plot and must not be described as raw measured data.

The grid contains 25 K increments plus the important source-range boundaries (373.15, 1273.15, 1473.15, and 1873.15 K) and the configured model ceiling (2000 K).

## Correlations used

Barsoum et al. (1999), *Thermal properties of Ti3SiC2*:

- Molar heat capacity:
  `Cp_molar(T) = 164.4 - 16419/T` J/(mol K), with T in kelvin.
- Total thermal conductivity:
  `k(T) = 38.6 - 0.0045*T` W/(m K). The paper reports a fit quality of R^2 = 0.86.
- Molar mass used for conversion: 0.195708 kg/mol (project configuration).
- Mass-specific heat capacity: `Cp_mass = Cp_molar / 0.195708` J/(kg K).
- Density: 4520 kg/m^3, held constant as a project assumption. No temperature-dependent density fit is claimed.

## Source ranges and interpretation

- Barsoum reports DSC heat-capacity measurements over 100–1000 °C (373.15–1273.15 K) and discusses extrapolation to approximately 1200 °C (1473.15 K).
- Qiu et al. (2020), DOI 10.1016/j.net.2019.07.030, repeats the Barsoum correlation and states a broader 25–1600 °C Cp range. This is a secondary-source range statement, not an independent raw high-temperature dataset.
- Barsoum's thermal-conductivity fit is reported over approximately 25–1200 °C (298.15–1473.15 K).
- Values above 1473.15 K for k are extrapolations. Cp values above the primary paper's stated extrapolation to about 1200 °C are not directly validated by Barsoum; the interval to 1873.15 K is included only to record the range stated by the secondary review, and 1873.15–2000 K is an additional extrapolation.
- Togo et al. (2010), DOI 10.1103/PhysRevB.81.174301, Wang and Zhou (2004), DOI 10.1103/PhysRevB.69.144108, and Son et al. (2016), DOI 10.1016/j.commatsci.2016.08.015, are cross-check references for calculated Cp/Cv trends at lower temperatures. Wang and Zhou's listed heat-capacity values are Cv, not Cp, and should not be substituted directly for Cp.
- None of these sources supplies direct Ti3SiC2 Cp or k measurements at 2000 K.

## Column definitions

- `temperature_K`, `temperature_C`: absolute temperature in K and corresponding Celsius temperature.
- `Cp_molar_J_mol_K`: molar Cp from the Barsoum correlation.
- `Cp_mass_J_kg_K`: mass-specific Cp after conversion using the project molar mass.
- `thermal_conductivity_W_m_K`: k from the Barsoum linear correlation.
- `density_kg_m3`: fixed 4520 kg/m^3 project assumption.
- `Cp_data_status`, `k_data_status`, `density_data_status`: provenance and extrapolation flags for each row.

## Model ceiling is not a stability claim

2000 K (1726.85 °C) is the configured numerical/property-model ceiling only. It is not evidence that Ti3SiC2 remains structurally stable at this temperature. The material's high-temperature decomposition/peritectic behavior must be assessed separately from the heat-capacity and conductivity correlations.

## Reproducibility

The companion script `generate_ti3sic2_barsoum_curve.py` regenerates the CSV using only Python's standard library.
