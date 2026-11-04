# Liquid Ag viscosity calibration and NIST uncertainty audit

## Result and scope

The nine liquid-Ag measurements in Dubberstein et al. (2016) are now independently checked against the original NIST ThermoML JSON, including its per-point expanded uncertainties. The reported Arrhenius equation is consistent with all nine expanded-uncertainty intervals. These observations cover **1273–1673 K only**; this does not validate viscosity at 2000 K or 2277 K, and it does not authorize dimensional melt-pool flow predictions.

Silver melts at about 1235 K. The often-used 2277 K value is a high-temperature melt-pool condition, not the melting point.

## Sources and provenance

T. Dubberstein, M. Schürmann, H. Chaves, H.-P. Heller, and C. G. Aneziris, “A Novel Vibrating Finger Viscometer for High-Temperature Measurements in Liquid Metals and Alloys,” *International Journal of Thermophysics* 37, 100 (2016), DOI `10.1007/s10765-016-2104-7`.

The original NIST ThermoML JSON was fetched from `https://trc.nist.gov/ThermoML/10.1007/s10765-016-2104-7.json`, archived at `results/p2_nist/nist_original_10.1007_s10765-016-2104-7.json` (17,347 bytes; SHA256 `80b4ab63347b4d5abc621d52e46b003a41069fe2061eaf474225532f9294f0fc`). It records viscosity in Pa·s and temperature in K. The Ag property metadata specifies “Vibrating wire viscometry”; uncertainty metadata gives a 95% confidence level and “Propagation of evaluated standard uncertainties.” The JSON contains Au, Ag, and Sn datasets; nine values in the Ag component were selected by the NIST component identifier.

The user-uploaded `nist_ag_9points_extracted_not_original.json` is explicitly a hand-extracted subset, not the original NIST file. Its nine Ag temperature/value/expanded-uncertainty tuples match the original NIST JSON exactly. It is preserved beside the original with that provenance warning. The uploaded `download_nist_original.sh` was inspected but not executed; the original JSON was already available for this audit.

The article's Table 4 `eta_sigma` column is preserved in `ag_viscosity_experimental.csv` without assigning it a statistical definition; the article text inspected here does not define that symbol. Numerically, NIST's expanded uncertainty is exactly twice the Table 4 `eta_sigma` value at all nine points, but that numerical relation is not enough to assign the paper's undefined symbol a statistical meaning. The paper's general statement of below-30% high-temperature measurement error also cannot be directly reconciled with NIST's compiler-assessed pointwise 95% expanded uncertainty (57.36% at 1673 K) from the available metadata, so both are reported as stated.

## Fit and uncertainty check

The paper's correlation is

\[
\ln\!\left(\frac{\mu}{\mathrm{mPa\,s}}\right)=-0.4631+\frac{2089}{T},
\]

with T in kelvin. `analyze_p1b_p2.py` regenerates the raw-source extraction and residual tables. The maximum absolute fit deviation is 7.7396%; the maximum absolute residual divided by the NIST 95% expanded uncertainty is 0.3662. All nine residuals lie within one reported expanded uncertainty. This is a comparison with the publisher-reported correlation, not a new regression or a claim that NIST uncertainty fully captures every systematic error.

| T (K) | μ (mPa·s) | NIST U95 (mPa·s) | U95 / μ (%) | Fit deviation (%) |
|---:|---:|---:|---:|---:|
| 1273 | 3.520 | 0.744 | 21.14 | -7.74 |
| 1323 | 2.904 | 0.450 | 15.50 | +5.11 |
| 1373 | 2.812 | 0.646 | 22.97 | +2.48 |
| 1423 | 2.649 | 0.330 | 12.46 | +3.12 |
| 1473 | 2.540 | 0.220 | 8.66 | +2.32 |
| 1523 | 2.498 | 0.458 | 18.33 | -0.69 |
| 1573 | 2.399 | 0.352 | 14.67 | -1.01 |
| 1623 | 2.332 | 0.390 | 16.72 | -2.24 |
| 1673 | 2.214 | 1.270 | 57.36 | -0.92 |

The full raw tables are `results/p2_nist/ag_viscosity_experimental.csv` and `results/p2_nist/ag_viscosity_fit_vs_data.csv`; the source-derived summary is `results/p2_nist/nist_ag_uncertainty_summary.json`. The existing top-level experimental and fit CSVs now also include the NIST uncertainty cross-check.

## Extrapolation and use decision

The equation evaluates to approximately 1.78854 mPa·s at 2000 K and 1.57513 mPa·s at 2277 K. Both are outside the 1273–1673 K measurements and are unsupported extrapolations, not experimentally validated values. The reported primary-paper statement of <30% high-temperature measurement error and the broad pointwise NIST uncertainty, especially 57.36% at 1673 K, should accompany any use of the dataset.

Treat the curve as a literature candidate only over its measured temperature interval. Ag/Ti₃SiC₂ wall friction and slip remain uncalibrated; this viscosity audit does not resolve them.
