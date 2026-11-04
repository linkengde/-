# Liquid Ag viscosity calibration

## Result and scope

The primary-paper data and its reported Arrhenius fit are transcribed and checked over the paper’s nine-point experimental interval, 1273–1673 K. This provides an auditable candidate viscosity curve **within that interval only**. It does not establish a validated viscosity at 2000 K or at 2277 K. Silver melts at about 1235 K; 2277 K is a high-temperature melt-pool condition considered in earlier project work, not the Ag melting point. No extrapolated viscosity authorizes a dimensional melt-pool velocity calculation.

## Source and provenance

T. Dubberstein, M. Schürmann, H. Chaves, H.-P. Heller, and C. G. Aneziris, “A Novel Vibrating Finger Viscometer for High-Temperature Measurements in Liquid Metals and Alloys,” *International Journal of Thermophysics* 37, 100 (2016), DOI [10.1007/s10765-016-2104-7](https://doi.org/10.1007/s10765-016-2104-7).

The user-provided article PDF was inspected locally (SHA256 `34bb1a81fdabfd81c638dc33be45f4c8ea8c956fae22f80c818af4457c1c4963`). Table 4 reports nine Ag points and the paper gives

\[
\ln\!\left(\frac{\mu}{\mathrm{mPa\,s}}\right)=-0.4631+\frac{2089}{T},
\]

with absolute temperature T in kelvin. The article states that the relative measurement error for the high-temperature metal measurements is below 30%. Its Table 4 includes an additional `ησ` column; the inspected text does not define that symbol, so the values are retained in the CSV without interpreting them as standard deviations or confidence limits.

The experimental points were transcribed from the primary paper’s Table 4. The user reports independently inspecting the NIST ThermoML JSON record (`https://trc.nist.gov/ThermoML/10.1007/s10765-016-2104-7.json`) in the current chat environment: its nine Ag values match the paper-table transcription, and the JSON contains per-point expanded-uncertainty fields. This is recorded as external cross-check evidence. The JSON and its numeric uncertainty fields are not yet archived in this repository: direct requests from this runner returned HTTP 403. Per the user's instruction, the user will source the data separately. The uncertainty values, units, and coverage definition can be calculated and checked once that local file is available.

## Fit check against the nine tabulated points

Residual is defined as 100×(fit−reported)/reported. The fit has MAPE 2.8479%, maximum absolute residual 7.7396%, mean signed residual +0.0469%, and log-space RMSE 0.03620. The full pointwise comparison is in `ag_viscosity_experimental.csv`; `ag_viscosity_fit.csv` contains the fitted curve only within 1273–1673 K.

| T (K) | Reported μ (mPa·s) | Fit μ (mPa·s) | Residual (%) |
|---:|---:|---:|---:|
| 1273 | 3.520 | 3.247565 | -7.7396 |
| 1323 | 2.904 | 3.052275 | +5.1059 |
| 1373 | 2.812 | 2.881716 | +2.4792 |
| 1423 | 2.649 | 2.731703 | +3.1221 |
| 1473 | 2.540 | 2.598915 | +2.3195 |
| 1523 | 2.498 | 2.480685 | -0.6931 |
| 1573 | 2.399 | 2.374853 | -1.0065 |
| 1623 | 2.332 | 2.279652 | -2.2448 |
| 1673 | 2.214 | 2.193624 | -0.9203 |

The fit is a compact correlation through this data set, not a claim that its percent residual is the full experimental uncertainty. The paper’s <30% high-temperature accuracy statement is broader and materially larger than the fit residuals. The paper describes its values as close to Kehr (2009) and plots comparisons with Assael et al. (2012), Gebhardt, Iida, Martin-Garin, and Schenck; no pointwise independent data from those sources were digitized in this audit, so no quantitative cross-source validation is claimed.

For clarity, formal evaluations of this equation are about 1.78854 mPa·s at 2000 K and 1.57513 mPa·s at 2277 K. Both lie outside the reported measurement interval and are **unsupported extrapolations**, not measurements or validated material-property values. They are excluded from the fit CSV and must not be used as calibrated inputs without independent evidence.

## Use decision

Use the correlation as a literature candidate only inside 1273–1673 K, preserving its broad stated measurement error and atmosphere/sample limitations. The nine-point NIST value match is externally confirmed as reported by the user; archiving the JSON and extracting its per-point uncertainty values remain open. Do not extrapolate it to the full melt-pool temperature range or use it to produce a real Ag70 flow prediction.
