# Static potential comparison decision

## Answer

On the single matched 135-atom relaxed interface geometry, **MACE-MP-0b3-medium, V12, V14 and V17 give substantially more reasonable static force agreement than calibration-v0 EAM+TSCBOP+Morse**. Their full force-vector RMSEs are 0.239–0.314 eV/Å versus 2.255 eV/Å for the hybrid. Effective interaction-force RMSEs are 0.423–0.455 versus 1.351 eV/Å, and their phase-subtracted interaction-energy errors are 3.68–4.35 eV rather than the old hybrid's 25.69 eV overbinding. This is a useful single-geometry improvement, not a production validation.

**V18 is not a suitable selection from this evidence.** Its full force RMSE is 2.224 eV/Å, close to the old hybrid, and its TSC-only component force RMSE is 2.158 eV/Å. The already recorded V18 48-atom pure-backbone RMSE of 2.889 eV/Å versus 0.0382 for MACE-MP remains decisive evidence of backbone regression.

No single MACE version wins across the interface metrics. V17 has the lowest total force RMSE here; V14 is closest in the composition-subtracted energy increment; the foundation model has the lowest effective interaction-force RMSE by a very small margin. All five MACE candidates still predict an Ag normal interaction force that is much stronger in magnitude than the DFT increment (DFT is −0.171 eV/Å on Ag), and the MACE energy increments are 3.28–4.35 eV less attractive than DFT. The old hybrid's force and interaction energy are far worse. A unique winner should not be chosen from one configuration.

## 135-atom matched static scores

| Model | Full force RMSE | Ag | Ti | Si | C | Eint model (eV) | ΔEint vs DFT total (eV) | Interaction-force RMSE | Ag normal interaction force |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| MACE-MP-0b3-medium | 0.243 | 0.336 | 0.165 | 0.315 | 0.239 | -18.41 | 3.82 | 0.423 | -3.255 |
| MACE-V12 | 0.314 | 0.318 | 0.307 | 0.366 | 0.291 | -18.20 | 4.03 | 0.455 | -4.041 |
| MACE-V14 | 0.281 | 0.293 | 0.278 | 0.339 | 0.243 | -18.55 | 3.68 | 0.428 | -3.354 |
| MACE-V17 | 0.239 | 0.243 | 0.230 | 0.365 | 0.154 | -17.88 | 4.35 | 0.427 | -2.249 |
| MACE-V18 | 2.224 | 0.483 | 2.520 | 1.150 | 2.693 | -18.95 | 3.28 | 0.452 | -8.914 |
| legacy_EAM_TSCBOP_Morse_calibration_v0 | 2.255 | 2.150 | 2.356 | 1.728 | 2.385 | -47.92 | -25.69 | 1.351 | -45.001 |

RMSE columns are eV/Å. The DFT effective interaction energy is `-22.229525 eV` using native total energies and `-21.378661 eV` using force-consistent free energies. MACE values are ordinary potential energies; the native-total comparison is primary and the free-energy comparison is reported separately. No absolute model energy is subtracted from an absolute DFT energy.

## 77-atom gap scan

All five frozen MACE models put 2.40 Å at the lowest relative energy of these three supplied gaps. At 2.16 Å the relative energies are about 2.51–3.30 eV above the 2.40 Å value; at 2.64 Å they are about 0.23–0.55 eV above it. At 2.40 Å, the maximum force norm is:

| Model | ΔE(2.40) eV | Fmax eV/Å | Ag normal sum eV/Å | min Ag–Si Å | min Ag–Ti Å | min Ag–C Å |
| --- | --- | --- | --- | --- | --- | --- |
| MACE-MP-0b3-medium | 0.000 | 0.652 | 3.58 | 2.400 | 4.753 | 5.814 |
| MACE-V12 | 0.000 | 0.481 | 2.61 | 2.400 | 4.753 | 5.814 |
| MACE-V14 | 0.000 | 0.531 | 2.99 | 2.400 | 4.753 | 5.814 |
| MACE-V17 | 0.000 | 0.564 | 3.39 | 2.400 | 4.753 | 5.814 |
| MACE-V18 | 0.000 | 3.472 | 1.85 | 2.400 | 4.753 | 5.814 |

No cross distance is below 1.5 Å. This curve only compares relative energies on one 0° strained parent family; the set has no DFT label, so it does not establish the true equilibrium distance, DFT accuracy, or generalization. It is not a validation set of three independent structures.

## Required status per candidate

| Model | STATIC_INTERFACE_RESULT | BACKBONE_RESULT | INTERACTION_FORCE_RESULT | DFT_REFERENCE_COMPARABILITY | SMALL_MD_READINESS | PRODUCTION_READY |
| --- | --- | --- | --- | --- | --- | --- |
| MACE-MP-0b3-medium | CONDITIONAL_IMPROVEMENT on one matched geometry; full RMSE 0.243 vs old 2.255 | Lower than old TSCBOP here (0.504 vs 2.562); single geometry only | Interaction-force RMSE 0.423 vs old 1.351; normal force too strong | Comparable matched DFT labels; family-level independence unverified | HOLD: no thermal/stability validation | NO |
| MACE V12 | CONDITIONAL_IMPROVEMENT on one geometry; full RMSE 0.314 | Lower than old backbone here (0.385 vs 2.562); single geometry only | 0.455 vs old 1.351; normal force too strong | Comparable matched DFT labels; V12 overall historical status UNDETERMINED | HOLD | NO |
| MACE V14 | CONDITIONAL_IMPROVEMENT on one geometry; full RMSE 0.281 | Lower than old backbone here (0.352 vs 2.562); single geometry only | 0.428 vs old 1.351; normal force too strong | Comparable matched DFT labels; used assessment-evaluated run-45 weights; historical screen fail | HOLD | NO |
| MACE V17 | CONDITIONAL_IMPROVEMENT; best full RMSE 0.239 on this geometry | Best of these TSC-only component RMSE 0.307 vs old 2.562 | 0.427 vs old 1.351; energy increment underbinds by 4.35 eV and normal force too strong | Comparable matched DFT labels; provisional model, prior dev gates failed, no blind test | HOLD | NO |
| MACE V18 | MIXED: full RMSE 2.224, near old 2.255 | FAIL/DEGRADED: TSC-only RMSE 2.158, retained 48-atom regression 2.889 | 0.452 vs old 1.351 but normal force is very strong | Comparable matched DFT labels; prior dev gates fail | NO: known backbone regression plus no thermal evidence | NO |
| Legacy EAM+TSCBOP+Morse v0 | NO: full RMSE 2.255, interaction E overbinds by 25.69 eV | NO: TSC-only RMSE 2.562 | NO: interaction force RMSE 1.351; Ag normal force −45.00 vs DFT −0.171 | Comparable matched DFT labels; source/run0 reproducible | HOLD for new use; static validation fails | NO |

## Decision and limitations

- MACE versions V12/V14/V17 and the MACE-MP foundation weights are better static candidates than the old hybrid for this single matched DFT geometry. V17 is the best full-force row here, but prior development gates failed and it has no untouched test; do not select it for production on this basis.
- V14 uses the exact run-45 model SHA that its existing validation assessment identifies as the evaluated model. The separately recorded selected epoch-71 checkpoint has another SHA; this comparison does not silently conflate them.
- V18 is rejected as a backbone candidate despite its lower effective cross-force error than the old hybrid; a cross-term improvement cannot hide backbone degradation.
- The 77-atom structures are only a frozen numerical gap-response check. No DFT data exist for them.
- `SMALL_MD_READINESS` is HOLD for every model because this request only authorized static inference; none received a thermal stability or short-MD check in this task. `PRODUCTION_READY=NO` for all models.
- The 135-atom geometry had no exact frame match in the scanned V12/V14/V17/V18 train/valid/available test partitions. Parent-family independence remains unverified, so this is not described as a blind benchmark.
