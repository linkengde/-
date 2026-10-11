# Ag64 step8000 candidate review

**Review status: HOLD.** This is an offline data-role audit of D commit `38ff5afb5c3cbe14ca0889212910f0ca6558e991`. No DFT, MD, MACE training, or potential changes were made.

## Identity and integrity

The candidate handoff is `Ag64_step8000_candidate_release/`. I fetched the specified D branch commit and independently recalculated the four hashes in its `SHA256SUMS.txt`; all match. The XYZ hash also matches both `MANIFEST.json` and `QA_report.json`. The XYZ has 64 Ag atoms with unique IDs 1–64, a 10.826824658074 Å cubic cell, full 3D PBC, task ID `C-D-AGL-1250K-rho097-seed01`, and `global_step=8000`. The frozen Ag potential SHA matches the project value `958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447`.

The source manifest archived in D's original warmup bundle identifies source root `PF-WMH2006-Ag64-NVT-melt-cool-1250K-six-pilot-source-campaign-20261011`, trajectory parent/source family `PF-AGL-1250K-rho097-seed01`, and PCG64 seed `125009701`. The step3000 PRELIQUID diagnostic and step8000 structure therefore belong to one trajectory family and must not be split. The step8000 release manifest records the checkpoint manifest SHA `44c5b966f52cbface56ed7a1fc135adc421288e91e76263c5ed53382e0230fa8`; its recorded state and array hashes agree with `QA_report.json`. However, the checkpoint files themselves are not included in this handoff or present in the current D branch tree, so I verified metadata consistency, not the checkpoint bytes. The new release also omits an explicit `source_family` field; the family above is derived from the original `source_manifest.json`. The shared campaign root and a different random seed do not by themselves prove independence from other D/E/F task families.

## QA review

`QA_report.json` and `MANIFEST.json` report `TWO_CONSECUTIVE_WINDOWS_ALL_VALID`. The step5000 and step6000 windows fail the temperature-mean tolerance; the step7000 and step8000 windows are the two consecutive full passes:

| Check step | Mean T (K) | T sample SD (K) | Local Q6 mean | Coordination mean | Coordination half-window drift | MSD slope (Å²/fs) | MSD fit R² | Crystalline fraction |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 7000 | 1637.233 | 155.979 | 0.37675 | 10.51042 | 0.12500 | 0.002201 | 0.9919 | 0 |
| 8000 | 1638.765 | 84.680 | 0.39634 | 9.98958 | 0.10417 | 0.001801 | 0.9887 | 0 |

Both windows meet the reported 1600 K mean tolerance (±50 K) and the finite-size sample-SD limit (329.180 K). MSD slopes are positive with high fit R². The local-Q6 trend rule passes in both windows; at step8000 the fitted slope is positive (+0.04510 per 1000 steps) with p=0.08264, which passes the report's p≥0.05 rule but is not strong evidence that local order is stationary. Coordination drift within each window passes, while the mean falls by about 0.521 between the two windows; the QA rule checks within-window half drift, so this cross-window change is a caution rather than a recorded failure. The minimum distances are 2.207 and 2.314 Å. Structural screening used local Q6 plus a coherent-bond graph as a CNA/PTM substitute; it did not report PTM.

These diagnostics describe a **1600 K hot-stage** configuration in the planned `rho097` cell, not an equilibrated 1250 K liquid. The structure looks disordered and mobile under the recorded screen, but the independent release marker below conflicts with that screen.

## Blocking inconsistency and role

`QA_report.json` and `MANIFEST.json` say the two-window screen passed, and the XYZ comment says `PRELIMINARY_LIQUID_CANDIDATE`; however, the same XYZ comment has `STRUCTURE_PASS=F`. Read as the usual T/F Boolean field, this is false and conflicts with the QA status. The original warmup source manifest also had `STRUCTURE_PASS=false`, but that earlier state does not resolve the new step8000 disagreement. Do not silently reinterpret or edit the frozen XYZ. D should provide a provenance-preserving clarification/correction that states the exact meaning of the field and ties it to the QA report.

| structure_id | Current role | Conditional DFT-label role | Source family | Reason |
|---|---|---|---|---|
| `Ag64_step8000` | **HOLD** | **DEV candidate**, after resolving the flag conflict and validating source checkpoint/family evidence | `PF-AGL-1250K-rho097-seed01` | Retaining this sole currently evidenced hot-liquid trajectory family for validation is more useful than using it for fitting. If DFT-scored as DEV, this family must never be moved to TRAIN. Liquid TRAIN labels must come from a different, independently audited source family. Never TEST. |

The fixed-family assignment also covers the step3000 PRELIQUID diagnostic and any restart or later frame from this trajectory. The step3000 frame itself remains excluded from labels because it did not pass liquid QA. No TEST data were accessed or changed.

## DFT plan if the HOLD is cleared

The proposed calculation is one fixed-geometry, fixed-cell GPAW/PBE/PAW static single point on the exact 64-atom structure, using the method and convergence rules already recorded in `FINAL_DFT_EXECUTION_PLAN.md`. Required outputs are the native total energy, force-consistent free energy, and all 64×3 forces in eV and eV/Å with input/source hashes. Archive the full stress tensor as well; use stress as a model label only if the liquid-cell k-point and stress numerical convergence checks pass. Do not use EAM energy, force, or pressure as a DFT label. The DFT label remains **DEV**, not TRAIN or TEST, if used to score candidate models.

`DFT_RELEASE=NO` until the XYZ flag is reconciled and checkpoint/source-family evidence is independently closed. This document records a data-role recommendation only; it does not authorize a calculation.
