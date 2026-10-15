# Ag–Ti₃SiC₂ cross-potential decision

## Required status fields

| Status | Decision | Basis |
|---|---|---|
| `TSCBOP_STATUS` | **RETAIN_UNCHANGED** | Backbone force errors remain documented limitations. They are not assigned to Ag cross terms or used as fit residuals. SHA-256 `729273a6488d75756ed7927432c77c7347a1b115cecbaefae4981bafcec80fd7`. |
| `AG_EAM_STATUS` | **RETAIN_UNCHANGED** | Ag intraphase reference remains fixed. SHA-256 `cb3d605cd27ec6efb50c4dd7462c2ab3d967e63775dc4400de9134315a44d807`. |
| `CROSS_POTENTIAL_DATA_READY` | **NO — partial evidence only** | One matched phase-resolved geometry and six reused-family local references; no pair-resolved independent distance/registry coverage or untouched validation family. |
| `MORSE_FIT_FEASIBLE` | **NO — not demonstrated** | Nine independent Morse parameters are not established identifiable; current Ag–C central-force basis also fails the stored geometry lower-bound gate. |
| `MORSE_TRIAL_FIT_COMPLETED` | **NO** | No optimizer or parameter search was run. Calibration-v0 remains unchanged. |
| `ADDITIONAL_DFT_REQUIRED` | **YES** | See `DFT_REQUEST_TO_C.md`; C must freeze structures, method, roles, and the Ag–C model/gate decision first. |
| `CROSS_POTENTIAL_PRODUCTION_READY` | **NO** | No candidate fit, untouched validation pass, or small-model stability evidence exists. Historical Stage67 stable-pool results do not validate this calibration independently. |

## Results and next action

The matched 135-atom phase subtraction gives a useful single-point cross interaction energy (−21.378661 eV in the force-consistent free-energy convention) and effective cross-force RMSE of 1.3513 eV/Å overall. Existing V19 candidate02 local tests show severe marked-pair errors and are not independent. The Ag–C radial-force lower bound warns that the present Morse functional form cannot meet the stored stringent force gate on two diagnostics, without implying all other models or configurations fail.

Keep the current calibration-v0 intact for historical reproducibility. E requests C's review of the proposed minimum DFT package and Ag–C acceptance decision. No fitted parameter file should be created until matched structures exist, the cross-force Jacobian is shown full-rank/well-conditioned, and a distinct family has been frozen outside fit selection. If those gates pass later, create a new versioned candidate beside calibration-v0 and validate only on small models before any heating study.

## Scope closure

This task did not modify Stage67 trajectories, heat source, potential files, DFT labels, structure files or role registry. It did not run DFT, LAMMPS, minimization, MD, TTM or fitting. The only new artifacts are the reports in this directory.
