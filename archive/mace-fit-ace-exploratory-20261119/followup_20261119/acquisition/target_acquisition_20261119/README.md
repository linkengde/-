# Acquisition follow-up — 2026-11-19

This independent archive contains three verified v24 DFT interface labels, the ACE and MTP exploratory refits, their training and scoring data, and a post-fit gate reevaluation under the updated force threshold.

## Current result

- All three registered v24 DFT labels passed SCF and archive verification.
- Current exploratory gates are energy RMSE ≤10 meV/atom, force-component RMSE ≤0.2 eV/Å, and Si hotspot center-vector error ≤0.5 eV/Å.
- The user changed the force limit from 0.1 to 0.2 eV/Å. `gate_reevaluation_force0p2.json` recomputes pass/fail from existing predictions; it does not retrain models.
- No ACE or MTP candidate passes all gates. The MTP candidates that meet energy and force limits still fail the Si hotspot center validation by a large margin.
- No candidate is qualified for heating. See `FINAL_REPORT.md` for per-model metrics.

## Contents

- `v24_verified_labels/`: DFT inputs, output labels, logs, summaries, and verification records.
- `ace_acquisition_refits/`: ACE candidates and metrics.
- `mtp_fit_candidate/`: prepared MTP training/development data and masks.
- `mtp_acquisition_refits/`: nine MTP models, logs, predictions, and metrics.
- `gate_reevaluation_force0p2.json`: updated-gate decision based on the saved predictions.
- `FINAL_REPORT.md`: current fit-gate results and next acquisition direction.
- `SHA256SUMS.txt`: SHA-256 for every archived file except this manifest.

The blind test remains unopened.
