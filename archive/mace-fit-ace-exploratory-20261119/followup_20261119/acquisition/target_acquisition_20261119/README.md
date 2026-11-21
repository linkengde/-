# Acquisition follow-up — 2026-11-19

This independent archive contains three verified v24 DFT interface labels, the ACE and MTP exploratory refits, and their full scoring inputs and outputs.

## Current result

- All three registered v24 DFT labels passed SCF and archive verification.
- The three ACE refits and nine MTP candidates were evaluated against the agreed exploratory gates: energy RMSE ≤10 meV/atom, force-component RMSE ≤0.1 eV/Å, and held-out Si hotspot center-vector error ≤0.5 eV/Å.
- No candidate passed all gates. The ACE candidates fit periodic energy and force metrics but miss the held-out Si center-force vector. The MTP sweep also misses force and hotspot gates.
- No candidate is approved for heating. See `FINAL_REPORT.md` for per-model metrics.

## Contents

- `v24_verified_labels/`: DFT inputs, output labels, logs, summaries, and verification records.
- `ace_acquisition_refits/`: ACE candidates and metrics.
- `mtp_fit_candidate/`: prepared MTP training/development data, masks, and manifest.
- `mtp_acquisition_refits/`: nine MTP models, logs, predictions, and metrics.
- `FINAL_REPORT.md`: fit-gate results and next data-acquisition direction.
- `SHA256SUMS.txt`: SHA-256 for each archived file except this manifest.

The blind test remains unopened. The v11 index 01 hotspot is a development label, not a blind test.
