# Dataset v10 audit: v25 ceramic supplement

Status: `passed_with_limits`. Shared training source for POD and SNAP.

- Training: 128 frames / 5491 atoms.
- Validation: 12 frames / 412 atoms.
- Held-out challenge: 2 frames / 96 atoms.
- New v25 additions: 15 training frames; 4 development frames in validation.
- Four exact v25 training duplicates match existing v9 training energies and force components exactly and were excluded.
- 6 historical interface challenge files exactly match v9 training geometries and are excluded from independent evaluation.
- 2 non-overlapping challenge structures remain (currently two Ti3SiC2 high-temperature/cooling structures).
- v25 native energies were selected; they match the existing label convention on all exact duplicate frames. Finite-width free energies are not fitted.
- v22/v23 archives are verified, but all 12 labels are already represented exactly in v9 training data.
- v24 has no verified completed DFT labels.
- There is no fresh, untouched test set. Challenge labels were used in prior v9 evaluation and are not blind for v10 selection.
- Parent v9 already supplies six small-cell Ti3SiC2 hot/cooling labels, seven 52-atom Ag/Ti3SiC2 melt/cool/solid-hold labels, two liquid-Ag training labels plus one Ag validation label, and additional small-cell active-learning interface labels.
- These inherited thermal/interface labels remain sparse small-cell samples; there is no DFT trajectory of a fully liquid Ag pool contacting a disordered ceramic skeleton.
- v25 adds 15 new 12-atom periodic Ti3SiC2 baseline/displacement/strain/random-distortion labels and four development labels; it adds no high-temperature or Ag-interface samples.
