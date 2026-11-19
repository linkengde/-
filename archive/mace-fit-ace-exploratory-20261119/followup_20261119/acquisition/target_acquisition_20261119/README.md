# Acquisition follow-up — 2026-11-19

This folder records the prepared data and active label-acquisition plan for the exploratory Ag–Ti₃SiC₂ ACE fit. No new fit result is claimed yet.

## Existing ACE baseline

- Periodic development: energy RMSE 4.530 meV/atom and force-component RMSE 0.05054 eV/Å; both meet the agreed exploratory limits (10 meV/atom and 0.1 eV/Å).
- Registered v11 Si-center development force error: 8.417 eV/Å; this fails the 0.5 eV/Å center-vector limit. The v11 development frame is not a blind test.

## Prepared data

- The 1,021,793-atom QSGS/MACE source was screened at 81 candidate centers. The four crop geometries here have 26 Å cells and 7 Å radius, with a 12 Å gap to periodic images. Crops contain 68–104 atoms; coordinate and source-ID audits passed. They are geometry candidates only; they do not yet have DFT labels.
- Crop center roles are three training centers (Ag 819467, Ti 35127, Si 783134) and one Si development holdout (558025). The holdout is exploratory, not blind.
- The independently archived v11 labels are retained with their source audit and split: indices 00/02/03 train; index 01 development.
- The ACE refit will compare Si-center force weights 1, 10, and 100 while assigning zero weight to crop energy and crop-edge forces. It will also score aggregate periodic train/development errors and direct registered center-force vectors.

## Active v24 acquisition

As of 2026-11-19T01:34:45.883289+00:00, sequence status is `running`. No selected v24 archive has completed. The first running label is `AgTi_JN1967_registryB_surfaceC_dz_p0p02_sigma0p10_k6x6x1_v24`; selected follow-ups are AgTi +0.02 Å and AgC −0.02/+0.02 Å interface response structures. Inputs, runtime scripts, and the sequence status snapshot are included under `v24_c_response_run/`. The live workspace supervisor will verify each output before the automatic fit stage starts.

## Provenance and next action

An automatic archive step waits for the DFT and fit status files, validates labels, writes the report and SHA-256 manifest, and pushes only this independent branch. The fitting candidate data is built only after all three registered v24 outputs pass SCF, input-hash, result, and archive checks. After the fit, the model weights that meet every exploratory gate will be recorded together with per-configuration metrics. If no candidate meets the Si-center gate, the four direct QSGS crop geometries are the next DFT labeling set; no DFT values are inferred from the MACE forces.

The manifest `SHA256SUMS.txt` covers every file in this folder except itself.
