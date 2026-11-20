# DFT acquisition and ACE refit report

Updated: 2026-11-20T23:16:07.168889+00:00
DFT sequence status: `failed_stop_no_retry`
ACE pipeline status: `blocked_by_dft_stop`

## Registered exploratory gates

- Energy RMSE: ≤ 10 meV/atom on full-label train and periodic development.
- Force-component RMSE: ≤ 0.1 eV/Å on full-label train and periodic development.
- Registered Si hotspot center-force vector error: ≤ 0.5 eV/Å for training centers and the v11 development center.
- The v11 development frame is not a blind test.

## Verified DFT labels

- `AgTi_JN1967_registryB_surfaceC_dz_p0p02_sigma0p10_k6x6x1_v24`: verification `PASS`, SCF converged `True`.
- `AgC_JN1967_registryB_surfaceC_dz_m0p02_sigma0p10_k6x6x1_v24`: verification `PASS`, SCF converged `True`.
- `AgC_JN1967_registryB_surfaceC_dz_p0p02_sigma0p10_k6x6x1_v24`: verification `missing`, SCF converged `False`.

ACE fit summary not present; inspect the archived fit status and logs before using these labels.

## Limits and next step

This is an exploratory candidate. Passing these gates does not establish transfer to the 1M-atom MACE-derived QSGS trajectory. The four direct QSGS crops remain geometry-only until their own DFT center-force labels are calculated. If no candidate meets every registered gate, use the reported residuals to choose the next acquisition; do not treat a failed fit as a valid production potential.
