# DFT acquisition and ACE refit report

Updated: 2026-11-21T11:34:19.030825+00:00
DFT sequence status: `selected_three_archives_verified`
ACE pipeline status: `fit_complete_needs_review`

## Registered exploratory gates

- Energy RMSE: ≤ 10 meV/atom on full-label train and periodic development.
- Force-component RMSE: ≤ 0.1 eV/Å on full-label train and periodic development.
- Registered Si hotspot center-force vector error: ≤ 0.5 eV/Å for training centers and the v11 development center.
- The v11 development frame is not a blind test.

## Verified DFT labels

- `AgTi_JN1967_registryB_surfaceC_dz_p0p02_sigma0p10_k6x6x1_v24`: verification `PASS`, SCF converged `True`.
- `AgC_JN1967_registryB_surfaceC_dz_m0p02_sigma0p10_k6x6x1_v24`: verification `PASS`, SCF converged `True`.
- `AgC_JN1967_registryB_surfaceC_dz_p0p02_sigma0p10_k6x6x1_v24`: verification `PASS`, SCF converged `True`.

## ACE candidates

Passing models: none
- `AgCSiTi_B1616_E10000_F1_a0.001_hotspot_center_weight1.ace.yaml`: gate_pass=False; periodic-dev energy=4.3989366732541715 meV/atom; force-component RMSE=0.044635648187655474 eV/Å; Si-center vector error=6.7672534418853205 eV/Å.
- `AgCSiTi_B1616_E10000_F1_a0.001_hotspot_center_weight10.ace.yaml`: gate_pass=False; periodic-dev energy=4.398809367614074 meV/atom; force-component RMSE=0.04463936047057714 eV/Å; Si-center vector error=6.767363785463602 eV/Å.
- `AgCSiTi_B1616_E10000_F1_a0.001_hotspot_center_weight100.ace.yaml`: gate_pass=False; periodic-dev energy=4.398795012252506 meV/atom; force-component RMSE=0.04463975490415603 eV/Å; Si-center vector error=6.767376049522513 eV/Å.

## Limits and next step

This is an exploratory candidate. Passing these gates does not establish transfer to the 1M-atom MACE-derived QSGS trajectory. The four direct QSGS crops remain geometry-only until their own DFT center-force labels are calculated. If no candidate meets every registered gate, use the reported residuals to choose the next acquisition; do not treat a failed fit as a valid production potential.
