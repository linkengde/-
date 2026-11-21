# DFT acquisition and ACE/MTP exploratory fit report

Updated: 2026-11-21T14:35:39.551766+00:00 UTC

DFT sequence: `selected_three_archives_verified`; all 3/3 registered labels verified.

## Current exploratory gates

- Energy RMSE ≤10 meV/atom on training and periodic development data.
- Force-component RMSE ≤0.2 eV/Å (updated from 0.1 by the user).
- Si hotspot center-vector error ≤0.5 eV/Å; v11 center 01 remains development, not blind test.
- `gate_reevaluation_force0p2.json` recomputes these gates from existing fit metrics; no model was retrained for this threshold change.

## DFT labels

- `AgTi_JN1967_registryB_surfaceC_dz_p0p02_sigma0p10_k6x6x1_v24`: archive verification `PASS`, SCF converged in 51 iterations.
- `AgC_JN1967_registryB_surfaceC_dz_m0p02_sigma0p10_k6x6x1_v24`: archive verification `PASS`, SCF converged in 61 iterations.
- `AgC_JN1967_registryB_surfaceC_dz_p0p02_sigma0p10_k6x6x1_v24`: archive verification `PASS`, SCF converged in 65 iterations.

## ACE candidates

| Model | Dev energy (meV/atom) | Dev force (eV/Å) | Si center 01 error (eV/Å) | All gates |
|---|---:|---:|---:|:---:|
| `AgCSiTi_B1616_E10000_F1_a0.001_hotspot_center_weight1.ace.yaml` | 4.399 | 0.0446 | 6.767 | no |
| `AgCSiTi_B1616_E10000_F1_a0.001_hotspot_center_weight10.ace.yaml` | 4.399 | 0.0446 | 6.767 | no |
| `AgCSiTi_B1616_E10000_F1_a0.001_hotspot_center_weight100.ace.yaml` | 4.399 | 0.0446 | 6.767 | no |

## MTP candidates

| Model | Dev energy (meV/atom) | Dev force (eV/Å) | Si center 01 error (eV/Å) | All gates |
|---|---:|---:|---:|:---:|
| `MTP10_v24acq_CENTER1000_F005.mtp` | 10.565 | 0.2030 | 4116.403 | no |
| `MTP10_v24acq_CENTER1000_F010.mtp` | 5.385 | 0.1593 | 1966.498 | no |
| `MTP10_v24acq_CENTER1000_F020.mtp` | 14.423 | 0.1998 | 1705.614 | no |
| `MTP10_v24acq_CENTER100_F005.mtp` | 4.861 | 0.2169 | 4545.734 | no |
| `MTP10_v24acq_CENTER100_F010.mtp` | 6.656 | 0.1415 | 1582.127 | no |
| `MTP10_v24acq_CENTER100_F020.mtp` | 25.473 | 0.1795 | 24894.098 | no |
| `MTP10_v24acq_F005_SF0p0.mtp` | 30.481 | 0.2692 | 4406.605 | no |
| `MTP10_v24acq_F010_SF0p0.mtp` | 7.734 | 0.1622 | 43071.956 | no |
| `MTP10_v24acq_F020_SF0p0.mtp` | 9.307 | 0.1240 | 97894.230 | no |

## Decision and next step

No model passes all three current gates. Several MTP candidates meet the revised energy/force gates, but their Si center 01 errors remain far above 0.5 eV/Å. ACE candidates also miss the Si center gate. None is qualified for heating.

The next fit needs more representative Si–Ag interfacial hotspot labels. The registered blind test remains unopened.
