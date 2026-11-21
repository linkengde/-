# DFT acquisition and ACE/MTP exploratory fit report

Updated: 2026-11-21T14:14:53.079636+00:00 UTC

DFT sequence: `selected_three_archives_verified`; verified labels: 3/3.

## Registered exploratory gates

- Energy RMSE ≤10 meV/atom.
- Force-component RMSE ≤0.1 eV/Å.
- Registered Si hotspot center-vector error ≤0.5 eV/Å.
- v11 index 01 is development, not blind test; the blind test remains unopened.

## DFT labels

- `AgTi_JN1967_registryB_surfaceC_dz_p0p02_sigma0p10_k6x6x1_v24`: archive verification `PASS`, SCF converged `True` in 51 iterations.
- `AgC_JN1967_registryB_surfaceC_dz_m0p02_sigma0p10_k6x6x1_v24`: archive verification `PASS`, SCF converged `True` in 61 iterations.
- `AgC_JN1967_registryB_surfaceC_dz_p0p02_sigma0p10_k6x6x1_v24`: archive verification `PASS`, SCF converged `True` in 65 iterations.

## ACE candidates

No ACE candidate passed all gates. Energy and force RMSE passed, but every candidate missed the held-out Si-center error gate.

| Candidate | Energy (meV/atom) | Force (eV/Å) | Si-center error (eV/Å) | Pass |
|---|---:|---:|---:|:---:|
| `AgCSiTi_B1616_E10000_F1_a0.001_hotspot_center_weight1.ace.yaml` | 4.3989 | 0.04464 | 6.7673 | no |
| `AgCSiTi_B1616_E10000_F1_a0.001_hotspot_center_weight10.ace.yaml` | 4.3988 | 0.04464 | 6.7674 | no |
| `AgCSiTi_B1616_E10000_F1_a0.001_hotspot_center_weight100.ace.yaml` | 4.3988 | 0.04464 | 6.7674 | no |

## MTP candidates

The 9-candidate force-weight/center-weight sweep completed. None passed the force and held-out hotspot gates.

| Candidate | Dev energy (meV/atom) | Dev force (eV/Å) | Si-center 01 error (eV/Å) | Pass |
|---|---:|---:|---:|:---:|
| `MTP10_v24acq_F005_SF0p0.mtp` | 30.481 | 0.2692 | 4406.605 | no |
| `MTP10_v24acq_F010_SF0p0.mtp` | 7.734 | 0.1622 | 43071.956 | no |
| `MTP10_v24acq_F020_SF0p0.mtp` | 9.307 | 0.1240 | 97894.230 | no |
| `MTP10_v24acq_CENTER100_F005.mtp` | 4.861 | 0.2169 | 4545.734 | no |
| `MTP10_v24acq_CENTER100_F010.mtp` | 6.656 | 0.1415 | 1582.127 | no |
| `MTP10_v24acq_CENTER100_F020.mtp` | 25.473 | 0.1795 | 24894.098 | no |
| `MTP10_v24acq_CENTER1000_F005.mtp` | 10.565 | 0.2030 | 4116.403 | no |
| `MTP10_v24acq_CENTER1000_F010.mtp` | 5.385 | 0.1593 | 1966.498 | no |
| `MTP10_v24acq_CENTER1000_F020.mtp` | 14.423 | 0.1998 | 1705.614 | no |

## Decision and next acquisition

No current ACE or MTP model is qualified for heating. Do not use these candidates for production trajectories. The dominant remaining error is transfer to the distinct Si hotspot environment. Add representative high-force Si local environments to the training set and retain an independent hotspot label for validation before another qualification attempt.
