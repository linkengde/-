# V18 AgSi force-error localization

## Scope

This analysis reproduces the completed V18 CPU screen on only the two archive-verified public AgSi CIF references and localizes its force residuals. It reads `mace_periodic_v18_clean_core/data/train.extxyz` only for local Ag–Si contact coverage. Development/test structures and sealed labels remain unopened. No training, model edits, DFT/MPI, dataset integration, or production/gap-scan changes are made.

Both references share one CIF source family, and the marked Ag–Si distance lies near a training contact. This is a single-system diagnostic, not independent validation or evidence of model readiness.

## Reproduce

From the repository root, using the existing pinned environment:

```bash
/workspace/.venvs/mace-cpu/bin/python research/mace-v12-transfer/coordination/reports/window-b/v19_AgSi_force_error_localization_window_b/localize_forces.py
```

The script checks the prior inference record, V18 model/checkpoint/train/manifest hashes, both result archives and their exact source geometries, and reproduces the published all-atom and marked-pair errors before writing localized outputs.

## Metrics

- Per-atom vector error is `norm(F_model - F_DFT)` in eV/Å; all-atom force-vector RMSE is `sqrt(mean_i(sum_xyz(delta_F_i**2)))`.
- Per-species rows include vector RMSE, maximum and mean vector error, component RMSE, and each species' share of total squared error.
- The signed marked-pair force is `dot(F_X - F_Ag, minimum_image_unit_vector_Ag_to_X)`. Net Ag-layer normal force is the sum of Fz over all Ag atoms, positive along +z.
- Top-error neighborhoods list species counts and distances within 4.5 Å, including nearby Ag–X pairs. The training contact coverage uses only marked pairs in `train.extxyz`.
- Geometry probes, if supported by a residual site recurring in the top eight of both references with aligned error direction, are ±0.03 Å single-atom perturbations. They are unowned, unlabeled, launch-disabled report artifacts requiring review and registration before any future acquisition.

Provisional gates remain force-vector RMSE ≤0.05 eV/Å and absolute marked-pair separation-force error ≤0.10 eV/Å. Any suggested geometry remains diagnostic-only and cannot be used as a label or training example.

## Findings

Both references reproduce the prior diagnostic and fail both gates. The largest force errors are on framework Ti atoms, not the Ag layer. In the compact pilot, Ti contributes 64.88% and C 33.18% of total squared force error; Si contributes 1.56% and Ag 0.38%. The vacuum control gives the same distribution within rounding. Residuals are almost entirely along z: top Ti errors are about 4.558 eV/Å with one sign on one Ti layer and about 4.431 eV/Å with the opposite sign on the adjacent layer. The overall force-vector RMSE is 2.63375 eV/Å for the pilot and 2.63351 eV/Å for the vacuum control; marked-pair separation-force absolute errors are 1.71714 and 1.71553 eV/Å.

| Reference | All-atom force-vector RMSE (eV/Å) | Ag–Si separation force, DFT → MACE (eV/Å) | Absolute separation-force error (eV/Å) | Ag-layer net Fz error, MACE − DFT (eV/Å) |
|---|---:|---:|---:|---:|
| Compact pilot | 2.63375 | −1.53713 → −3.25427 | 1.71714 | −2.28285 |
| 26 Å vacuum control | 2.63351 | −1.53874 → −3.25427 | 1.71553 | −2.28025 |

Both force gates fail in both references (limits: 0.05 eV/Å all-atom RMSE and 0.10 eV/Å absolute Ag–Si separation-force error). Compact-pilot species vector RMSEs are Ti 3.12270, C 2.73490, Si 0.83924, and Ag 0.58442 eV/Å, consistent with the framework-dominated residual distribution.

For top-error Ti centers, the 4.5 Å neighborhood contains 5 C, 4 Si, and 8 Ti atoms, with nearest distances about 2.088 Å (C), 2.681 Å (Si), and 2.970 Å (Ti). There is no Ag or Ag–X pair within 4.5 Å; the nearest Ag is farther away (about 9.4–13.5 Å for the listed top sites). Thus the largest residuals arise in the Ti₃SiC₂ framework away from the Ag contact. The marked Ag–Si pair still has a large signed separation-force error, so the interface force response also fails.

The script writes four label-free geometry probes: ±0.03 Å along the reproducible z-residual direction on one representative Ti from each of the two opposite-sign layers (IDs 19101002 and 19101004). They are stored only under this report directory with `launch_enabled=false`, no owner, no energy/force labels, and `training_eligible=false`. They do not modify the registered gap scan or authorize any DFT job.
