# V18 pure-phase force screen

## Scope

The frozen V18 epoch-79 model was evaluated on two existing, checksum-verified GPAW/PBE numerical-control archives: the Ag fcc baseline at 8×8×8 and the Ti₃SiC₂ baseline at 10×10×2. This was CPU-only inference. No new DFT/MPI, fitting, training, model edits, development/test structures, or sealed labels were used. Archive and input-geometry checks are recorded in `metrics.json`.

These controls diagnose transfer to pure phases. They are not independent validation structures or evidence of production readiness. The Ti₃SiC₂ DFT reference has a nonzero force RMS of 0.05648 eV/Å, so the meaningful comparison is the model-minus-DFT residual; the screen is not a relaxation or finite-temperature stability test.

## Results

| DFT control | Atoms | DFT force RMS (eV/Å) | V18 force RMS (eV/Å) | Force-vector RMSE (eV/Å) | 0.05 eV/Å screen |
|---|---:|---:|---:|---:|---|
| Ag, 8×8×8 | 32 | 0.0000173 | approximately 0 | 0.0000173 | PASS |
| Ti₃SiC₂, 10×10×2 | 48 | 0.05648 | 2.93774 | 2.88943 | FAIL |

For Ti₃SiC₂, the per-species force-vector RMSE is 3.47676 eV/Å for Ti, 2.62959 eV/Å for C, and 0.00135 eV/Å for Si. The error is almost entirely in the z component (2.88943 eV/Å); x and y component RMSEs are 0.00125 and 0.00089 eV/Å. Ag passes this one frozen fcc geometry, while the model fails badly on the Ti₃SiC₂ backbone. Together with the AgSi residual localization, this shows the failed AgSi force screen is not explained by the Ag–Si contact alone. Correcting only that cross term would leave a large framework-force error.

The Ti₃SiC₂ force screen does not prove the microscopic cause or demonstrate dynamical instability. It does mean this V18 model should not be used to assess framework stability or to interpret interface forces until the Ti/C framework response is corrected and retested. No gap-scan input was enabled by this work.

## Reproduction

Run from the repository root in the pinned CPU MACE environment:

```bash
/workspace/.venvs/mace-v12/bin/python research/mace-v12-transfer/coordination/reports/window-a/v19_V18_pure_phase_force_screen/infer_pure_phase.py
```

The script verifies the selected V18 model lineage, archive checksum manifests, SCF convergence and recorded archive verification, source CIF-derived geometry hashes, and exact atom ordering before inference. It rewrites `metrics.json` only after both force predictions are finite.

The selected model SHA-256 is `757586671174db85acbf7971ccfe92de19cc58378ee2431a797efa6d6163cb62`. DFT result hashes, input hashes, method settings, per-species force errors, and runtime versions are in `metrics.json`.
