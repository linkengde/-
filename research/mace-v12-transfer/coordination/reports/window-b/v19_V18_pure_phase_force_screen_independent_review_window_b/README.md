# Independent reproduction of the V18 pure-phase force screen

## Scope

This is an independent CPU-only reproduction of A's V18 epoch-79 force screen on exactly two existing public GPAW numerical controls: `Ag_baseline_k8x8x8_v19` and `Ti3SiC2_baseline_k10x10x2_v19`. It verifies archive checksums, convergence and archive receipts, source/result identity, the selected model and checkpoint hashes, then recomputes MACE forces and compares the published aggregate, species and Cartesian-component errors.

The analysis did not open `train.extxyz`, development/test/holdout structures or sealed labels. It did not run DFT/MPI or training. These pure-phase controls are numerical diagnostics, not independent validation structures or evidence of production readiness. The Ti₃SiC₂ DFT geometry itself has a nonzero force RMS; the reported model error is the residual against that archived force field.

## Reproduce

From the repository root in the pinned MACE CPU environment:

```bash
/workspace/.venvs/mace-cpu/bin/python research/mace-v12-transfer/coordination/reports/window-b/v19_V18_pure_phase_force_screen_independent_review_window_b/reproduce_screen.py
```

The script reads only the two named public DFT archive directories, their exact source structures, A's published pure-phase `metrics.json`, the V18 selection record/artifact manifest, and the selected model/checkpoint. It checks the results against A's metrics using relative tolerance `1e-7` and absolute tolerance `1e-8` eV/Å. Machine-readable values and hash evidence are in `audit.json`; the compact comparison is in `reference_comparison.csv`.

## Results

| DFT numerical control | Atoms | DFT force RMS (eV/Å) | MACE force RMS (eV/Å) | Force-vector RMSE (eV/Å) | Maximum atom error (eV/Å) | x / y / z component RMSE (eV/Å) | 0.05 eV/Å screen |
|---|---:|---:|---:|---:|---:|---:|---|
| Ag, 8×8×8 | 32 | 0.00001728 | approximately 0 | 0.00001728 | 0.00002704 | 0.00000998 / 0.00000998 / 0.00000998 | PASS |
| Ti₃SiC₂, 10×10×2 | 48 | 0.05647788 | 2.93774295 | 2.88943327 | 4.25818832 | 0.00124752 / 0.00089351 / 2.88943286 | FAIL |

All recorded aggregate, per-species and Cartesian-component comparisons reproduce A's report within the stated tolerance. For Ti₃SiC₂, species force-vector RMSE is **3.47676 eV/Å for Ti**, **2.62959 eV/Å for C**, and **0.00135 eV/Å for Si**. The residual is almost entirely in z. Ag passes this single frozen fcc control; Ti₃SiC₂ fails the force screen by a wide margin.

This confirms the published numerical comparison and points to a strong Ti/C framework response error in this one Ti₃SiC₂ geometry. It does not establish a microscopic cause or dynamical instability. Combined with the separately localized AgSi residuals, it shows that Ag–Si contact forces alone cannot account for the failed AgSi force screen.

## Lineage and archive checks

- Selected model SHA-256: `757586671174db85acbf7971ccfe92de19cc58378ee2431a797efa6d6163cb62`
- Epoch-79 checkpoint SHA-256: `d9faab45f32df7366b32f44a9f91c1aab9c72f2df422f24c9996fa3243b2c137`
- Both named archive `SHA256SUMS.txt` manifests passed; each summary reports SCF convergence and each recorded archive verification is `PASS`.
- Source structure hashes and exact atom numbers/order, IDs, positions, cell and PBC match their archived results and A's report.
- All V18 lineage checks for the selected model/checkpoint passed. The training file itself was not accessed.

