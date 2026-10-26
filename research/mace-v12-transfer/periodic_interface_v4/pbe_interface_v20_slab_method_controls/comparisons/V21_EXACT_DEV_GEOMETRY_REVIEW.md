# V21 exact-development-geometry slab method review

All four fixed-geometry GPAW PBE/PW500 controls converged and passed archive verification. Each was compared with the approved `AgSi_v19_independent_dev_01` baseline using the exact same ordered atoms, coordinates, cell, PBC, and geometry-file SHA256. These are numerical-method controls only; they are not training labels or independent validation.

Baseline: `AgSi_v19_independent_dev_01`, PBE/PW500, 6×6×1, sigma 0.10 eV. Shared input geometry SHA256: `bca0a8fc1a2ef4e83eebe62d18f0b77480e33692be8886d9910dfd8622dd4b92` (c = 34.2058765 Å; pbc = true,true,false).

| Control | k-points | Sigma (eV) | SCF iterations | Native ΔE (meV/atom) | Free-energy Δ (meV/atom) | Force RMS Δ (eV/Å) | Ag–Si projection Δ (eV/Å) | Budget result |
|---|---:|---:|---:|---:|---:|---:|---|
| `AgSi_dev34_k6x6x1_sigma0p05_v21` | 6×6×1 | 0.05 | 53 | -0.04073 | -5.65683 | 0.01495 | 0.02619 | `SENSITIVITY_EXCEEDS_BUDGET` |
| `AgSi_dev34_k6x6x1_sigma0p20_v21` | 6×6×1 | 0.20 | 55 | 0.60050 | 22.05047 | 0.03001 | -0.04723 | `SENSITIVITY_EXCEEDS_BUDGET` |
| `AgSi_dev34_k4x4x1_sigma0p10_v21` | 4×4×1 | 0.10 | 58 | 0.64383 | 0.53755 | 0.01006 | -0.01574 | `SENSITIVITY_EXCEEDS_BUDGET` |
| `AgSi_dev34_k8x8x1_sigma0p10_v21` | 8×8×1 | 0.10 | 63 | 0.18029 | 0.13687 | 0.00282 | -0.00174 | `PASS_WITHIN_PROVISIONAL_BUDGET` |

Provisional limits: 2 meV/atom for native and free energies, 0.01 eV/Å for force-vector RMS, and 0.02 eV/Å for the marked Ag–Si pair projection.

## Interpretation

The 8×8×1, sigma 0.10 control is within all four provisional limits relative to the 6×6×1 baseline. The 4×4×1 control exceeds the force-RMS limit narrowly (0.01006 versus 0.01000 eV/Å); its energy and marked-pair checks pass. Sigma 0.05 and 0.20 each exceed the free-energy, force-RMS, and marked-pair limits, although their native-energy differences remain within the energy limit.

Retain 6×6×1, sigma 0.10 eV as the practical reference for this frozen development geometry: the higher 8×8×1 mesh adds little change within the selected budgets, while the lower mesh and alternative smearing widths do not meet all limits. This is a provisional, geometry-specific method choice, not proof of absolute convergence or universal optimality.

The original V20 controls used a different c = 44.2058765 Å cell; these V21 controls use the exact V19 development cell. All controls remain excluded from training and validation, and no approved dataset labels were changed.

## Machine-readable evidence

- `comparison_AgSi_dev34_k6x6x1_sigma0p05_v21__AgSi_v19_independent_dev_01.json`
- `comparison_AgSi_dev34_k6x6x1_sigma0p20_v21__AgSi_v19_independent_dev_01.json`
- `comparison_AgSi_dev34_k4x4x1_sigma0p10_v21__AgSi_v19_independent_dev_01.json`
- `comparison_AgSi_dev34_k8x8x1_sigma0p10_v21__AgSi_v19_independent_dev_01.json`
