# v9 Ag–Ti existing-point relative-energy diagnostic

This is an uncalibrated, post-hoc diagnostic on three existing periodic PW-PBE structures. No model parameters or energy zero were fitted. All structures have the same 49 atoms, composition C11Ag28SiTi9, 28 Å cubic cell, PBC, and DFT method (PBE, 500 eV, 1×1×1 k-points).

| Separation (Å) | DFT E (eV) | v9 E (eV) | DFT ΔE from 2.70 Å (eV) | v9 ΔE (eV) | ΔE error (meV/cell) |
|---:|---:|---:|---:|---:|---:|
| 2.3871 | -176.730775269 | -173.018175846 | 0.284045693 | 0.298411580 | 14.366 |
| 2.5500 | -176.987528749 | -173.290852864 | 0.027292214 | 0.025734562 | -1.558 |
| 2.7000 | -177.014820962 | -173.316587425 | 0.000000000 | 0.000000000 | 0.000 |

Caveats:
- 2.3871 Å was included in v9 force training, while its energy label was excluded; this is not an independent structural holdout.
- 2.70 Å is a PW-PBE method-transfer point on a geometry represented in LCAO force data.
- 2.55 Å is a new distance in the same motif; one interval is insufficient to validate a whole energy curve.
- These relative energies only probe the sampled Ag–Ti coordinate in one motif. They do not validate Ag–Si/Ag–C binding energies, thermal disorder, energy conservation under MD, or formation energies across compositions.
- No calibration was applied. Element-reference calibration would require a training-only fit; validation energies must stay blind.

Model SHA256: `7d183ed99dae18bfb0f84eaf33be38458671d18508ff59e9413b7339bd97db86`
