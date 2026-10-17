# Ti3SiC2 P0 fixed-geometry DFT diagnostics

All six authorized structures completed with SCF convergence. They remain `DIAGNOSTIC_ONLY`; TRAIN/DEV/TEST roles were not assigned or changed. Geometry and atom ordering in each labelled output exactly match its frozen input.

Method lock: GPAW 26.7.0, ASE 3.29.0, PBE/PAW, 500 eV plane-wave cutoff, Fermi-Dirac 0.10 eV, Monkhorst-Pack 6×6×2, SCF thresholds 1e-5, fixed geometry, native MSR1, MPI 4, OMP 1.

Execution lock SHA-256: `81d9ce869d153b3bde909f3d8259a3b89a437a8150158dcd32211e9a911afd0f`.

The positive Si displacement completed in a separate attempt directory after its first launcher invocation stopped at preflight because `GPAW_SETUP_PATH` was unset. That first failure produced no SCF data and is preserved. The successful attempt used the already installed PAW setup path and the same scientific settings.

| Structure | Status | SCF steps | Total energy (eV) | Free energy (eV) | Max force (eV/Å) |
|---|---|---:|---:|---:|---:|
| TSC48_VOL_m02 | SUCCESS | 16 | -389.386851657 | -389.627349554 | 0.123827 |
| TSC48_VOL_p02 | SUCCESS | 17 | -389.939807472 | -390.192354596 | 0.177615 |
| TSC48_BASAL_XY_SHEAR_m005 | SUCCESS | 16 | -389.811744458 | -390.058464482 | 0.084977 |
| TSC48_BASAL_XY_SHEAR_p005 | SUCCESS | 17 | -389.810461087 | -390.056795128 | 0.086613 |
| TSC48_SI_ID7_Z_m020 | SUCCESS | 24 | -389.815059678 | -390.061516238 | 0.266821 |
| TSC48_SI_ID7_Z_p020 | SUCCESS | 24 | -389.815059520 | -390.061516140 | 0.266820 |

Stress tensors are recorded in GPa in `RESULT.json` and the CSV. Total energy and force-consistent free energy are reported separately. Each task directory contains its input, labelled output, GPAW log, result JSON, command/status records, and verified `SHA256SUMS.txt`.

Limitations: these are diagnostic single points on the six frozen P0 structures. They do not establish model generalization or production suitability; strained-cell stress convergence and independent validation families remain separate checks. No TRAIN/DEV/TEST assignment, fitting, MD, or additional DFT was performed.
