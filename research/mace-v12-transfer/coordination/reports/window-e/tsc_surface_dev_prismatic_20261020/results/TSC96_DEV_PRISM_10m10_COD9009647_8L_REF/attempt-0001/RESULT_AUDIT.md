# Prismatic Ti3SiC2 DEV candidate reference point

## Execution result

The fixed-geometry reference single point completed successfully on 2026-10-20. GPAW reported `Converged in 46 steps`; the runner wrote `RESULT.json`, `labelled.extxyz`, and `RUN_STATUS.json` only after checking the convergence marker and validating the unchanged coordinates, cell, PBC, and atom order.

| Quantity | Result |
|---|---:|
| Structure | `TSC96_DEV_PRISM_10m10_COD9009647_8L_REF` |
| Candidate role | `DEV_GEOMETRY_CANDIDATE_REFERENCE` (role not assigned or changed) |
| Source family | `surface_prismatic_10m10_cod9009647_8layer` |
| Atoms / composition | 96 / Ti48 Si16 C32 |
| Input SHA-256 | `60096a122270e57ac2fb1cbbc2c0d1b087eab12127bcee087f30d244dc8900e1` |
| Geometry SHA-256 | `30a4b161cd9a912c90006f3c0d35978e4d417306734f3fc982c9c8684fa600d3` |
| PBE total energy (extrapolated) | −752.5779427569561 eV |
| Force-consistent free energy | −753.3027126030335 eV |
| Maximum force magnitude | 2.0111018439658994 eV/Å |
| SCF | converged, 46 iterations |
| Stress | available; see `RESULT.json` (eV/Å³ and GPa) |
| Peak cgroup memory | 16,583,155,712 bytes |
| Free workspace at completion | 14,219,640,832 bytes |

The nonzero maximum force is retained as the fixed-geometry surface response. The geometry was not relaxed and the force magnitude is not treated as evidence of a failed SCF.

## Locked calculation settings

GPAW 26.7.0, ASE 3.29.0, Python 3.12.14, PBE, plane waves at 500 eV, Fermi–Dirac width 0.10 eV, standard Monkhorst–Pack 12×2×1 with no extra Gamma flag or shift, non-spin-polarized, fixed coordinates and cell, SCF thresholds 1e-5, maximum 160 iterations, GPAW default Davidson (`niter=2`), MSR1 (`beta=0.05`, `nmaxold=10`, `weight=20`, remaining locked values in `SCIENTIFIC_CONFIG.json`), 4 MPI ranks and 1 OMP thread. Ti/Si/C PAW hashes and the configuration hash are recorded in the committed scientific configuration and `RESULT.json`.

Configuration SHA-256: `c1b19b90a091d735bbbc6ebb184fb31ec1f7d3985c75f952879bb43bfda85fba`.

## Validation and limits

- `RUN_STATUS.json` and `RESULT.json` both report `SUCCESS` and `scf_converged=true`.
- The GPAW log contains the matching `Converged in 46 steps` marker.
- Exactly 96 finite force vectors and a finite 3×3 stress tensor were written.
- The labelled structure has 96 atoms with Ti48 Si16 C32; coordinates, cell, PBC, and ordering exactly match the frozen input.
- The pair is a geometry-family DEV candidate, not an independent crystallographic-source blind test: it shares COD 9009647 with the originating bulk family. No TRAIN/DEV/TEST manifest role was altered and no TEST identity or label was accessed.
- This is one fixed-geometry reference point only. It does not by itself establish an independent DEV family; the paired deterministic displacement point remains to be completed.

## Checkpoint disposition

The last full-wavefunction checkpoint was written at SCF step 40 and its on-disk SHA-256 matched the runner-recorded value `e503e04a3ac48066de3e66612297c768368b04d4792073af9938622173ecf8e0`. SCF subsequently converged and all labels were written, so the checkpoint was deleted to recover temporary workspace. It was not read back or used for a restart; no restart is required for this successful point. The hash and disposition are recorded in `CHECKPOINT_STATUS.json` and `CHECKPOINT_DISPOSITION.json`.
