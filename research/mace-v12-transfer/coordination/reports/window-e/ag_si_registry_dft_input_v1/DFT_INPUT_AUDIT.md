# Ag–Si registry DFT input audit

Status: `READY_FOR_DFT_REGISTRY_SCREENING`  
DFT execution: `HOLD` (not authorized by the source manifest and not run in this task).

## Frozen source and role

The inputs were read byte-for-byte from source commit `9a5482d93a8756850ab457a605cf75dfa1f8a3d4`, directory `research/mace-v12-transfer/coordination/reports/window-c/ag_si_dev_v2_structure_generation`. The source manifest SHA-256 is `a194e13815271ce1c99b11d41019c6db0d08adfe28da621bfb6c5f4883c62724`. Each original XYZ is copied byte-for-byte into its candidate directory; the source registry and all source files remain unchanged.

The source manifest says `GEOMETRY_ONLY_NOT_ASSIGNED_TRAIN_DEV_TEST`, `independent_dev_family=false`, and `dft_execution_authorized=false`. These screening candidates are not promoted to DEV or any other data role by this package. All four derive from the COD 9009647/Kisi 1998 parent family.

## Candidate integrity and geometry

ASE 3.29.0 parsed the exact XYZ bytes. Composition is 7 Ag, 4 C, 3 Si, 6 Ti (20 atoms; formula C4Ag7Si3Ti6). The cell matrix, PBC, vacuum values and source hashes below are recorded from those files. Vacuum is measured from the fixed coordinates to the z cell boundaries; no atom or cell was moved. A warning is required only when either value is strictly below 18 Å.

| Registry | Atoms | Composition | Source XYZ SHA-256 | Cell (Å) | PBC | Bottom z vacuum (Å) | Top z vacuum (Å) | Vacuum check |
|---|---:|---|---|---|---|---:|---:|---|
| R1 | 20 | Ag7 C4 Si3 Ti6 | `b758b39571a6ce8b24f687e9b341c4328336c108d761b4adaa25e3fb51932890` | 3.071 × 3.071 × 71.268266; γ=120° | T,T,F | 18.000000000000 | 18.000000000174 | none |
| R2 | 20 | Ag7 C4 Si3 Ti6 | `cb2e16376887f7cdeb05c7ddf4aacfacacbea0a000837a3888a0f3ea91785f70` | 3.071 × 3.071 × 71.268266; γ=120° | T,T,F | 18.000000000000 | 18.000000000174 | none |
| R3 | 20 | Ag7 C4 Si3 Ti6 | `edc7a00668be7fb0f76616463d3eacf2ac53d0ea8837b21fb9fb3c66eafa6f76` | 3.071 × 3.071 × 71.268266; γ=120° | T,T,F | 18.000000000000 | 18.000000000174 | none |
| R4 | 20 | Ag7 C4 Si3 Ti6 | `6465b1e2ea8c675a4d368a6f385cc31302a5c4b017d9b39842e82c93c196210c` | 3.071 × 3.071 × 71.268266; γ=120° | T,T,F | 18.000000000000 | 18.000000000174 | none |

All candidates have `cell=[[3.071,0,0],[-1.5355,2.6595640150220112,0],[0,0,71.26826600017428]] Å` and PBC `[true,true,false]`. The z vacuum is 18.000000000000 Å at the bottom and 18.000000000174 Å at the top, so none triggers the `<18 Å` warning rule.

## Locked single-point preparation

All four `input.py` files implement the same GPAW/PBE PAW plane-wave setup and preserve the frozen geometry:

- GPAW 26.7.0, ASE 3.29.0, gpaw-data 1.2.1; local Ag, C, Si, Ti PAW hashes are locked in `PARAMETER_LOCK.json`.
- PBE; plane-wave cutoff 500 eV; Fermi-Dirac smearing 0.10 eV.
- SCF thresholds: energy, density and eigenstates each `1e-5`; max 160 iterations; mixer beta 0.05, nmaxold 8, weight 100.
- Single point only: atom relaxation and cell relaxation are both disabled. The scripts assert source hash, atom count, element counts, cell and PBC before any calculation; they assert positions, cell and PBC remain byte-identical in numeric value before writing results.
- The requested output labels are total energy, force-consistent free energy, atomic forces, and the full 3×3 ASE tensile-positive stress tensor in eV/Å³. Slab stress includes the vacuum cell volume; stress is marked `stress_supervision_allowed=false` pending C review.
- Coarse k-point mesh: 6×6×1. Medium mesh: 8×8×1. Both are prepared only; neither is run or marked as production-approved. The scripts use the GPAW default Monkhorst-Pack interpretation of integer mesh tuples, with no explicit gamma override or extra offset; this convention is included in any later C release gate.

Default invocation of each `input.py` performs only frozen-input validation and prints the plan. A GPAW calculation additionally requires an explicit `--run`, a matching SHA-256 supplied for a C release JSON, an affirmative execution authorization, C approval identity, and exact matching task, geometry and method fields. No release is included here.

## Verification and non-actions

The four source SHA values match both the source `SHA256SUMS.txt` and `MANIFEST.json`. The generated `R1/structure.xyz` through `R4/structure.xyz` are exact copies. `SHA256SUMS.txt` covers every generated `input.py`, structure copy, this audit and the parameter lock.

No DFT, geometry optimization, registry edit, cell edit, coordinate edit, or MACE training was performed. No DFT output or computed label is included.
