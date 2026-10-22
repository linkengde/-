# (11-20) prismatic Ti₃SiC₂: frozen TSCBOP force comparison

## Scope and result

This report compares the two converged 120-atom PBE/PAW fixed-geometry DFT force sets against the unchanged TSCBOP potential using LAMMPS `run 0`. It performs no minimization, MD, DFT, or parameter fitting. The two geometries are the reference slab and the same slab with source atom IDs 4 and 113 displaced by −0.02 and +0.02 Å along z, respectively.

The absolute TSCBOP force disagreement is large: vector RMSE is 1.9452 eV/Å for REF and 1.9448 eV/Å for RESP. It is concentrated in Ti and C; Si force errors are much smaller. For this one paired displacement, the *change* in TSCBOP forces tracks the DFT force change more closely than the absolute forces do, but the global differential-response vector RMSE remains 0.01117 eV/Å, and the local value around the displaced Ti atoms is 0.02678 eV/Å. This single directional finite response does not demonstrate general surface accuracy.

## Input and implementation checks

- REF input SHA-256: `fd5ff5c1cae1490bcbde9ffd08f329c7ec8cac687555f29568d725f4230cc5c3`; geometry SHA-256: `e876359dc190632d02ca2630bbd145afb4b866ddd19de65c8e1019c2e6f6c1a4`; DFT SCF converged in 44 steps.
- RESP input SHA-256: `6f91573bcb4a0feb7e9f05329f6eca422d3e7aac723763e961b52c7b83ddf9b6`; geometry SHA-256: `817e65667a2ab51194cc69c9800b1efc4dcadfc49e7a3dbf3867a8577df54fb3`; DFT SCF converged in 49 steps.
- Both structures contain 120 atoms (Ti₆₀Si₂₀C₄₀). Input and labelled-output symbols, positions, cell, PBC, and source atom IDs match. DFT forces were read from the converged `RESULT.json` in source-atom order and checked for shape and finite values.
- LAMMPS 22 Jul 2025 Update 4; `units metal`; `pair_style tersoff`; `pair_coeff * * TSCBOP.tersoff Ti Si C`; `OMP_NUM_THREADS=1`; exactly `run 0`. The force unit is eV/Å. The source potential SHA-256 is `729273a6488d75756ed7927432c77c7347a1b115cecbaefae4981bafcec80fd7`, and each run used a byte-identical copy.
- The maximum difference between LAMMPS dump coordinates and the wrapped source coordinates was `1.42×10⁻¹⁴ Å`. Atom IDs and elements matched exactly. No absolute DFT-versus-LAMMPS energy comparison is made.

Errors below are `F_TSCBOP − F_DFT`. Component RMSE is over all x/y/z error components. Vector RMSE is `sqrt(mean_i(sum_xyz(error_i²)))`; P95 and maximum are percentiles/maxima of each atom's vector-error magnitude.

## Absolute force comparison

| Geometry | Element | N | Component RMSE (eV/Å) | Vector RMSE (eV/Å) | P95 vector error (eV/Å) | Maximum (eV/Å) | Source ID at maximum |
|---|---:|---:|---:|---:|---:|---:|---:|
| REF | all | 120 | 1.1231 | 1.9452 | 2.8112 | 2.8113 | 29 |
| REF | Ti | 60 | 1.2442 | 2.1551 | 2.8112 | 2.8113 | 29 |
| REF | Si | 20 | 0.1250 | 0.2165 | 0.4073 | 0.4073 | 8 |
| REF | C | 40 | 1.2057 | 2.0884 | 2.1742 | 2.1742 | 46 |
| RESP | all | 120 | 1.1228 | 1.9448 | 2.8086 | 2.8129 | 30 |
| RESP | Ti | 60 | 1.2443 | 2.1552 | 2.8118 | 2.8129 | 30 |
| RESP | Si | 20 | 0.1251 | 0.2166 | 0.4120 | 0.4120 | 8 |
| RESP | C | 40 | 1.2050 | 2.0871 | 2.1741 | 2.1754 | 46 |

The frozen slab is not relaxed: its DFT maximum force is about 1.713 eV/Å. The large pointwise force discrepancy should therefore be read as a same-geometry force mismatch for this candidate, not as an equilibrium property or a calibrated validation score.

## Paired displacement response

For each source ID, response error is `(F_TSCBOP_RESP − F_TSCBOP_REF) − (F_DFT_RESP − F_DFT_REF)`. The local subset comprises 18 atoms within 3.3 Å (minimum-image distance in the REF cell) of either displaced terminal Ti.

| Region | N | Component RMSE (eV/Å) | Vector RMSE (eV/Å) | P95 vector error (eV/Å) | Maximum (eV/Å) | Source ID at maximum |
|---|---:|---:|---:|---:|---:|---:|
| All atoms | 120 | 0.00645 | 0.01117 | 0.02425 | 0.06359 | 22 |
| Within 3.3 Å of displaced Ti | 18 | 0.01546 | 0.02678 | 0.06359 | 0.06359 | 22 |

The DFT/TSCBOP response-vector RMS magnitudes are 0.03094/0.03273 eV/Å over all atoms, and 0.07964/0.08411 eV/Å in the 3.3 Å local subset. By element, the differential-response error vector RMSE is 0.00773 eV/Å for Ti, 0.00564 eV/Å for Si, and 0.01640 eV/Å for C. These are finite differences for one antisymmetric displacement mode, not independent derivatives for the two Ti atoms.

## Interpretation and limits

The near agreement in this paired force response coexists with roughly 1.94 eV/Å absolute force-vector RMSE. This may indicate that a local response slope is captured while the baseline force field is offset on the unrelaxed geometry; the present pair cannot establish the cause or support a parameter change. The largest absolute errors occur in Ti and C, whereas Si is substantially closer on these two structures.

Both structures descend from COD 9009647 and are one geometry family. The facet is a distinct modeled prismatic orientation, but it is not an independent crystallographic-source DEV/TEST family. No data roles changed. TSCBOP, Ag EAM, and any Ag cross interaction were not modified. This result does not validate Ag interfaces, other surface terminations, finite-temperature behavior, or production use.

## Artifacts

- `SUMMARY.json`: machine-readable metrics and input bindings.
- `REF_FORCE_COMPARISON.csv`, `RESP_TIZ020_FORCE_COMPARISON.csv`: per-atom DFT/TSCBOP forces and errors keyed by source atom ID.
- `RESPONSE_DIFFERENCE.csv`: same-ID finite response comparison.
- `runs/REF/` and `runs/RESP_TIZ020/`: exact LAMMPS input, data, potential copy, force dump, and logs.
- `run_pair.py`: reproducible comparison runner.
- `SHA256SUMS.txt`: hashes for the complete report package.
