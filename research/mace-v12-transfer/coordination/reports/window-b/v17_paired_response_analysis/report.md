# V17 paired response analysis

## Scope and status

This review compares the four matched positive/negative V17 diagnostic pairs and their three V14 training parents. It evaluates the frozen V16 checkpoint on those eight diagnostic geometries and the parent structures. It does not read any V17 withheld-test label, train or select a model, or launch DFT, MD, or TTM.

All eight diagnostic archives pass their published verification checks. The checks include input and archive hashes, four MPI ranks, SCF convergence, finite energy and force values, atom identity, and input/output geometry identity. The four pairs use the same PBE plane-wave settings: GPAW 26.7.0, 500 eV, Gamma-only, Fermi-Dirac smearing 0.1 eV, and the recorded 1e-5 SCF thresholds. Native extrapolated and free energies are both retained.

Frozen model: V16 checkpoint `MACE_periodic_v16_interface_energy_run-45.model`, SHA256 `918aad24341c6892d65ec0d31b19baf610ac0efde1cfd9d2e1c0567320e0cf98`. The reviewed selection record says 80 epochs, selected epoch 80 (zero-based checkpoint index 79), training/validation-only selection, and no blind-label selection. Inference used CPU, PyTorch 2.5.1 and mace-torch 0.3.16.

## Matched response results

The central secant is `(E_plus - E_minus)/(q_plus + q_minus)`. The generalized force is the sum of the forces on displaced atoms projected onto their common positive displacement axis. A small `secant + mean generalized force` residual indicates consistency over that finite interval under the stated sign convention; it does not establish an infinitesimal derivative.

| Pair | Step (Å) | DFT native secant | DFT free secant | V16 secant | DFT residual native / free | V16 residual | V16 force-response RMS / max (eV/Å) |
|---|---:|---:|---:|---:|---:|---:|---:|
| AgC framework neighbor | 0.04 | -0.5309 | -0.5297 | -0.4960 | +0.0006 / +0.0018 | -0.0079 | 0.0235 / 0.0740 (Ti#604589) |
| AgSi transverse Ag | 0.15 | -0.8436 | -0.7994 | -0.8030 | -0.0716 / -0.0273 | -0.0215 | 0.0380 / 0.0919 (C#367111) |
| AgSi framework neighbor | 0.04 | -0.7247 | -0.7558 | -0.7650 | +0.0325 / +0.0014 | +0.0022 | 0.0248 / 0.1397 (Si#366934) |
| AgTi framework neighbor | 0.04 | +2.7478 | +2.6879 | +2.6819 | +0.0583 / -0.0015 | -0.0025 | 0.0186 / 0.0591 (C#480805) |

All quantities in the slope/residual columns are in eV/Å. The force-response RMS is over persistent atom IDs in the `plus - minus` force-vector difference. The maximum atom is listed with its element and ID.

The free-energy secant has a smaller DFT force-consistency residual than the native-energy secant for the AgSi transverse, AgSi framework, and AgTi pairs. AgC is close under both energy definitions. This is evidence to preserve and compare both energy definitions; it is not evidence that the DFT forces or archives are wrong. The 0.15 Å AgSi collective displacement has the largest finite-interval ambiguity, and its native/free secants differ by about 0.044 eV/Å.

V16's model-energy secants are close to its mean projected forces for all four pairs. For AgSi framework and AgTi, its secants are closer to the DFT free-energy secants than to the native-energy secants. This points to an energy/force convention that must be made explicit when building V17. The old V14 parent records store `REF_energy` and `REF_forces`, but their `source_method` does not specify smearing or native/free energy convention; parent-relative energy offsets are therefore not interpreted as absolute errors.

## Force-vector and spatial response

The per-structure V16 force-vector RMSE on these eight targeted diagnostics ranges from 0.0328 to 0.0464 eV/Å. This is a diagnostic-set observation, not a validation pass: the structures were chosen for targeted training acquisition, and these points do not replace the frozen development-validation or withheld-test roles.

Signed marked-pair separation-force errors (`V16 - DFT`) range from -0.0332 to +0.0399 eV/Å across the eight structures. The projection uses the direct Cartesian vector from marked Ag to contact atom, matching the historical V16 evaluator. This scalar alone does not describe the full force error.

Average per-atom vector-error RMSE for the 3.5 Å MIC framework shell was 0.0592 eV/Å for AgC framework motion, 0.0486 for AgSi framework motion, 0.0432 for the AgSi transverse shift, and 0.0561 for AgTi framework motion. Moved-atom group errors averaged 0.0533, 0.0699, 0.0249, and 0.0423 eV/Å, respectively. The force response is not only a marked-pair distance issue:

- AgC's largest per-structure residual is Ti#604589 in the framework-neighbor shell (0.1076 eV/Å on the positive displacement); the paired response error also peaks there (0.0740 eV/Å).
- AgSi framework motion has a 0.1397 eV/Å maximum paired response error on moved Si#366934. Ti#368827, outside the defined 3.5 Å shell, has per-structure errors up to 0.1100 eV/Å.
- For the collective AgSi transverse shift, the largest paired response error is C#367111 in the framework shell (0.0919 eV/Å); Ti#368827 outside that shell reaches 0.1047 eV/Å per structure.
- AgTi framework motion has its largest per-structure residual on shell Ti#479021 (0.0961 eV/Å); the largest paired response error is C#480805 in the shell (0.0591 eV/Å).

These local and outside-shell residuals support adding both signs of the targeted framework and transverse acquisitions to the V17 training pool. They do not establish that changing a loss weight or optimizer parameter would fix the remaining V16 screening failures. A controlled ablation would be needed to attribute those failures to training parameters.

## Provenance and energy-convention findings

The V14 training dataset SHA256 is `3323cdbd3081fc8d91673662c5f83c8da93ba3b4d7a23428f6b70d0e7a4c9faf`; the relevant parent rows are AgTi frame 28, AgC frame 29, and AgSi frame 30. Dataset hash, frame index, configuration type, formula, atom count, control IDs, symbols by persistent atom ID, cell, and PBC agree with the acquisition manifests. Each child is also confirmed to be the declared signed displacement from that exact row within 3e-8 Å using MIC, with no undeclared moved atoms.

There is one metadata issue: the upstream `frame_geometry_sha256` does not match a fresh hash of the referenced row using the canonical function in the V17 candidate generator. The computed and declared values are:

| Parent | Declared upstream frame hash | Recomputed row hash |
|---|---|---|
| AgC, frame 29 | `72d649542158aade39a3cf50dc73759ad25df8a7ff7b2b117f1aa0c6ffc513a3` | `76fdf4d13972b2ebed7b0476f5166607f003ffa690b36d53707afec196aaccf8` |
| AgSi, frame 30 | `4a5540c7fd824a01df8fcb6155e62c7ff045ec36d637ac5f405247c9b5dbe4e5` | `15d350743e243e83f26abb27101e983b8970bb63a1569d56b6b65d56fcd43226` |
| AgTi, frame 28 | `18dc1be733a506b748386175293d2f094e720e3e1c9dc96e259623088c005b7c` | `0ad9bc5d0d6636b02c31048d0c33aa198786e5d4fa784a668701c1e744bcd764` |

The input and archive hashes still verify, and the row-level identity and displacement checks above pass. The hash discrepancy is therefore a source-metadata problem rather than a failed DFT archive check, but it must be reconciled in the V17 integration/provenance record rather than described as an exact geometry-hash match.

## Recommendation to A

1. Add all four complete matched pairs (all eight verified labels) as V17 training acquisitions after recording the parent-hash discrepancy and retaining the exact V14 dataset SHA, frame index, IDs, and direct displacement evidence. No pair is excluded by archive or convergence checks.
2. Keep native extrapolated and free energies as separate fields. Before training/evaluation, state which energy definition is used with the finite-smearing forces; do not compare or combine native and free energies in a single metric. The three non-AgC pairs show a clearer finite-difference alignment of forces with the free-energy secant.
3. Keep both signs of each pair together. Prioritize checking the indicated local C/Ti/Si response atoms in V17 common-geometry evaluation. Do not change loss weights or other training parameters solely from these four pairs; the evidence supports coverage and convention checks but does not isolate a parameter cause.
4. Preserve A's frozen role assignment: use only approved development-validation labels for selection, and keep withheld-test labels sealed until the V17 model and checkpoint are frozen.

## Reproduction and files

The analysis script and exact command are documented in `README.md`. `analysis_final/analysis.json` contains the model/selection hashes, source-manifest hashes, archive hashes, methods, parent identity checks, and force/energy conventions. CSVs provide pair metrics, per-structure errors, response groups, persistent-ID force vectors, and paired directional responses. `SHA256SUMS.txt` covers the published report, script, README, JSON and CSV artifacts.
