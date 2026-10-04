# MACE interface model versions

This note summarizes the changes supported by the files preserved in this repository. A model version name alone does not prove which single training variable changed.

| Version | Preserved change | What it was intended to address |
|---|---|---|
| v9 `forceonly` | Initial four-element interface potential, trained from force-focused data. The model is preserved; its complete original training command and dataset are not in this handoff package. | Establish a first force-capable Ag/Ti/Si/C baseline. Absolute energy calibration was not established. |
| v10 `energycal` | Added PW-PBE reference cells for fcc Ag, hcp Ti, diamond Si, and diamond C; added the Ag–Ti 2.55 Å periodic energy/force point and an energy label at the 2.3871 Å frame. Kept Ag–Ti 2.70 Å out as a same-motif distance holdout. Energy-composition rank is 4/4. | Make four-element total-energy offsets identifiable and calibrate Ag–Ti periodic energies. |
| v11 `energyfocus` | A v11 checkpoint is preserved. Its data directory contains a manifest and README identifying the v10 energy-calibration dataset and 23/2/3 split; no separate v11 training command/log is preserved, so the exact v10-to-v11 parameter change cannot be reconstructed from this package. | Improve the energy-focused model after v10; the precise isolated experimental factor is undocumented. |
| v12 `interface_energy` | Preserved the v11 split, replacing four Ag–Si/Ag–C force-only LCAO frames with same-geometry PW-PBE energy/force labels: AgSi 2.2/2.8 Å and AgC 2.0/2.4 Å. Train/valid/test remain 23/2/3; train energy labels increase from 15 to 19. Trained 80 epochs, selecting epoch 76. | Replace inconsistent or incomplete interface labels with periodic plane-wave PBE energies and forces. |
| v13 `interface_energy` | Preserved every v12 training frame and added three verified PW-PBE distance points: AgC 2.30 Å, AgSi 2.60 Å, AgTi 2.60 Å. Frozen valid/test files are byte-identical; split is 26/2/3 with energy-composition rank 4/4. Trained 80 epochs, selecting epoch 72. | Add distance-response information without replacing v12 results. |

## Comparable metrics

The current evaluator scores all saved models on the same three frozen structures, two external midpoint structures, and (for v9–v12) three additional distance points. The midpoint geometries used below are related to the v13 distance-scan training structures, so they are not independent registry or thermal validation.

| Model | Frozen-test force-vector RMSE (eV/Å) | External midpoint force-vector RMSE (eV/Å) | External midpoint energy MAE (meV/atom) |
|---|---:|---:|---:|
| v9 | 0.07339 | 0.09239 | 312.60 |
| v10 | 0.06287 | 0.10076 | 257.35 |
| v11 | 0.06323 | 0.09033 | 155.70 |
| v12 | 0.05821 | 0.09224 | 11.69 |
| v13 | 0.05716 | 0.02412 | 3.18 |

On the Ag–C midpoint alone, force-vector RMSE is 0.10418, 0.11390, 0.10512, 0.11659, and 0.02546 eV/Å for v9 through v13. v13 also reduces the separating-force error relative to v12, but it remains 0.110 eV/Å on that one geometry. The v13 midpoint check shares its source motif with a distance point added to v13 training.

The current validation status remains **UNDETERMINED**. Keep long MD and TTM on hold until v13 is checked on genuinely different Ag–C and Ag–Si interface registries and thermal-like geometries, with agreed acceptance tolerances. The per-structure data and hashes are in `mace_periodic_v13_interface_energy/results/v9_v10_v11_v12_v13_comparison.json`.
