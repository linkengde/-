# MACE v10 energy calibration report

## Decision

**v10 does not pass the energy calibration gate. Do not use it for production MD or TTM.** Force errors improved on the small frozen test, but energy errors increased relative to v9 and the independent interface checks disagree by environment.

## Training

- 80 epochs completed. The exported stage-one model was selected at epoch 76 by the trainer's combined validation loss.
- Training set: 23 structures, 15 with DFT energy labels. Composition matrix rank: 4/4.
- New PBE/500 eV reference structures: fcc Ag, hcp Ti, diamond Si, diamond C. Also included AgTi 2.3871 Å and 2.55 Å periodic PBE energies; AgTi 2.70 Å stayed in the frozen test set.
- Objective weights: energy 1, forces 1000. Best-checkpoint validation RMSE: energy 66.7 meV/atom; force 39.7 meV/Å.
- Re-evaluation over only the 15 energy-labeled training frames gives energy MAE 50.4 and RMSE 69.2 meV/atom. The MACE end-of-run table prints 2812 meV/atom for `train_Default` because it aggregates force-only frames without energy labels; that aggregate is not a valid labeled-training energy metric.

## Frozen test: v9 versus v10

Both models were evaluated on the identical 3-frame, 89-atom split.

| Metric | v9 | v10 |
|---|---:|---:|
| Energy MAE | 41.8 meV/atom | 83.4 meV/atom |
| Maximum absolute energy error | 75.5 meV/atom | 129.5 meV/atom |
| Force vector RMSE | 0.0734 eV/Å | 0.0629 eV/Å |

At the held-out AgTi 2.70 Å geometry, the absolute energy error increased from 75.5 to 129.5 meV/atom. The energy-focused split also contains Si-gap and C-gap structures; v10 improves Si-gap energy, but its C-gap energy error rises from 25.8 to 110.9 meV/atom.

## Independent interface structures

- **Ag–Si, 51 atoms:** absolute energy error increased from 129.7 to 214.5 meV/atom. Full force-vector RMSE changed from 0.0797 to 0.0866 eV/Å. The central pair separating-force error changed from 0.127 to 0.119 eV/Å.
- **Ag–C, 48 atoms:** absolute energy error improved from 495.5 to 300.2 meV/atom. Full force-vector RMSE changed from 0.1042 to 0.1139 eV/Å. The central pair separating-force error improved from 0.0623 to 0.0216 eV/Å. This is a method-transfer check on a previously studied geometry, not an unseen geometry.

These results show environment-dependent tradeoffs; they do not establish four-element cross-potential accuracy.

## Rejected post-fit diagnostic

A four-element residual correction was fitted using only the 15 energy-labeled training frames (rank 4/4). It lowered frozen-test energy MAE from 83.4 to 31.2 meV/atom, but it worsened the external Ag–C error from 300.2 to 401.2 meV/atom and did not resolve Ag–Si (166.1 meV/atom). It also degraded the pure Ag reference. The derivative model is retained only as a diagnostic and is **not** the selected model.

## Artifacts and next step

- Frozen comparison: `v9_v10_frozen_test_comparison.json` and `.csv`.
- External comparison: `v9_v10_external_holdouts.json` and `.csv`.
- Rejected residual correction: `energy_reference_correction/energy_reference_correction_metrics.json`.
- Selected raw v10 model: `checkpoints/MACE_periodic_v10_energycal_run-43.model`.

The next controlled trial is v11 on the identical splits with `energy_weight=100` and `forces_weight=1000`. Only if v11 improves the frozen and independent energy checks without unacceptable force regression should it proceed to LAMMPS static parity. No LAMMPS run, MD, or TTM was performed with v10.

`CROSS_POTENTIAL_VALIDATION = NOT_ESTABLISHED`
`NEXT_STEP = TRAIN_V11_ENERGYFOCUS_AND_EVALUATE`
