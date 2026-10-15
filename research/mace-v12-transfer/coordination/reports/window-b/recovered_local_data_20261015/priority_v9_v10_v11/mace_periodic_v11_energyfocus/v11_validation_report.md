# MACE v11 validation report

## Decision

v11 passes the current **small-cell energy/force smoke gate** and the two-rank LAMMPS/ASE parity check. It is suitable for further controlled validation. It is **not** yet established for long MD, TTM, or production use: independent AgSi and AgC total-energy offsets remain large, and each contact currently has only one small external PW-PBE structure.

## Training

- 80 epochs completed; the trainer exported its selected checkpoint from epoch 62.
- v11 uses the same 23/2/3 train/valid/test structures as v10, with identical copied split files and unchanged external AgSi/AgC holdouts.
- The only objective change from v10 was energy weight 100 instead of 1; force weight remained 1000.
- Best-checkpoint validation errors: 10.9 meV/atom energy RMSE and 40.3 meV/Å force RMSE.
- Model: `checkpoints/MACE_periodic_v11_energyfocus_run-44.model`.

## Frozen test comparison

The same 3 DFT structures (89 atoms total) were evaluated with v9, v10, and v11.

| Metric | v9 | v10 | v11 |
|---|---:|---:|---:|
| Energy MAE | 41.8 meV/atom | 83.4 meV/atom | 10.6 meV/atom |
| Maximum absolute energy error | 75.5 meV/atom | 129.5 meV/atom | 16.8 meV/atom |
| Force vector RMSE | 0.0734 eV/Å | 0.0629 eV/Å | 0.0632 eV/Å |

v11's AgTi 2.70 Å energy error is 14.3 meV/atom; C-gap is 16.8 and Si-gap is 0.9 meV/atom.

## Independent external structures

- **AgSi, 51 atoms:** v11 absolute energy error 107.9 meV/atom versus v9 129.7; full force-vector RMSE 0.0737 versus 0.0797 eV/Å. Central separating force is 1.815 versus DFT 1.708 eV/Å.
- **AgC, 48 atoms:** v11 absolute energy error 203.5 meV/atom versus v9 495.5; full force-vector RMSE 0.1051 versus 0.1042 eV/Å. Central separating force is −1.424 versus DFT −1.415 eV/Å.

The contact forces are close on these two structures, but their energy offsets remain too large to claim calibrated AgSi/AgC interface energies. AgC is a method-transfer check on a previously studied geometry, not a new motif.

## Two-rank LAMMPS smoke

- CPU-only Kokkos/ML-IAP, 2 MPI ranks, 49 atoms, periodic AgTi validation cell.
- `run 0`: LAMMPS PE −177.69739 eV; ASE with the same v11 model −177.69738763 eV; absolute energy difference 2.37×10⁻⁶ eV.
- Force-vector RMSE against ASE: 2.46×10⁻⁶ eV/Å; maximum component difference 4.82×10⁻⁶ eV/Å.
- One unheated NVE step: final step 1, atom count 49 throughout, no lost atoms or nonfinite values. Printed total energy is unchanged (−177.69739 eV); at step 1 KE is 0.0004383 eV, PE is −177.69783 eV, temperature is 0.0706 K, and Fmax is 4.36007 eV/Å.
- The `PackForwardComm has only been tested on Kokkos devices` warning remains. The two-rank CPU run nevertheless matched ASE to dump precision.
- This smoke did not invoke TTM or heat the system. It is a LAMMPS interface check, not a physical trajectory validation.

## Evidence files

- v9/v10/v11 frozen comparison: `../mace_periodic_v10_energycal/v9_v10_v11_frozen_test_comparison.json` and `.csv`.
- v9/v10/v11 external comparison: `../mace_periodic_v10_energycal/v9_v10_v11_external_holdouts.json` and `.csv`.
- LAMMPS/ASE parity: `lammps_validation/run0_ase_parity.retry4.json`.
- One-step NVE summary: `lammps_validation/run1_nve_summary.json`.
- Full LAMMPS logs and dumps are in `lammps_validation/`.

## Next step

Add same-method PW-PBE energy-and-force labels for AgSi and AgC interface environments while preserving new independent holdouts, retrain v12, and repeat the frozen/external comparison. Only after those checks should an unheated stability segment or TTM plumbing be considered. Do not start long MD or production TTM from this result.

`FROZEN_TEST_GATE = PASS`
`LAMMPS_MLIAP_PARITY = PASS`
`AGSI_AGC_ENERGY_VALIDATION = INCOMPLETE`
`CROSS_POTENTIAL_VALIDATION = NOT_ESTABLISHED`
`NEXT_STEP = ADD_INDEPENDENT_PW_PBE_AGSI_AGC_DATA_AND_TRAIN_V12`
