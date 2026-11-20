# Round7 POD/SNAP v12p Ag-contact candidate

## Outcome

v12p is an exploratory candidate trained with four new, converged DFT contact crops from a no-laser 10 fs SNAP v12o cold-screen endpoint. The fit and full-model short checks are complete. **Neither POD nor SNAP is cleared for heating** because both show strong no-laser self-heating; SNAP also collapses many Ti-C contacts.

## Data and fit

- The original split remains 128 training / 12 validation configurations.
- Si01 center atom 104359 remains excluded from training and is a development holdout.
- Four local crops (two Ag-Si, one Ag-C, one Ag-Ti) add 8 central DFT force rows. The local set has 30 frames and 158 total training configurations. Crop-edge forces use POD v10 anchors; local energy weights are zero.
- The 4 GPAW PBE/PW350 Gamma labels converged in 52, 40, 37, and 42 SCF iterations, respectively. All labels, logs, JSON records, extraction scripts, and contact provenance are included.

| Model | 12-frame validation energy RMSE | 12-frame validation force RMSE | Si01 development holdout vector error |
|---|---:|---:|---:|
| POD v12p | 0.078328 eV/atom | 0.574072 eV/A | 6.554 eV/A |
| SNAP v12p | 0.058697 eV/atom | 0.670792 eV/A | 15.284 eV/A |

The ~3 eV/A force goal is exploratory. Both Si01 errors remain above it. The 12-frame validation set is small and does not establish transferability. SNAP's weighted local-energy diagnostic is NaN because those configurations intentionally have zero energy weight; local-force and validation metrics are finite.

## Full 458,248-atom static screen

| Model | Fmax | F99 | Atoms above 500 eV/A | Nonfinite force rows |
|---|---:|---:|---:|---:|
| POD v12p | 14.923 eV/A | 5.095 eV/A | 0 | 0 |
| SNAP v12p | 48.295 eV/A | 13.894 eV/A | 0 | 0 |

The 500 eV/A limit is only an emergency extreme-force screen. Passing it does not establish MD stability.

## No-laser cold stability screens

Both use dt=1 fs, 10 steps (10 fs), a fixed bottom skeleton at z=0-20 A, a 300 K Langevin bath at z=20-40 A, and NVE for remaining mobile atoms.

### POD v12p

- Atoms retained: 458,248/458,248; no atom loss.
- Mobile temperature: 300 to 662.2 K without a heat source.
- Potential energy change: -19.958 keV; kinetic energy change: +19.777 keV; total-energy change: -180.74 eV.
- Maximum displacement: 0.338 A; 1 Ti-C pair below 1.8 A; minimum Ti-C distance 1.783 A.

### SNAP v12p

- Atoms retained: 458,248/458,248; no atom loss.
- Mobile temperature: 300 to 2951.0 K without a heat source.
- Potential energy change: -146.154 keV; kinetic energy change: +144.744 keV; total-energy change: -1410.13 eV.
- Maximum displacement: 0.939 A; 1701 Ti-C pairs below 1.8 A; minimum Ti-C distance 1.448 A.

POD cold-start relaxation is still substantial; SNAP is severely unstable. Do not use either candidate for heating until the initial-structure relaxation and subsequent no-laser cold stability are resolved. No heating was started in this work.

## Archive map and reproduction

- `candidate/v12p_10fs_ag_contacts_candidate/`: training and validation frames, local labels, POD/SNAP inputs and fitted coefficients, logs, Si01 evaluation.
- `active_labels/`: DFT crop provenance and raw outputs.
- `full_model/`: static and 10 fs cold-screen inputs, logs, audits, and gzip-compressed force/end-state dumps/data.
- `scripts/`: fitting, runner, audit, and continuation scripts used.
- Shared starting structure is archived at `../v12n/full_model/Ag_Ti3SiC2_bicontinuous_200x200x300.data.gz`.

The archived input files retain cloud-workspace absolute paths. For a new environment, update paths to the local LAMMPS executable, FitSNAP/POD files, candidate directory, and shared starting model. Decompress `.gz` data/dumps before direct LAMMPS or OVITO loading.
