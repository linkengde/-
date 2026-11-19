# Round7 POD/SNAP v12o candidate (2026-11-19)

## Status

v12o was fitted from the v12n dataset plus two converged local DFT hotspot configurations. The candidate **does not pass full-model cold stability** and is **not cleared for heating**. See the 10 fs 1 fs-step reports below. This is an exploratory candidate and not a validated production potential.

## Data and fit

- Baseline split retained: 128 fitting configurations / 12 validation configurations.
- The separate Si01 center (source atom 104359) remains excluded from training and was evaluated as a development holdout.
- Added 26 local force configurations; four central atom force rows carry new DFT labels (two for the Ti-C contact, two for the C-C contact). Crop-edge force rows use POD v10 anchors. The local energy weight is 0.
- Full fit data, model inputs/coefficients, logs, DFT structures/labels, training manifest, holdout evaluation inputs/outputs are in `candidate/v12o_50fs_hotspots_candidate/`.

## Fit and independent checks

| Candidate | 12-frame validation energy RMSE | 12-frame validation force RMSE | Si01 center vector error |
|---|---:|---:|---:|
| POD v12o | 0.076972 eV/atom | 0.567965 eV/Å | 7.921826 eV/Å |
| SNAP v12o | 0.065099 eV/atom | 0.686343 eV/Å | 25.482571 eV/Å |

The Si01 value is a development holdout, not a blind test. Both exceed the exploratory ~3 eV/Å local-error goal. The validation split is only 12 configurations and does not establish transferability by itself.

## Full 458,248-atom static screen

Model: 200 x 200 x 300 Å, Ag/Ti3SiC2 approximately 1:1 by mass; bottom skeleton atoms in z=0-20 Å fixed for the force screen.

| Candidate | Fmax | F99 | >500 eV/Å | Nonfinite forces |
|---|---:|---:|---:|---:|
| POD v12o | 18.015 eV/Å | 5.437 eV/Å | 0 | 0 |
| SNAP v12o | 102.294 eV/Å | 21.397 eV/Å | 0 | 0 |

500 eV/Å is only the emergency extreme-force screen; passing it does not establish accuracy or MD stability.

## POD full-model cold screen

No laser; dt=1 fs; 10 steps=10 fs; fixed bottom skeleton z=0-20 Å; 20-40 Å 300 K Langevin bath; remaining mobile atoms NVE.

- 458,248 atoms retained; zero atom loss.
- Mobile temperature: 300.0 -> 702.0 K.
- Potential energy change: -22.151 keV; kinetic energy change: +21.950 keV; total-energy change: -201.48 eV for the full system.
- Maximum displacement: 0.343 Å; one Ti-C pair enters below 1.8 Å (minimum 1.783 Å).
- Interpretation: strong initial structural energy release occurs without a laser. The model/initial structure is not cold-stable enough for heating qualification.

## SNAP full-model cold screen

Completed with the same settings; 458,248 atoms retained, zero atom loss.

- Mobile temperature: 300.0 -> 4,172.2 K in 10 fs.
- Potential energy change: -213.527 keV; kinetic energy change: +211.422 keV; total-energy change: -2,104.91 eV for the full system.
- Maximum displacement: 1.346 Å; 2,922 Ti-C pairs below 1.8 Å (minimum 1.359 Å).
- LAMMPS loop time: 674 s for 10 steps; total command wall time: 801 s on 4 MPI ranks.
- Interpretation: severe no-laser structural relaxation; SNAP v12o fails the cold screen.

## Runtime comparison

For the same 458,248-atom model and 10-step cold screen on 4 MPI ranks, POD took 226 s in the LAMMPS loop (276 s command wall time), while SNAP took 674 s (801 s command wall time), about 3.0x slower here. These short-run rates are environment-specific; longer-run extrapolations are only rough estimates.

## Reproduction

The unchanged full model is archived at `../v12n/full_model/Ag_Ti3SiC2_bicontinuous_200x200x300.data.gz`. Unpack it in the source model directory expected by the scripts, then adjust the software paths at the top of each script to the local LAMMPS executable and candidate model directory. All training configurations and outputs required to inspect the v12o fit are included with the candidate.
