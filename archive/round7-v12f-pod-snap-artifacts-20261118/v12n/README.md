# Round7 v12n POD/SNAP C–C cold-hotspot candidate

This archive contains the complete v12n POD and SNAP fits, exact 128/12 data split, local DFT force-label crops, fitting logs, Si01 holdout evaluation, the 458,248-atom input geometry, full-model static screens, and 1 fs no-source cold trajectories. Compressed LAMMPS state/dump files are included for reproducibility.

## Fit and independent checks

- Training: 128 baseline configurations plus 24 local frames (152 total). Validation: the unchanged 12-frame set. Si01 center atom 104359 remains excluded from training.
- Added two GPAW PBE/PW350 Gamma labels for C128985–C134363 at separations 1.386 and 1.114 Å (four central C force targets total).
- Validation RMSE: POD energy 0.076566 eV/atom, force 0.594333 eV/Å; SNAP energy 0.061876 eV/atom, force 0.718041 eV/Å. Both are within the exploratory aggregate targets (~0.1 eV/atom and ~3 eV/Å).
- Si01 held-out center vector error: POD 7.6449 eV/Å; SNAP 36.7102 eV/Å. Both miss the ~3 eV/Å single-point target.
- Full-model static Fmax/F99: POD 17.220/5.479 eV/Å; SNAP 102.321/19.198 eV/Å. Neither exceeds the emergency 500 eV/Å stop line. This is not a stability pass.

## Cold stability outcome

All runs use dt=1 fs, no laser source, fixed skeleton atoms in z=0–20 Å, and a 300 K Langevin bath in z=20–40 Å; all other non-fixed atoms use NVE.

- POD 10 fs: no atom loss; temperature 300→710 K, maximum displacement 0.371 Å.
- POD 50 fs: the originally targeted C–C pair ends at 2.553 Å; the v12m collapse at 31 fs is not reproduced. However, the full model is not cold-stable: temperature reaches 3381 K, maximum displacement 4.090 Å, 12,622 Ti–C pairs are shorter than 1.8 Å, minimum Ti–C is 1.248 Å, and minimum C–C is 1.014 Å. Atom count is preserved and final Fmax is 87.76 eV/Å, but the severe structural rearrangement fails the cold-stability gate.
- SNAP 10 fs: no atom loss, but temperature reaches 3523 K, maximum displacement 1.046 Å, 2,143 Ti–C pairs are shorter than 1.8 Å, and minimum Ti–C is 1.441 Å. It fails the cold-stability gate.

**Decision:** v12n remains an exploratory training candidate. Do not use it for heating/cooling or solidification until the unresolved Si01 error and full-model spontaneous structural collapse are addressed. The 500 eV/Å number is an emergency force stop, not a model accuracy or cold-stability criterion.
