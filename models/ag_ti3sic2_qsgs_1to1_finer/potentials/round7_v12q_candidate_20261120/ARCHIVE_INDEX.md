# Round7 Ag–Ti3SiC2 four-element v12q potential archive

This is the fitted exploratory v12q candidate pair used for the current large-model TTM/MD attempt. POD is the active model; SNAP is retained as an independent comparison. The directory includes both coefficient sets, model descriptors/parameters, fit inputs and logs, 128 training configurations, 12 validation configurations, local training data, the Si01 development-holdout inputs/results, and the original v12q README.

## LAMMPS files

- POD: `POD/AgTiSiC_POD_v12q.pod` plus `POD/pod_v12q_coefficients.pod`; fit recipe and logs are beside them.
- SNAP: `SNAP/AgTiSiC_SNAP_v12q.mod`, `.snapcoeff`, and `.snapparam`; FitSNAP input, training/validation data, and metrics are in `SNAP/`.
- Element order for both files: Ag, C, Si, Ti.

## Recorded validation

- POD: 12-frame validation RMSE, energy 0.080395 eV/atom and force 0.568209 eV/Å.
- SNAP: 12-frame validation RMSE, energy 0.063991 eV/atom and force 0.710311 eV/Å.
- Si01 development holdout vector-force error: POD 6.594 eV/Å; SNAP 12.017 eV/Å.

These are exploratory fitted candidates. The recorded errors and holdout do not establish high-temperature transferability or production readiness. The active heating run remains separate from this archive and is not included here.

## Refit note

`POD/fitpod.data` and `run_v12q_fits.sh` preserve the original cloud paths and executable locations. The fitted potential files can be used directly; to refit elsewhere, update those paths for the checkout and installed FitPOD/FitSNAP tools. `SHA256SUMS.txt` records hashes for all archived files except itself.
