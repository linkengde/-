# v15 local-registry blind holdout inputs (Window A)

These three fixed-geometry inputs are frozen from Window B’s geometry-only audit. The Ag-C proposal is the prior unlabeled geometry reassessed with species-aware C-Ti framework context. No energy or force labels are included.

They test local lateral Ag registry changes over existing small-cluster framework motifs. They are **not** independent morphology or extended-interface tests. Keep all three out of v15 training, validation, checkpoint selection and hyperparameter tuning; score only after freezing the model.

`input_manifest.json` records source hashes, candidate hashes, IDs, cell/PBC, central pairs, pair minima and split role. This directory currently contains no DFT output. Planned reference: GPAW 26.7.0, PW-PBE, 500 eV, Gamma, 0.1 eV Fermi smearing, four MPI ranks.
