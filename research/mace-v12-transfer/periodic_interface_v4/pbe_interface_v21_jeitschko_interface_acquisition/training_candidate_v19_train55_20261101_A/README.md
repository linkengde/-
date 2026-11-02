# V19 train55 Ti-response experiment

This experiment appends six verified v23 topmost-Ti displacement labels to the byte-identical train49 candidate. The original three-frame dev3 is unchanged. This is an experiment-only fit: official train36/dev3 were not modified, and test or sealed data was not opened.

## Candidate data and provenance

The six fixed-geometry labels cover AgTi, AgSi and AgC registry-B slabs. The four symmetry-equivalent topmost Ti atoms were displaced together by −0.02 Å and +0.02 Å normal to the surface. The DFT method is GPAW 26.7.0 / ASE 3.29.0, PBE, 500 eV, 6×6×1 k points, 0.10 eV Fermi width, xy dipole correction and four MPI ranks. The result archives were verified and uploaded in source commit `94801a2ad963bdd32fcfcd871104953f0ca14228`.

The candidate has 55 training frames and the same three frozen development frames. Dataset construction checks source archive and input hashes, energy/force mapping, exact duplicates, and the byte-identical train49 prefix. Train SHA256: `5a9eabc8c31d83670a81f73ba8da1ff18d3fe2311e6bec062aca3be1c8b094aa`. Dev3 SHA256: `e9e86e3a353c50c8173cefede370e2653ec2e81506f4b8122b01b7a8efe76e4c`.

## Fixed fit and verification

The fit repeats the established recipe: seed 45, CPU, four threads, batch size 1, 80 epochs, learning rate 1e−4, weight decay 5e−7, energy weight 100, force weight 1000, and `REF_energy` / `REF_forces` targets. The fixed epoch-79 model was evaluated only on frozen dev3. The verifier confirmed all 80 epoch records and checkpoint files before pruning, and confirmed tensor identity between the epoch-79 checkpoint and exported model.

| Interface | Train49 force-vector RMSE (eV/Å) | Train55 force-vector RMSE (eV/Å) | Force gate | Separation-force error (eV/Å) | Separation gate |
| --- | ---: | ---: | --- | ---: | --- |
| AgC | 0.08906 | 0.05725 | Fail | 0.02158 | Pass |
| AgSi | 0.09129 | 0.05827 | Fail | 0.08379 | Pass |
| AgTi | 0.09081 | 0.06579 | Fail | 0.03753 | Pass |

All three force-vector RMSE values remain above the ≤0.05 eV/Å gate; all separation-force errors meet the ≤0.10 eV/Å gate. Aggregate force-vector RMSE is 0.06056 eV/Å, improved from train49's 0.09039 eV/Å and close to train43's 0.06003 eV/Å. The status is `TRAIN55_EXPERIMENT_DEV3_SCREEN_FAIL`; this model is not promoted to official training or production use.

The frozen-dev3 residual diagnostic shows improvements over train49 for Ti in all three interfaces and for C in AgC/AgSi. Remaining force error is distributed across species and interfaces; AgC Si residual increased, so the six topmost-Ti labels did not address every error mode. These are diagnostic correlations, not a causal attribution. The dev3 geometries and labels were not copied into training.

## Next acquisition

Keep train55 as an experiment-only baseline. The next DFT acquisition should broaden force-response coverage on registered training-interface parent geometries, sampling subsurface Ti and C internal modes in addition to the already sampled topmost-Ti orbit. Use the same DFT method and exact source/input/hash audit, then evaluate one fixed recipe against frozen dev3. Do not tune epochs or hyperparameters against dev3, and do not alter the official split or inspect sealed test data.

## Artifacts

- `dataset/`: train55, frozen dev3 and provenance manifest with checksums.
- `run/results/train55_verification.json`: fixed-epoch evaluation, gate results and checkpoint/model identity audit.
- `run/results/train55_residual_diagnostics.json`: per-interface, per-species and distance-resolved dev3 error comparison against train43/train49.
- `run/checkpoints/`: only the latest epoch-79 checkpoint is retained after full-coverage verification.
- `run/models/`: exported and compiled epoch-79 models; `run/logs/` and `run/train.stdout` preserve training evidence.

To rebuild the dataset and verify the training/evaluation artifacts, use `prepare_train55.py`, `preflight_train55.py`, `run_train55.sh`, `verify_train55.py` and `diagnose_train55_residuals.py`. The fit script restarts from the latest saved checkpoint if interrupted before completion.
