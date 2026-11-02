# V19 train55 Ti-response experiment

This experiment appends six verified v23 topmost-Ti displacement labels to the byte-identical train49 candidate. The original three-frame dev3 is unchanged. This is an experiment-only fit: official train36/dev3 were not modified, and test or sealed data was not opened.

## Candidate data and provenance

The six fixed-geometry labels cover AgTi, AgSi and AgC registry-B slabs. The four symmetry-equivalent topmost Ti atoms were displaced together by −0.02 Å and +0.02 Å normal to the surface. The DFT method is GPAW 26.7.0 / ASE 3.29.0, PBE, 500 eV, 6×6×1 k points, 0.10 eV Fermi width, xy dipole correction and four MPI ranks. The result archives were verified and uploaded in source commit `94801a2ad963bdd32fcfcd871104953f0ca14228`.

The candidate has 55 training frames and the same three frozen development frames. Dataset construction checks source archive and input hashes, energy/force mapping, exact duplicates, and the byte-identical train49 prefix. Train SHA256: `5a9eabc8c31d83670a81f73ba8da1ff18d3fe2311e6bec062aca3be1c8b094aa`. Dev3 SHA256: `e9e86e3a353c50c8173cefede370e2653ec2e81506f4b8122b01b7a8efe76e4c`.

## Fixed fit and verification

The fit repeats the established recipe: seed 45, CPU, four threads, batch size 1, 80 epochs, learning rate 1e−4, weight decay 5e−7, energy weight 100, force weight 1000, and `REF_energy` / `REF_forces` targets. The fixed epoch-79 model was screened against frozen dev3. The verifier confirmed all 80 epoch records and checkpoint files before pruning, and confirmed tensor identity between the epoch-79 checkpoint and exported model. Because the 55-frame training set includes both Ti–Si–C substrate configurations and Ag-containing interface slabs, this is an Ag–Ti–Si–C interface cross-potential candidate, not a Ti–Si–C-only skeleton fit. Its Ag interaction coverage is limited to the sampled interface configurations; it is not a general Ag alloy potential.

| Interface | Train55 force-vector RMSE (eV/Å) | Historical ≤0.05 gate | Current ≤0.07 gate | Separation-force error (eV/Å) | ≤0.10 gate | Max force error (eV/Å) |
| --- | ---: | --- | --- | ---: | --- | ---: |
| AgC | 0.05725 | Fail | Pass | 0.02158 | Pass | 0.10548 |
| AgSi | 0.05827 | Fail | Pass | 0.08379 | Pass | 0.08373 |
| AgTi | 0.06579 | Fail | Pass | 0.03753 | Pass | 0.14411 |

On 2026-11-02, the working per-interface force-vector RMSE tolerance was set to ≤0.07 eV/Å for the intended practical use; under this criterion, all three interfaces pass, and all separation-force errors remain ≤0.10 eV/Å. Aggregate force-vector RMSE is 0.06056 eV/Å. The original ≤0.05 eV/Å assessment remains a historical fail in `run/results/train55_verification.json`; it has not been overwritten. This is a post-hoc reclassification, and frozen dev3 was monitored as validation data during fitting, so it is not an independent blind validation. The result passes the current working precision screen only; the model remains `train_acquisition_candidate_pending_audit` and is not promoted to the official training set or production use. The force criterion is RMSE; it does not bound every atom's maximum force error (see the table).
The frozen-dev3 residual diagnostic shows improvements over train49 for Ti in all three interfaces and for C in AgC/AgSi. Remaining force error is distributed across species and interfaces; AgC Si residual increased, so the six topmost-Ti labels did not address every error mode. These are diagnostic correlations, not a causal attribution. The dev3 geometries and labels were not copied into training.

## Next acquisition

Keep train55 as the already-fitted experiment-only cross-potential baseline. The v24 acquisition is currently collecting six additional Ag–Ti–Si–C interface carbon-response labels. After all six archives pass verification, build the next candidate from those labels and rerun the fixed recipe against frozen dev3 to check robustness under the current 0.07/0.10 working screen. Do not tune epochs or hyperparameters against dev3, and do not alter the official split or inspect sealed test data.

## Artifacts

- `dataset/`: train55, frozen dev3 and provenance manifest with checksums.
- `run/results/train55_verification.json`: original fixed-epoch evaluation using the historical 0.05 eV/Å gate, plus checkpoint/model identity audit.
- `run/results/train55_revised_gate_assessment_20261102.json`: explicit reclassification under the current 0.07 eV/Å working tolerance; records scope and limitations.
- `run/results/train55_residual_diagnostics.json`: per-interface, per-species and distance-resolved dev3 error comparison against train43/train49.
- `run/checkpoints/`: only the latest epoch-79 checkpoint is retained after full-coverage verification.
- `run/models/`: exported and compiled epoch-79 models; `run/logs/` and `run/train.stdout` preserve training evidence.

To rebuild the dataset and verify the training/evaluation artifacts, use `prepare_train55.py`, `preflight_train55.py`, `run_train55.sh`, `verify_train55.py` and `diagnose_train55_residuals.py`. The fit script restarts from the latest saved checkpoint if interrupted before completion.
