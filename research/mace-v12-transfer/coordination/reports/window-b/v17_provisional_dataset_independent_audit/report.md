# Independent provisional V17 dataset audit

43 training frames verified:35 inherited labels exactly preserved and8 signed diagnostic labels exactly copied from source output. Exact composition rank4;29 directly supported and14 partial-source rows. Newly computed frame0 replacement is not integrated: source frame0 remains unchanged. Finite training energies/forces, ordered species, coordinates, cell/PBC and manifest provenance counts pass. Three development geometries were parsed without label values; no exact/near train overlap under published thresholds. Frozen withheld geometry hashes do not match training; withheld inputs and outputs were not opened.

Source snapshot changes: []

Any stale partial aggregate manifest reflects subsequent B aggregate reconciliation rather than a changed label. Per-source details in audit.json retain expected/current hashes; other unexplained changes require A review. Do not silently rewrite the production manifest. Original inherited method/convention uncertainties remain: this pool is PROVISIONAL, not fully PW-PBE certified. Training-source checks and full rank do not certify model force quality. Builder/preflight were inspected; preflight was not run because it reads development reference values. No production dataset changes or model/DFT calculations during audit.
