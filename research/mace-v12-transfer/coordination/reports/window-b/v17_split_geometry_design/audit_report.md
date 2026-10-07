# V17 split geometry design — audit summary

The first proposed V16-derived validation geometries were rejected: each exactly matched the corresponding V16 blind-holdout input. That overlap was discovered by extending the comparison to all frozen probe inputs. The preserved geometry-only evidence is under `attempt_diagnostics/attempt_01/`.

The revised set in `final_set/` contains one development-validation geometry and two withheld-test geometries for each interface. The screen covered all 197 frames across 18 historical model-split files (including all 35 V16 training rows) and 38 published input geometries, including all eight hash-verified planned V17 training inputs. Each final candidate has zero exact/near hits against that inventory, zero pairwise role-set overlaps, and passed its species-specific minimum-distance guard. Closest-geometry and element-pair details are in `final_set/screening_manifest.json`.

Candidate structures preserve atom IDs/order, atom count, cell and PBC and contain no energy/force labels. No DFT output, blind label, model inference/training, MD or TTM was read or run. The candidates remain correlated with existing small-cluster parent motifs and do not establish morphology independence or thermal coverage. A must review the role manifest before training; withheld-test labels must remain sealed from selection.
