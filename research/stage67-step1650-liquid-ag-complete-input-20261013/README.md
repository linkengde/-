# Stage67 step1650 liquid-like Ag post-processing input

## Source frame

- Source: Stage67, 200 Å transverse stable-meltpool reference branch (Stage67_200A_4x_meltpool_to_0p300_v4_repaired).
- Frame: step 1650, t = 0.165 ps.
- Cell: 150.1 × 200 × 200 Å; 303,372 atoms; periodic flags pp pp pp.
- Intended use: cloud-side post-processing and extraction of liquid-like Ag.
- No LAMMPS run or MD restart was performed while preparing this package. Source files were copied without modification.

## Available inputs

- input/high_frequency_step_1650.dump: source atom frame; includes CellTl and CellTl_valid columns.
- input/q6_hf_step_1650.json: existing step1650 Q6/PTM-and-cutoff-cluster summary.
- input/main_cluster_ids_step_1650.txt: IDs for the existing maximum cluster.
- input/thermal_history_stable_window.csv: existing 11-sample thermal/energy summary for 0.155–0.165 ps; not a per-atom temperature trajectory.
- input/reference_manifest.json and input/SHA256SUMS_source.csv: source provenance and source-package hashes.
- reports/: existing stability, lineage, manual review, input-audit, and extraction reports.

## Frozen analysis provenance and limitations

The existing implementation's candidate mask is Ag (type=4) with PTM Structure Type == 0, CellTl_valid == 1, and CellTl >= 975 K; the existing cluster cutoff is 3.74625 Å. The Q6 threshold 0.446766 is recorded in the summary metadata; the available implementation uses the PTM Other proxy for the per-atom candidate mask. Do not silently change this criterion.

The existing summary reports 8,728 liquid-like Ag atoms and a 1,032-atom maximum cluster. The available package does not contain per-atom PTM/Q6 labels or a LiquidCandidate table for all 8,728 atoms. There is no separate CellTl.csv; the step1650 dump itself contains CellTl and CellTl_valid. Thus this package gathers all matching source artifacts found, but reproducing the full 8,728-atom selection requires cloud-side analysis with the frozen method and compatible tools. Do not treat the 1,032-atom maximum cluster as the full liquid-like population.

## Layout

- input/: source frame and existing machine-readable inputs.
- reports/: existing reports and supporting audit tables.
- SHA256SUMS.txt: SHA256 for every packaged file except this checksum file itself.