# P1P artifact and checkpoint index

The committed result set contains the reports, summary and time-series CSV files, validation scripts, logs, and phase-colored figures needed to review the P1P interface/replay analysis and A/B/C comparison.

The raw restart archives are intentionally not included in this commit: they add about 2.5 GB and are not required to inspect the reported results. They remain in the current workspace and were not deleted:

- Full-run restart fields and checkpoints: `abc_full/A/`, `abc_full/B/`, `abc_full/C/`
- Short-run restart fields and checkpoints: `abc_short/A/`, `abc_short/B/`, `abc_short/C/`
- Instrumented C replay input checkpoint: sibling `../p1p_cell_failure_diag_20261107/restart.dump` and its sidecars

Each run folder also contains the CSV histories, text logs, run settings, and checkpoint hashes selected for commit. The diagnostic replay source and per-cell CSV are committed; the large restart payload is kept locally for recovery if a later investigation needs it. The P1P production solver and Ag70 geometry are unchanged.
