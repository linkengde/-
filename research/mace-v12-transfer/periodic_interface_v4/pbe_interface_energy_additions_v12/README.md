# AgSi / AgC PW-PBE additions for v12

This run adds energy-and-force labels to existing AgSi and secondary AgC scan geometries. Each input is recalculated with GPAW 26.7.0, PBE, PW(500 eV), Gamma, 0.1 eV Fermi smearing, and 1e-5 energy/density/eigenstate convergence. The existing external midpoint structures and AgTi 2.70 A frozen test remain excluded.

The four selected geometries replace the matching LCAO-only force frames in the v12 training split, rather than being duplicated. Compact extxyz labels, summaries, progress records, and logs are archived; large rolling GPAW states remain under `/tmp`.

The independent distance-point inputs in `inputs/validation_distance_scans/` are separate from training. They are central-pair displacements of the AgC, AgSi, and AgTi midpoint geometries, with all other atoms and the cell fixed. Their input hashes and target distances are in `validation_distance_input_manifest.json`. Run them with four MPI ranks using `run_validation_distance_scans.sh`; completed compact outputs are checked by `archive_validation_pw.py`.
