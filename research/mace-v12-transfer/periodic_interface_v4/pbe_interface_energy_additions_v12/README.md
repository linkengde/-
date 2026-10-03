# AgSi / AgC PW-PBE additions for v12

This run adds energy-and-force labels to existing AgSi and secondary AgC scan geometries. Each input is recalculated with GPAW 26.7.0, PBE, PW(500 eV), Gamma, 0.1 eV Fermi smearing, and 1e-5 energy/density/eigenstate convergence. The existing external midpoint structures and AgTi 2.70 A frozen test remain excluded.

The large rolling GPAW states are staged under `/tmp` because the project filesystem had only 4.6 GB free. Small extxyz labels, summaries, progress snapshots, and logs will be copied back after successful convergence; state files will not be archived. The four selected geometries replace the matching LCAO-only force frames in the future v12 training split, rather than being duplicated.
