# MACE small interface comparison data

This bundle collects the existing small Ag/Ti/Si/C interface development data available in the workspace. It contains three PBE reference structures (Ag–C, Ag–Si, Ag–Ti; 48–51 atoms), two provisional MACE models (V17 and V18), and their development-evaluation JSON files plus a compact metrics table.

These are development-validation motifs used during model diagnosis. They are not fresh blind tests, are not the four D4O planar Ag/Ti3SiC2 slab candidates, and do not establish that either MACE model is qualified for production interface dynamics. The recorded evaluations report `FAIL`; use them only as diagnostic evidence.

The JSON gates are 10 meV/atom absolute energy error, 0.05 eV/Å force-vector RMSE, and 0.1 eV/Å separating-force absolute error. `metrics.csv` gives the recorded errors and per-gate comparisons for each structure.

Original source: `/workspace/stage67-300x300-model-build-worktree/research/mace-v12-transfer/periodic_interface_v4/`.
