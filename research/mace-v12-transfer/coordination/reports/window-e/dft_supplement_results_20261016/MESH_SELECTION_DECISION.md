# Mesh selection for the six Ti3SiC2 diagnostic points

The first 10x10x4 attempt was stopped before a valid label after about 41.75 minutes and two partial SCF iterations. The chosen 6x6x4 mesh is based on existing converged DFT controls for the same 48-atom parent geometry and the same PBE/PW/500 eV/0.10 eV method.

Direct comparison of the archived parent output at 6x6x4 and 10x10x4 confirms identical cell, PBC, ordered atoms and positions. The native-energy difference (10x10x4 minus 6x6x4) is -0.1416505 meV/atom. The force-vector RMS difference is 0.0062581 eV/Angstrom, and the maximum per-atom force-vector difference is 0.0097241 eV/Angstrom. This is well below the existing Ti3SiC2 model-force errors (order 2 eV/Angstrom) and supports 6x6x4 for this diagnostic bulk family.

All six new points will use exactly 6x6x4, so their relative labels remain internally consistent and match the existing 6x6x4 parent reference. The parent-only mesh test does not prove identical error for strained or displaced cells; this is retained as a limitation. The structures remain `DIAGNOSTIC_ONLY`, and no TRAIN/DEV/TEST role is assigned.
