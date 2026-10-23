# Geometry-only Ti3SiC2 `(10-11)` facet candidate

This package contains a 69-atom symmetric Ti-terminated pyramidal slab and one deterministic terminal-Ti response geometry. Both derive from COD 9009647. They form a new, unscored geometric family, not an independent crystallographic-source family.

The package is for offline source/family planning only. It contains no DFT labels and assigns no TRAIN, DEV, or TEST role. The pair may be considered for same-source facet-held-out development diagnostics after its calculation method and resources are reviewed. It must not be described as a source-independent blind validation set. TEST structure selection and labels remain outside E's package and must be controlled by an independent evaluator.

The extxyz structures use PBC on all axes, with 19.365 A minimum vacuum projected along the surface normal. The slab is as-cleaved from the CIF, non-bulk-stoichiometric at its symmetric Ti termination, and has not been relaxed. No energy or surface stability claim is made.

`build_candidate.py` uses ASE 3.29.0 and Pymatgen 2026.9.24 in the isolated `/tmp/e-surface-ase-venv`; it checks source CIF SHA, stoichiometry, termination, inversion symmetry, vacuum, and minimum pair distance. It does not run DFT or a potential.
