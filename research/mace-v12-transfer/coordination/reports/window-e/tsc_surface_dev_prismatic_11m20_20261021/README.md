# Ti3SiC2 (11-20) prismatic DEV geometry-family candidate

A symmetric 120-atom, 10-layer slab from COD 9009647, with 18 A vacuum and an inversion-paired terminal-Ti response geometry. This is a new facet geometry family compared with the existing (10-10) prismatic candidate, but shares the COD crystallographic source. It is not an independent crystallographic-source blind test.

The geometry is unrelaxed and unstrained. Candidate role only; no TRAIN/DEV/TEST assignments are changed. The plan is two sequential PBE/PAW fixed-geometry single points, 500 eV, 7x2x1 standard Monkhorst-Pack, Fermi-Dirac 0.10 eV. The 7x2x1 mesh is selected to preserve reciprocal spacing relative to the previous 12x2x1 mesh; this is not a convergence study.

Run only after the generated manifests, PAW hashes, resources, disk headroom and absence of competing compute jobs pass the supervisor preflight.
