# Ti3SiC2 prismatic-surface DEV candidate (2026-10-20)

This directory contains one newly constructed 96-atom prismatic Ti3SiC2 slab, its symmetric reference/displacement pair, the locked fixed-geometry GPAW calculation settings, and a resource-guarded runner.

The surface orientation is conventional hexagonal `(100)`, equivalent to four-index `(10-10)`. The slab is cut from the archived COD 9009647 CIF after a half-`a` translation, with eight layers and 18 Å vacuum. No atom, cell, or strain relaxation was applied. The two input geometries share an identical cell, composition, atom order, and source IDs. The displacement case moves the two terminal Ti atoms in opposite z directions by 0.02 Å each.

The prismatic orientation and geometric family differ from the basal (0001) TRAIN slabs, but both derive from COD 9009647. These calculations can support a geometry-family DEV candidate; they are not an independent crystallographic-source blind test. `TRAIN/DEV/TEST` assignments are not changed here, and no TEST identity or label is included or accessed.

The runtime already contains GPAW 26.7.0, ASE 3.29.0, GPAW data 1.2.1, and a four-rank GPAW C-MPI backend. No software download is required. Each calculation is a fixed-geometry PBE/PAW plane-wave single point at 500 eV, standard Monkhorst-Pack `12x2x1`, Fermi-Dirac width 0.10 eV, and the recorded SCF tolerances. No relaxation or MD is performed.

The runner uses four MPI ranks and one OpenMP thread. It writes a full-wavefunction, single-precision checkpoint every 20 SCF iterations, retaining only the latest complete file. The checkpoint size estimate is scaled from a comparable archived GPAW checkpoint; actual memory and runtime are measured during each point. A safe stop or failed checkpoint produces no final energy/force/stress label. No automatic parameter changes or retries are enabled. Results, logs, resource profiles, and checkpoint hashes are retained under `results/<structure_id>/attempt-0001/`; large temporary checkpoints are not part of the source package and must not be committed to ordinary Git.

The two points are scheduled sequentially: reference first, then the displacement response. If a point stops without convergence, keep its latest checkpoint and logs for a later same-configuration resume; do not start the second point while the first is still running or unresolved. After successful result checks, archive the labels and evidence with SHA-256, push only this package's files to `cloud-e-b1-dft`, and verify by a fresh fetch.
