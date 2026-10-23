# A-approved V19 dataset

This package contains the frozen train36/dev3 split approved by window A after B's read-only audit. `train.extxyz` preserves the existing 36 training rows byte-for-byte; `valid.extxyz` contains the three fixed-geometry AgC, AgSi, and AgTi development labels with native extrapolated energies and forces. There is no test file, and no sealed/test labels were opened.

## Independence boundary

The B audit passed its geometry, hash, method, and archive checks, with a shared-source limitation: all three development structures inherit the same uploaded Ti3SiC2 COD 9009647 CIF. They are distinct from the screened train36, six V19 residual geometries, and inspected V17/V18 development geometries, but they are one crystallographic source family, not three independent source families. The upstream COD bytes remain externally unverified. This split supports development screening, not source-diverse or blind-test generalization claims. Slab k-point, smearing, and vacuum convergence remain unresolved.

See `dataset_manifest.json`, `../entry_review_approval.json`, and `SHA256SUMS.txt` for roles, lineage, approval evidence, and reproducibility hashes.
