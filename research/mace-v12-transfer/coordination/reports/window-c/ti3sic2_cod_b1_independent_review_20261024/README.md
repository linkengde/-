# Independent review of the six COD-derived Ti₃SiC₂ geometry proposals

The six B1-TSC files pass an independent geometry and file-integrity review. Their two unperturbed supercells reconstruct directly from the SHA-verified COD 9009647 CIF, and the four displaced files preserve the corresponding parent cell and atom ordering. All six are periodic, finite, stoichiometrically correct, label-free, and explicitly DFT-disabled. Their shortest distinct-atom distances range from 1.973 to 2.088 Å.

The review does **not** freeze roles or authorize labels. Every candidate inherits COD 9009647. A 72-atom cell, a 96-atom cell, or a random displacement seed does not create a new crystallographic source family. Keep all six together as one lineage group and exclude them from independent development/test evidence. They may only be considered as same-source training proposals after the source and role policy is explicitly approved.

This set is separate from the other six Ti₃SiC₂ geometries derived from Bai et al. (2010) in `independent_ti3sic2_bai2010_candidate_20261024/`. That candidate set has a distinct computational-paper lineage, but the exact numerical parameters still need full-text verification. Neither set is currently a frozen B1 role assignment.

No DFT, MPI, MACE inference/training, MD, or TTM was run for this review. The active Window-A task card permits only its registered V20 slab-method controls; these six geometries are not registered DFT inputs. Consequently, this review does not start the proposed B1 pilot.

## Reproduction

From the repository root, run:

```bash
/workspace/.venvs/mace/bin/python research/mace-v12-transfer/coordination/reports/window-c/ti3sic2_cod_b1_independent_review_20261024/audit_cod_tsc_candidates.py
```

The script verifies the source CIF and six candidate file hashes, rebuilds both parent supercells from the CIF, checks species/cell/coordinate identity, checks displacement lineage and minimum pair distances, and confirms that no energy/force labels or DFT launch flag are present. It writes `independent_review.json` in this directory.

The candidate inputs are reproduced byte-for-byte from the six local B1-TSC proposal files. The manifest records each input SHA256 and the source CIF SHA256. They remain geometry-only proposals, not dataset rows.
