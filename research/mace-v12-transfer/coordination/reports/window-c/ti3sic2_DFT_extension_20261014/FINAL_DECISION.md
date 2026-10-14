# Stage C final decision

## Status

```text
DFT_PLAN_READY=NO
STRUCTURE_VALIDATION_PASS=PARTIAL
RESOURCE_AVAILABLE=NO
START_DFT=NO
```

## Findings

- E6 requests 14 static points: 4 bulk volume/shear, 2 Si z responses, 4 clean-surface points, and 4 independent DEV/TEST points.
- Ten source-backed geometries were prepared and passed the geometry preflight: six 48-atom bulk structures and four 100-atom symmetric clean slabs. Their source/parent hashes, local atom IDs, cell/PBC, perturbation metadata, minimum distances and geometry SHAs are recorded in `DFT_STRUCTURE_LIST.csv`, `geometry_validation.json` and `SHA256SUMS.txt`.
- The independent DEV and evaluator-held TEST parent structures cannot be sourced from the reviewed E6 bundle. They are listed as four uncreated, source-pending slots. They must not be derived from the TRAIN parent, relabeled from different literature years, or fabricated as random perturbations.
- All ten prepared structures are grouped as one baseline-derived family and only proposed for TRAIN pending explicit role freeze. Existing numerical-control labels and E archives retain their current roles. No TRAIN/DEV/TEST label assignment is frozen by this package.
- DFT code and PAW setup versions match the archived method, and setup hashes match. Actual slab k convergence and resource feasibility have not been demonstrated.
- The current cgroup has ~6.01 GB live headroom and a recorded peak within 73,728 B of its 32 GiB limit. Resource availability for a DFT job is therefore not certified.

## Required before a future DFT release

1. Approve the proposed family-level role handling without changing existing numerical-control records, or keep all ten candidates on hold.
2. Identify a genuine DEV parent and a separate evaluator-controlled TEST parent, each with auditable source lineage and no overlap with TRAIN or each other; then build reference/response descendants and audit their cells, IDs, mesh and hashes.
3. Complete a separate surface k-point convergence and resource pilot decision. Do not change the E6 numerical method to work around resources.
4. Obtain a fresh, suitable DFT allocation with enough memory and persistent scratch/backup for the exact structures.
5. Issue a separate explicit DFT execution authorization. This package itself is not authorization.

No DFT, MACE/TSCBOP fitting, MD, TTM, potential edit, or E-archive modification was performed. Stop here and wait for further authorization.

