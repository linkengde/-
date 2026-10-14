# Stage C-D external run preparation

This directory packages instructions and evidence for a future external compute allocation. It does not contain a DFT execution authorization and does not change the 15 frozen inputs.

Status:

```text
DFT_INPUT_READY=YES
DFT_LABEL_READY=NO
WAITING_FOR_COMPUTE_RESOURCE=YES
SCF_EXECUTION_AUTHORIZED=NO
```

Use `DFT_RUN_INSTRUCTIONS.md` for the pinned environment, method, pilot meshes, and required outputs. `MEMORY_REQUIREMENT_ESTIMATE.md` records the known cgroup limit and unresolved peak-memory requirement. `STRUCTURE_CHECKSUMS.txt` lists the exact 15 geometry file hashes relative to the repository root; run its verification from the repository root.

No input geometry was copied or rewritten for this document package; checksums resolve to the original fixed input files in the DFT input package. The DFT package configs continue to have `execution_authorized=false`.
