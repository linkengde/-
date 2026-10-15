# DFT execution decision — M77 2.40 Å resource pilot

```text
DFT_EXECUTION_READY=NO
DECISION=HOLD
DFT_STARTED=NO
```

The D6 input, structure identity, ASE/GPAW versions, PAW setup hashes, and 4×4×1 k-point descriptor pass read-only checks. The current cgroup has only 6,026,514,432 bytes (~5.62 GiB) available under a 32-GiB limit, a historical peak within 73,728 bytes of that limit, and 15 memory-limit events. No M77 peak-memory benchmark exists. Therefore this environment fails the resource gate for even a single SCF pilot.

SCF initialization, force/stress production, and full output integrity could not be tested without crossing the failed resource gate. Static output fields are present in the D6 runner, but that does not count as a completed calculation.

**Decision:** do not execute DFT here. Wait for a fresh resource preflight on a host with sufficient measured headroom and durable scratch/archive capacity, then request a separately bounded single-point release. Keep the D6 input parameters, structure, lattice, PBC, and atom order unchanged.
