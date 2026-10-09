# B assignment: Ti₃SiC₂ four-point k-mesh decomposition

**Status: assigned. Analysis only; no DFT or MPI launch.**

B should start by auditing the three published points: 6×6×2, 6×6×4, and 8×8×4. Verify their manifests and archive records, confirm the same frozen geometry and method, and independently recompute the available pairwise metrics. Prepare a reproducible script for the complete 2×2 mesh grid.

The fourth point, **8×8×2**, is running under A. Treat it as missing until its archive appears on `main` and passes verification. Then run the script and report:

- in-plane changes 6×6→8×8 at kz=2 and kz=4;
- c-axis changes kz=2→4 at 6×6 and 8×8;
- the difference-in-differences interaction for native energy, free energy, and per-atom force vectors;
- force RMS, maximum per-atom difference, Cartesian component and species breakdowns;
- pass/fail against 2 meV/atom energy budgets and 0.01 eV/Å force-vector RMS.

Keep the conclusion limited to this frozen pure-phase numerical control. Do not launch another DFT job, training, gap scan, or dynamics, and do not change A's ownership. Keep the gap scan disabled. Publish a brief status before starting and a hash-verified report when finished.

The exact assignment is in [`assignment.json`](assignment.json).
