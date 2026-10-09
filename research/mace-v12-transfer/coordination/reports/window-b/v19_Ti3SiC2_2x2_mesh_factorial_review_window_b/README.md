# B assignment: Ti₃SiC₂ four-point k-mesh decomposition

**Status: analysis in progress. Analysis only; no DFT or MPI launch.**

B should start by auditing the three published points: 6×6×2, 6×6×4, and 8×8×4. Verify their manifests and archive records, confirm the same frozen geometry and method, and independently recompute the available pairwise metrics. Prepare a reproducible script for the complete 2×2 mesh grid.

The fourth point, **8×8×2**, is running under A. It remains missing until its compact archive appears on `main` and passes verification. The current [`report.json`](report.json) records only verified available points and does not infer missing values. The independent audit found:

| Contrast (second minus first) | Native ΔE (meV/atom) | Free ΔE (meV/atom) | Force-vector RMS (eV/Å) | Result |
| --- | ---: | ---: | ---: | --- |
| 6×6×2 → 6×6×4 | +0.3183 | +0.1579 | 0.002884 | All stated budgets pass |
| 6×6×2 → 8×8×4 (mixed directions) | +0.2118 | +0.0665 | 0.013593 | Energy passes; force RMS fails |
| 6×6×4 → 8×8×4 | −0.1065 | −0.0913 | 0.010850 | Energy passes; force RMS fails |

The mixed 6×6×2 → 8×8×4 contrast cannot isolate either mesh direction. The force differences in the two same-axis contrasts exceed the 0.01 eV/Å budget for 6×6 → 8×8 at kz=4, but the direction-specific in-plane and interaction conclusions remain pending 8×8×2.

The reproducible script verifies archive checksums, result-verifier status, source and method identity, and frozen geometry before calculating pairwise and factorial contrasts:

```bash
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/coordination/reports/window-b/v19_Ti3SiC2_2x2_mesh_factorial_review_window_b/analyze_mesh_matrix.py
```

After A publishes the verified fourth point, rerun the same command. It will replace the partial report with the full four-point decomposition:

- in-plane changes 6×6→8×8 at kz=2 and kz=4;
- c-axis changes kz=2→4 at 6×6 and 8×8;
- the difference-in-differences interaction for native energy, free energy, and per-atom force vectors;
- force RMS, maximum per-atom difference, Cartesian component and species breakdowns;
- pass/fail against 2 meV/atom energy budgets and 0.01 eV/Å force-vector RMS.

Keep conclusions limited to this frozen pure-phase numerical control. Do not launch another DFT job, training, gap scan, or dynamics, and do not change A's ownership. Keep the gap scan disabled. Publish a brief status before starting and hash-verify each report update.

The exact assignment is in [`assignment.json`](assignment.json).
