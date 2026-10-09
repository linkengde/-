# Ti₃SiC₂ c-axis k-point increment check

This report compares the same frozen 48-atom Ti₃SiC₂ structure at **6×6×2** and **6×6×4**. The in-plane mesh is fixed, so this pair isolates the change in the c-axis sampling. Both GPAW archives passed their published checksums and archive verification; geometry, atom IDs, input hash, and all method settings other than `kpts` match exactly.

| Metric, 6×6×4 minus 6×6×2 | Difference | Budget | Result |
|---|---:|---:|---|
| Native energy | 0.31834 meV/atom | ≤ 2 meV/atom | Pass |
| Free energy | 0.15785 meV/atom | ≤ 2 meV/atom | Pass |
| Force-vector RMS | 0.002884 eV/Å | ≤ 0.01 eV/Å | Pass |
| Maximum per-atom force-vector difference | 0.004881 eV/Å | Reported, no separate gate | — |

The z-mesh increment is small for this fixed structure and 6×6 in-plane mesh. This does not establish overall mesh convergence or physical/model accuracy. A's registered 8×8×2 run is still needed to isolate in-plane sensitivity; after it completes, the four mesh points can be considered together.

Machine-readable metrics and archive checksums are in [`report.json`](report.json). Recompute from the archived files with:

```bash
/workspace/.venvs/gpaw-mpi/bin/python compare_k6x6x2_k6x6x4.py --output /tmp/kz_mesh_report.json
```

The script refuses to overwrite an existing output file.
