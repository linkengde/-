# Model B: bounded cold-state diagnostic (0.25 ps)

Date: 2026-10-25. This is the single bounded no-source gate requested before any further heating. It uses the current fixed-scaffold geometry-obstacle Model B, with the same minimized Ag coordinates, 300 K velocity seed, 0.5 fs timestep, and Ag EAM/WCA potential as the TTM pilots. No thermostat, TTM, or external source is applied during the 500-step run.

## Result

| Quantity | Start | 0.25 ps end | Change |
|---|---:|---:|---:|
| Mobile-Ag kinetic temperature | 300.00 K | 227.29 K | −72.71 K |
| Total potential energy | −537983.480 eV | −536301.644 eV | +1681.836 eV |
| Total kinetic energy | 6950.770 eV | 5268.993 eV | −1681.778 eV |
| Total energy | −531032.709 eV | −531032.651 eV | +0.058 eV |
| Maximum mobile-Ag force | 1.452 eV/Å | 4.685 eV/Å | increased |

The total-energy drift is small and matches the separately measured no-source NVE integration drift. However, the large potential-to-kinetic energy conversion and rising maximum Ag force show that this starting structure is still undergoing substantial mechanical/structural evolution; it does not meet the cold-state gate for interpreting subsequent heating as Ag melting.

The initial-to-final mobile-Ag RMS displacement is `0.286 Å` (p90 `0.431 Å`, p99 `0.912 Å`, maximum `2.236 Å`). In the initial heat-source core (1275 Ag atoms within lateral radius 30 Å and surface depth 0–60 Å), RMS displacement is `0.380 Å` (p90 `0.590 Å`, p99 `1.050 Å`). The mobile-Ag per-atom potential-energy sum rises by `1668.033 eV`, while the same hot-core sum changes by `−0.996 eV`; most of the global potential-energy change is therefore outside that geometric hot core in this diagnostic.

CNA unknown-code fractions are already high at the starting snapshot: 71.06% over mobile Ag and 90.20% in the hot core. They change by less than 0.1 percentage point over this run. Because this is a composite geometry with many interfaces, CNA unknown alone cannot identify liquid Ag. The provisional composite liquid flag also needs a 1 ps displacement window, longer than this cold diagnostic.

At the provisional 3.98 Å geometric Ag-connectivity cutoff, the largest component contains 179,019/179,246 atoms initially and 179,032/179,246 at the end. This is a connectivity statistic only, not evidence of a liquid pool.

## Decision

**Cold-state gate: not passed.** Following the agreed plan, do not start the 1 ps nominal-source heating pilot from this initial state. Do not continue an open-ended FIRE/NVT relaxation sequence. The next step must address the B initial structure once with explicit stopping criteria, or mark this B preparation as unusable for a melting interpretation. No melting, stable-pool, heat-off, or solidification result was produced by this diagnostic.

## Reproduction and archived data

Run from the repository root with the repository MPICH-linked LAMMPS Python environment:

```bash
LD_LIBRARY_PATH=/workspace/.venvs/lammps/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH} \
PYTHONPATH=/tmp/ag-ttm-mpich:/workspace/.venvs/lammps/lib/python3.12/site-packages \
OMP_NUM_THREADS=4 python \
models/ag_ti3sic2_qsgs_1to1_finer/potential_validation/geometry_constraint_B/b_cold_diagnostic_0p25ps.py
python models/ag_ti3sic2_qsgs_1to1_finer/potential_validation/geometry_constraint_B/analyze_b_cold_diagnostic.py
```

The captured run is `log.b_cold_diagnostic_0p25ps`; the reproducible metrics are in `cold_diagnostic_metrics.json`. Start/end mobile-Ag snapshots are `cold_diagnostic_start.dump.gz` and `cold_diagnostic_final.dump.gz`.
