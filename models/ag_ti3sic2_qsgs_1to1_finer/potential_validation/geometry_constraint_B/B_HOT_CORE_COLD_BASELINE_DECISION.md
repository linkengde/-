# Source-core-only Model B cold-baseline decision

Date: 2026-10-25. After the whole-layer and 75 Å local-ROI cold checks failed, one last bounded test allowed only atoms inside the nominal source core to move. The goal was to determine whether distant Ag rearrangement was masking a usable local starting state.

## Setup

- Start: archived Ag-only FIRE-100 coordinates.
- Mobile region: 1,275 Ag atoms within lateral radius 30 Å and surface depth 0–60 Å.
- All other Ag and all Ti/Si/C atoms were held fixed under the geometry-only WCA obstacle model.
- NVE at 0.5 fs for 200 steps (0.1 ps), no external source and no TTM.
- Reproduction from the finer-model directory:

```bash
OMP_NUM_THREADS=4 /workspace/.venvs/lammps/bin/lmp -sf omp \
  -in potential_validation/geometry_constraint_B/in.hot_core_cold_nve_0p1ps \
  -log potential_validation/geometry_constraint_B/log.hot_core_cold_nve_0p1ps
python potential_validation/geometry_constraint_B/analyze_hot_core_cold_nve.py
```

## Results

- Ag kinetic temperature: 300.00 K → 223.38 K.
- Potential energy: +12.62 eV; kinetic energy: −12.62 eV; total-energy drift: below printed precision (`<0.01 eV`).
- Maximum force: 1.090 → 1.594 eV/Å; force p95: 0.757 → 0.872 eV/Å.
- Source-core displacement: 0.197 Å RMS, 0.296 Å p90, 0.434 Å p99, 0.588 Å maximum.

The small energy drift confirms stable integration. The core still cools substantially while converting potential energy into kinetic energy, and its force tail increases. It is therefore not demonstrated to be a stationary 300 K cold baseline.

## Decision

**Stop the current geometry-only Model B route. Do not heat this state.** Restricting motion first to a 75 Å ROI and then to the 1,275-atom source core does not remove the cold-state relaxation. No stable melt pool or local solidification result has been produced by these B tests.

Further FIRE/NVT cycles on these coordinates are not warranted. Continuing the research objective requires a materially different local atomic geometry or a better-supported Ag–scaffold interaction model, followed by a new finite cold-state gate. All three bounded B variants and their data are retained for reproducibility.
