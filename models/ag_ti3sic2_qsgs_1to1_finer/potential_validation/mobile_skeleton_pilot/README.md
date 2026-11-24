# Mobile Ti₃SiC₂ skeleton exploratory pilot

This independent pilot addresses whether the Ti₃SiC₂ skeleton can remain mobile above a fixed bottom slab. It does not modify the accepted QSGS model or any earlier Model B files.

## Setup

- Geometry: the separate TSCBOP-relaxed-parent QSGS candidate, using the Ag-only FIRE-100 endpoint from the geometry-obstacle Model B trial.
- Fixed atoms: all atoms initially at `0 ≤ z < 20 Å`.
- Mobile atoms: every Ag, Ti, Si, and C atom above that slab. In particular, the Ti₃SiC₂ skeleton is not frozen.
- Boundary: `p p f`.
- Dynamics: 300 K Maxwell velocities followed by 0.025 ps NVE, `dt = 0.25 fs`, without a heat source or thermostat.
- Interface: the existing exploratory WCA geometric wall (ε=0.5 eV, no attraction). It is not a calibrated Ag/Ti₃SiC₂ interface potential and cannot support quantitative adhesion or interface-conductance claims.

## Purpose and limits

This is a mobility/stability smoke test. It checks whether the newly mobile ceramic immediately undergoes large structural motion and whether the run completes without force or energy blow-up. It is deliberately shorter than an equilibration or heat-transfer simulation. The fixed bottom is a mechanical anchor, not a thermal reservoir. No melting/cooling, TTM, or heat-flux measurement is included.

Exploratory work does not need publication-level thresholds, but the initial state still needs to remain numerically bounded over this short interval. If it does not, that is a useful failure result and the run should not be extended into heating.

## Reproduction

Run from `models/ag_ti3sic2_qsgs_1to1_finer/`:

```bash
OMP_NUM_THREADS=4 /workspace/.venvs/ag-qsgs/bin/lmp -sf omp \
  -in potential_validation/mobile_skeleton_pilot/in.mobile_skeleton_nve_0p025ps \
  -log potential_validation/mobile_skeleton_pilot/log.mobile_skeleton_nve_0p025ps
```

The input writes start/end atom snapshots and a log. The analysis report will state whether the short run is stable enough to justify a separate temperature-gradient test; it will not claim that the skeleton is a validated heat sink.

## Result

The 25 fs NVE run completed for 916,534 atoms. The bottom group contained 113,499 atoms; 623,789 skeleton atoms and 179,246 Ag atoms were mobile. Total energy drift was 0.60 eV (about `1.3e-7` relative to the magnitude of the initial total energy), so the integration did not show an energy blow-up.

The structure did not begin from a stationary cold state: potential energy increased by 3,944.5 eV while kinetic energy fell by 3,943.9 eV. The mobile-group temperature fell from 300 K to 262 K; the final skeleton and Ag temperatures were 265 K and 252 K. Maximum moving-atom force rose from 3.42 to 6.14 eV/Å. The largest displacement was 0.61 Å (on C); no mobile atom moved more than 1 Å in 25 fs. The skeleton therefore had mobile vibrational degrees of freedom and remained bounded over this short screen, but continued structural relaxation dominates the energy change.

**Decision:** this answers the mechanics question provisionally: fixing only the bottom while allowing the upper skeleton to move is implementable, and it avoids the rigid-scaffold assumption. It does **not** establish a stable 300 K baseline, a physically valid Ag/ceramic heat-transfer rate, or a melt/cool result. Do not start a temperature-gradient or TTM run from this endpoint. The next useful step is to improve the cold starting geometry/interface model, then repeat a short no-source screen before adding a controlled hot-Ag/cool-skeleton gradient.

Reproduce from the model directory with the command above. The archived `.dump.gz` snapshots are the gzip-compressed output of the same `write_dump` commands. Run the analysis after compression:

```bash
gzip -f potential_validation/mobile_skeleton_pilot/start.dump
gzip -f potential_validation/mobile_skeleton_pilot/final.dump
/workspace/.venvs/ag-qsgs/bin/python potential_validation/mobile_skeleton_pilot/analyze_mobile_skeleton_pilot.py
```
