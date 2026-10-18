# Reduced MACE V17 local-heating run

## Result

- Model atoms: 1084 (Ti₃SiC₂=600, Ag=484); Ag mass fraction 72.73%. 
- Time step: 0.5 fs; local-heating duration: 1.000 ps; cooling duration: 0 ps.
- Heating wall time including static force and 300 K baseline: 5.34 h (19223.7 s).
- Ag melt established: NO; maximum connected low-q₆ Ag cluster: 9 atoms at 0.970 ps.
- Maximum low-q₆ fraction: 14.9%; minimum hot-region median q₆ 0.422 (300 K baseline 0.594); maximum hot-region RMS displacement 1.108 Å.
- Maximum mobile force 14.03 eV/Å; minimum TSC connectivity fraction 1.000.

The trajectory reached 1 ps under a 1900 K local Langevin target, but the measured Ag remained predominantly ordered: q₆ did not fall below 70% of its baseline for three samples, fewer than half of hot Ag atoms had q₆ < 0.35, the largest connected low-q₆ cluster was below 20 atoms, and RMS migration stayed below 1.5 Å. Therefore this run shows heating and local disorder only, not a demonstrated melt pool. Cooling was deliberately not run because no liquid state was established.

## Boundary and model limits

The reduced geometry keeps the original 30.575 × 26.479 Å lateral coincidence cell and ideal Ag/Ti₃SiC₂ interface. It contains a 7.54 Å six-plane stoichiometric Ti₃SiC₂ slab, four Ag(111) planes, 18 Å top vacuum, XY periodic boundaries, only the bottom Ti plane fixed, and all other TSC atoms mobile. It is a layered interface patch rather than the original bicontinuous 39,486-atom model. Ag/TSC mass balance shifts to 72.73 wt% Ag.

The first 0.30 ps used a 300-atom TSC heat sink. To reduce thermostat clamping in this thin slab, the remaining 0.70 ps used only the two mobile TSC planes immediately above the fixed support (200 atoms) at 300 K; other atoms remained dynamically mobile. The same V17 weights and 0.5 fs step were used throughout.

## Files

- `OVITO_HEAT_TRACE.extxyz`: baseline plus 1 ps heating trajectory; no cooling frames.
- `TEMPERATURE_AND_STRUCTURE.csv`: time series of temperature, force, q₆, coordination, displacement proxy, contacts, and TSC connectivity.
- `AG_LIQUID_CLUSTER_STATS.csv`: per-frame Ag disorder and connected-cluster indicators.
- `HEATING.extxyz`, `HEATING_METRICS.csv`, and checkpoint files: raw staged heating outputs.
- `REDUCED_MODEL.extxyz` / `REDUCED_MODEL.data`: independent reduced slab coordinates.

This result is a MACE V17 exploratory trajectory. No DFT, retraining, potential mixing, or MD with EAM/TSCBOP/Morse was used. The potential remains provisional; non-melting at this finite time does not establish an equilibrium melting temperature.
