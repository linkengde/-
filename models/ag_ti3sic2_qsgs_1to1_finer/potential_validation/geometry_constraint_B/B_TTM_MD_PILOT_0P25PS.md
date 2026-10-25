# Model B: extended reference TTM–MD pilot (0.25 ps)

Date: 2026-10-25. This extends the bounded 0.05 ps Model B pilot to test the early lattice-temperature response under the same source. It is still an exploratory rigid-obstacle model, not a production TTM calculation or a melt-pool result.

## Setup

- Geometry and initialization are the same as the [0.05 ps pilot](B_TTM_MD_PILOT_0P05PS.md): 916,534 QSGS atoms, Ag EAM, fixed Ti₃SiC₂/WCA geometric obstacle, relaxed-parent coordinates, and seeded 300 K mobile Ag velocities.
- The TTM reference solver uses the 10 Å grid, `R=30 Å`, `l_skin=20 Å`, `I0=5 eV/(ps·Å²)`, candidate Ag tables, and the opt-in sub-300 K continuation. Only Newton-layer Ag is represented in the electron grid and receives lattice energy exchange.
- The source runs for 500 MD steps at 0.5 fs each (0.25 ps). There is no liquid-cluster classifier, live heat-off controller, or cooling stage.
- Local lattice temperature is reported as the source-power-weighted mean over cells receiving at least 1% of the maximum cell source, alongside the maximum coarse-grained temperature among active Ag cells. The latter is a sparse-grid extreme and is not used as a melt criterion.

## Results

| Time | Newton-layer Ag kinetic temperature | Source-weighted Ag lattice temperature | Maximum active-cell Ag lattice temperature | Maximum electron temperature | Cumulative source | Electron-to-lattice transfer | Energy residual |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.050 ps | 158.951 K | 170.080 K | 437.398 K | 985.638 K | 146.737 eV | 4.644 eV | 0.101330 eV |
| 0.100 ps | 178.846 K | 193.874 K | 814.342 K | 1132.064 K | 293.474 eV | 16.147 eV | 0.042184 eV |
| 0.150 ps | 209.520 K | 225.008 K | 592.836 K | 1253.580 K | 440.211 eV | 29.800 eV | 0.059560 eV |
| 0.200 ps | 232.461 K | 249.232 K | 688.306 K | 1359.100 K | 586.948 eV | 46.052 eV | 0.044488 eV |
| 0.250 ps | 229.457 K | 253.204 K | 746.766 K | 1454.081 K | 733.685 eV | 64.309 eV | 0.058903 eV |

The Newton-layer Ag kinetic temperature reached a minimum of 129.277 K at 0.070 ps, then rose to 240.160 K at 0.220 ps before ending at 229.457 K. The source-weighted local lattice temperature ended at 253.204 K. The run therefore shows gradual energy transfer but no thermal evidence of a stable Ag melt pool. The maximum active-cell value is spatially sparse and cannot establish melting by itself.

The corrected final combined-energy residual is `0.05890298 eV`. A no-source NVE control with the same initial coordinates, velocity seed, timestep, and 500-step duration drifted by `0.05795181 eV`, a difference of `0.00095117 eV`. This indicates that the observed residual is dominated by the LAMMPS MD integration drift, while the TTM-only one-step ledger continues to conserve energy to numerical precision.

## Decision and limits

This longer pilot did not reach a stable melt pool under the current source and B setup. It does not decide whether stronger or longer heating could melt Ag: the B model still lacks an accepted cold baseline, and the geometry-only WCA wall and rigid ceramic are not physical Ag–Ti₃SiC₂ interface parameters. The velocity-rescaling coupling is deterministic mean-field exchange, not a validated production TTM method. No stable-pool controller, heat-off point, 100 ps trajectory, or resolidification result is present.

## Reproduction

The 500-step run was made with the archived driver and environment override:

```bash
LD_LIBRARY_PATH=/workspace/.venvs/lammps/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH} \
PYTHONPATH=/tmp/ag-ttm-mpich:/workspace/.venvs/lammps/lib/python3.12/site-packages \
OMP_NUM_THREADS=4 AG_TTM_PILOT_STEPS=500 python \
models/ag_ti3sic2_qsgs_1to1_finer/potential_validation/geometry_constraint_B/b_ttm_reference_md_pilot.py
```

The matching 500-step no-source control uses `AG_NVE_BASELINE_STEPS=500` with `b_nve_energy_baseline.py`. Captured output is archived in `log.b_ttm_reference_md_pilot_0p25ps` and `log.b_nve_energy_baseline_0p25ps`; the original 100-step records remain separate.
