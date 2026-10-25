# Pure Ag liquid-flag calibration

This is a reference calibration of the provisional per-atom Ag liquid flag used by `melt_metrics.py`:

```text
CNA unknown code == 5 AND displacement over the preceding 1 ps >= 1.0 Å
```

It uses the repository's unmodified `potentials/Ag_u3.eam`, 4,000-atom FCC Ag (initial lattice parameter 4.09 Å), and a 0.5 fs timestep. Each input runs 10 ps and writes a frame every 0.5 ps. The solid reference is NVT at 300 K; the liquid reference is NPT at 1600 K and 0 bar. The 1500 K fixed-volume run is retained as a cautionary disordered-but-not-classified case.

The reported displacement is reconstructed from wrapped dump coordinates using fractional-coordinate minimum images over a 1 ps window, with affine NPT box expansion removed. The thermodynamic analysis integrates the same affine-corrected frame-to-frame fractional displacements to obtain MSD. RDFs are normalized frame-by-frame using the instantaneous volume and averaged from 5–10 ps. These are pure Ag EAM calibration cases, not validation for Ag in the composite or near the Ti₃SiC₂ interface.

## Results

| Reference | Time | CNA unknown fraction | Displacement p50/p90/p99 (Å) | Liquid-flag fraction |
|---|---:|---:|---:|---:|
| 300 K NVT solid | 1 ps | 0.0000 | 0.140 / 0.228 / 0.305 | 0.0000 |
| 300 K NVT solid | 10 ps | 0.0015 | 0.190 / 0.302 / 0.419 | 0.0000 |
| 1500 K NVT, fixed volume | 10 ps | 0.7408 | 0.488 / 0.783 / 1.096 | 0.0173 |
| 1600 K NPT, 0 bar | 1 ps | 0.9988 | 0.515 / 0.881 / 1.246 | 0.0493 |
| 1600 K NPT, 0 bar | 10 ps | 1.0000 | 1.655 / 2.855 / 4.023 | 0.8240 |

The combined criterion rejects the 300 K solid and the strained 1500 K fixed-volume state while identifying the 1600 K zero-pressure liquid after the initial expansion/transformation. The 1500 K result shows why CNA disorder alone is insufficient. In the 1600 K NPT run, the 5–10 ps mean temperature is 1594.7 ± 22.4 K, mean pressure is −141 ± 1172 bar, density is 8.2848 ± 0.0325 g/cm³, and volume is 86482 ± 339 Å³. The affine-corrected MSD is linear from 2.5–10 ps (slope 3.888 Å²/ps, R²=0.9992; apparent self-diffusion coefficient 0.648 Å²/ps). The mean RDF from 5–10 ps has a first peak at 2.825 Å and first minimum at 3.975 Å. This supports a liquid pure-Ag reference and a provisional 3.98 Å geometric connectivity cutoff.

The measured NPT density is lower than expected from experimental liquid Ag, so it is a property of this EAM calibration and a reminder that the potential is only a structural reference here. The 1.0 Å displacement threshold and 3.98 Å cutoff remain provisional; connected-pool behavior and thresholds must be checked on composite trajectories.

The separate [pure-Ag TTM–MD reference pilot](TTM_MD_PILOT.md) checks short real LAMMPS atom timesteps with source-driven electron-grid updates at two grid resolutions. It verifies numerical energy handoff only; it does not produce a melt pool or validate the composite.

The research target is a **microscopic local Ag melt pool and local resolidification** inside the Ag–Ti₃SiC₂ microstructure. The continuous Ti₃SiC₂ scaffold is intended to confine liquid Ag migration and pool growth, then constrain the local freezing morphology after heat-off. This homogeneous, fully periodic pure-Ag run is only used to check phase metrics; it is not a model of whole-composite melting. Given the known limitations of the Ag–Ti/Si/C cross potential and the Ti₃SiC₂ potential for weld-like thermal cycling, later composite results must be described as conditional on the current potentials, not as validated quantitative confinement or interface-welding predictions.

Plots and processed data: [NPT curves](calibration_1600K_npt.png), [time series CSV](analysis_1600K_npt_series.csv), and [RDF CSV](analysis_1600K_npt_rdf.csv).

## Reproduce

Run from this directory with a LAMMPS build that includes EAM, CNA, and NPT:

```bash
lmp -in in.solid_300K -log log.solid_300K
lmp -in in.liquid_1500K -log log.liquid_1500K
lmp -in in.liquid_1600K_npt -log log.liquid_1600K_npt
python analyze_dumps.py dump.solid_300K dump.liquid_1500K dump.liquid_1600K_npt
python analyze_npt.py dump.liquid_1600K_npt log.liquid_1600K_npt
python plot_npt.py analysis_1600K_npt_series.csv analysis_1600K_npt_rdf.csv
```

For the NPT run, RDF analysis uses SciPy and plotting requires Matplotlib. The dump and large LAMMPS data files are generated outputs and are intentionally not checked in; compact processed CSVs and the plot are included. Random seeds are fixed in the inputs; exact floating-point trajectories can still vary with LAMMPS build and parallel decomposition.
