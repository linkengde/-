# Ag 70:30 vol% micro–nano comparison design

Design date: 2026-10-28  
Repository baseline: `7035293a06a70c2515d3ab2c973417e48e5a1958`  
This is a study design only. No new geometry or transient was generated.

## Question and present evidence

The intended question is whether changing the Ag/Ti₃SiC₂ bicontinuous feature scale changes local thermal localization and Ag melt response under the same contact-closing load. The current 70:30 vol% calculation is a single QSGS realization, 69.9996 vol% Ag. Its Case A 1× run remained fully solid through 1.0 ps (maximum 906.21 K under the candidate Ag melting reference 1235.08 K); the available 2× field checkpoint ends at 0.35 ps and also remains solid. The older 41.6654 vol% calculation is a different composition and geometry. Neither result is a micro/nano comparison or evidence of a material anti-welding ranking.

The comparison must first separate three effects that can otherwise be confounded: feature size, stochastic morphology, and source absorption/partition. A single repeated tile or a single QSGS seed does not quantify morphology variability.

## Model pairs and controlled factors

### Factors held equal in the primary pair

- Ag:Ti₃SiC₂ volume fraction = 70:30, with porosity stated separately and the same volume-fraction convention.
- Same physical outer specimen dimensions, initial temperature, material property functions, interface law, contact/outer boundary conditions, and physical simulation interval.
- Same experimentally constrained electrical/contact or arc loading history and the same source footprint in physical units. Do not retain the 3 nm × 2 nm source merely because it is in the current script; use it only as a named exploratory case until calibrated.
- Same liquid threshold, six-face connectivity convention, cluster tracking, sampling interval, energy ledger, and numerical convergence criteria.
- Same QSGS algorithm family and phase-generation rules, while changing only the correlation/feature scale to independently measured nano and micro target distributions. Record every random seed and input rank field.
- Same per-area interface conductance in the isolation study; measure or separately justify any scale-dependent interface chemistry/contact change rather than silently changing it with feature scale.

### Two source protocols to avoid a hidden confound

1. **Physical-load comparison (primary):** impose the same measured or macroscopically calibrated contact/arc event in both models. Compute electrical/thermal deposition by the same coupled source model. Differences in phase-local absorption are part of the physical response and must be reported.
2. **Equal-Ag-deposited-energy control (diagnostic):** normalize each mask's Ag-deposited energy to an identical value and preserve the same spatial/time profile. This isolates thermal storage/transport and topology from differences in how much prescribed source overlaps Ag. Label this as a controlled thermal-response comparison, not the physical event.

The current source is not calibrated and its integrated power changes with the QSGS mask; the recent 41.67% versus 70% contrast cannot be interpreted as a clean morphology or composition experiment. For scale comparisons at fixed 70:30 composition, both protocols make the source-overlap effect explicit.

## Geometry characterization and representative volumes

Use measured phase morphology to choose the target scales: polished-section SEM for microstructure, and TEM/serial-section or an appropriate 3D method for nanoscale channels and ligaments. Determine Ag-channel and Ti₃SiC₂-ligament distributions, phase connectivity, interface area per volume, anisotropy, porosity, and the smallest necks that must be resolved. Do not assign a “micron” or “nano” voxel scale from the QSGS grid alone.

For each scale, prepare multiple independent realizations at the target 70:30 volume fraction and compare distributions, not only one seed. Select representative volume size by increasing the sampling volume until the key morphology descriptors and effective thermal response plateau within the realization spread. A 3×3 repeat of one tile is a periodic-size check for that tile, not an independent realization or RVE proof.

Preserve measured feature-size distributions where possible; do not create the micro model by simply relabeling the existing nanometer coordinates. Record phase-mask hashes, voxel size, physical box, seed, interface-area estimator, and any morphology filtering. Resolve the narrowest relevant neck by several cells and demonstrate that the chosen resolution does not change the descriptors or transport response materially.

## Thermal and electrical model hierarchy

1. **Property/source closure:** attach traceable Ag and Ti₃SiC₂ thermal data and calibrate the contact load; distinguish intrinsic properties from porous bulk values. Assess bulk Ag conductivity at nanoscale: surface/interface scattering may invalidate applying a bulk constant unchanged when feature widths approach electron transport lengths. Add size-dependent properties only with evidence.
2. **Interface limits:** retain current Case A as an artificial insulating-scaffold limit. Retain current Case B as perfect Ag/Ti₃SiC₂ thermal contact with both phases conducting and storing energy. Implement a finite interfacial conductance model between those limits, with equal/opposite heat exchange and a parameter range supported by measurement or literature. Do not fit that range to achieve a preferred melt outcome.
3. **Electrical load:** for a resistive contact, solve or map `J·E` using the same measured waveform, contact geometry, and electrical property model. Treat arc deposition separately. Compare with a source-off thermal stage under a specified external electrode/sink path. If a directly prescribed thermal source is retained as a reduced model, report the energy delivered, footprint, phase partition, and calibration basis.
4. **Thermal time-scale check:** before using direct single-temperature lattice heating for a picosecond-scale load, compare load and observation times against a sourced electron–lattice relaxation time and compare gradients against the validity of local Fourier transport. If those conditions are not met, use the independent electron/lattice model as a benchmark in a separate workstream; this design does not alter or run that project.

## Scale-aware computational strategy

The strict direct comparison uses the same physical specimen extent and same physical hot-source footprint. It must resolve the smallest nano-scale feature everywhere, so cell count grows approximately as `(L/Δx)^3`; feasibility must be estimated before launching the full pair.

If that full direct mesh is unaffordable, use a staged comparison and state exactly what is being compared:

- First compute transport and morphology statistics on independently representative micro- and nano-RVEs.
- Derive scale-specific effective conductivity tensors, heat capacity, and interface-area response with uncertainty.
- Apply the same macroscopic contact heat history and external boundaries in a common specimen-scale model.
- Where a local melt prediction is required, embed a resolved source-centred microstructure submodel under matched temperature/flux histories from the specimen model, or run equal-size physical local domains under the same calibrated load.

An RVE-homogenized specimen model and a fully resolved source-centred model answer different questions. Do not present the homogenized field as a directly resolved melt-pool topology; do not present a small local RVE as the whole contact unless the boundary history and scale separation are justified.

## Required comparisons and acceptance checks

For each independent realization and each declared boundary/interface case, report:

- Ag fraction, porosity, interface-area density, channel/ligament chord distributions, connectivity, anisotropy, and narrow-neck resolution.
- Measured/calibrated input power and deposited energy versus time, with a phase-resolved and whole-domain energy ledger.
- Peak temperature and hotspot location, while reporting the spatial-resolution sensitivity of the maximum.
- Source energy at first nonzero Ag liquid fraction and at first multi-voxel six-face connected liquid cluster.
- Total Ag liquid volume fraction and largest connected cluster's physical volume, centroid/depth, and lateral extents.
- Cluster merge/split identity and residence interval using overlap tracking; report first-liquid and connected-pool times separately.
- Source-off liquid evolution and time to full solid only under the named external boundary/sink condition. Under an insulated boundary, call it internal energy redistribution and do not call it environmental cooling time.
- Time-step and spatial refinement, RVE-size and realization spread, and source/interface/property sensitivity. A mesh that agrees only while both models are fully solid does not validate a later melt-pool size or life.

The 0.3 ps and 0.5 ps continuity screens previously used are exploratory persistence criteria, not universal stability thresholds. A liquid cluster alone does not prove wetting, metallurgical joining, mechanical adhesion, or a successful/failed weld. Those require additional interfacial and mechanical physics and experimental validation.

## Comparison table template

| Metric | Nano distribution | Micro distribution | Matched protocol / uncertainty |
|---|---:|---:|---|
| Ag / Ti₃SiC₂ volume fraction and porosity | pending | pending | same basis; measured uncertainty |
| Interface area per volume | pending | pending | same estimator and resolution study |
| Ag-channel / Ti₃SiC₂-ligament distributions | pending | pending | measured 3D target and independent realizations |
| Effective conductivity tensor and volumetric capacity | pending | pending | same constituent data; interface sensitivity |
| Energy to first liquid / connected liquid | pending | pending | physical-load and equal-energy protocols distinguished |
| Maximum liquid fraction / largest cluster volume | pending | pending | physical volume, not voxel count |
| Cluster centroid, depth, extents, and persistence | pending | pending | same overlap tracker and output cadence |
| Source-off full-solid time | pending | pending | only with named heat-sink boundary and energy ledger |
| Mesh / time-step / RVE / realization uncertainty | pending | pending | show convergence and ensemble spread |

## Go/no-go sequence

1. Do not start micro/nano transient runs until the heat-source-to-contact mapping, candidate material sources, and target morphology distributions are documented.
2. Verify each morphology's fraction/connectivity/statistics and estimate direct-grid cost. Keep existing QSGS and checkpoints immutable.
3. Validate source integration, total energy, interface exchange, and limiting cases on small tests before any long calculation.
4. Run a short same-load matched pair only after mesh/time-step/RVE preflight. Continue to melt/cool branches only if the pair remains numerically valid and the constitutive data cover the reached temperatures.
5. Make paper-level claims only after physical-input, material, spatial, temporal, and realization uncertainties are bounded. A difference between two single seeds or between the old 41.67% and current 70% models is not a scale effect.

