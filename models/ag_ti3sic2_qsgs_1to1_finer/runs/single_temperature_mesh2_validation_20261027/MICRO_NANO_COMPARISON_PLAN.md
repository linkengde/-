# Fair microstructure-size comparison plan

## Fixed factors

- Ag:Ti₃SiC₂ = 70:30 by volume in both structures, with the same density convention and explicit porosity accounting.
- Identical macroscopic specimen dimensions, heat-source space/time profile and deposited energy per heated area (unless the experiment supplies a better calibrated waveform).
- Same single-temperature enthalpy model, phase properties, outer boundary cases, interface-contact assumption, time step criteria, output cadence and liquid-cluster definition.
- Preserve each QSGS realization and record its seed, rank field, phase mask, voxel spacing, periodicity and SHA-256.

## Geometry preparation before heat calculations

The current 70:30 thermal mask is a thresholded realization of the existing QSGS growth ranks, not a new atomistic sample. A micron-scale version should be generated from measured Ag-channel and Ti₃SiC₂-ligament dimensions or serial-section/3D microscopy, not by choosing an arbitrary voxel scale. Match the *physical* volume fractions and compare interface-area density, chord-length distributions, directional anisotropy, connected fractions, Euler/Betti descriptors where robust, and minimum resolved neck distributions. Use multiple independent realizations or measured representative volumes to quantify morphology variability.

For a fair scale test, retain the same outer dimensions while scaling the internal two-phase feature lengths. Refine the voxel grid until native voxel size is comfortably smaller than the smallest measured neck and carry out matched-resolution transport tests. Choose an RVE large enough that the reported descriptors and source-centred thermal response are insensitive to further enlargement; a repeated 3×3 tile alone checks periodic finite-size effects for one realization, not stochastic representativeness.

## Required matched outcomes

Compare (1) source energy to first liquid, (2) first nonzero and first multi-voxel connected liquid times, (3) Ag-liquid volume fraction, (4) main connected-cluster physical volume, centroid/depth and x/y extents, (5) overlap-tracked residence and merge/split history, (6) energy balance, and (7) source-off time to disappearance under each explicitly named boundary condition. Report the whole 3D temperature and liquid fields as well as central ROI metrics. Do not call an adiabatic source-off disappearance an environmental cooling time.

## Comparison table template

| Quantity | Nano | Micro | Matched protocol / uncertainty |
|---|---:|---:|---|
| Ag volume fraction | 0.7000 target | 0.7000 target | measured density/porosity convention |
| Interface area density | pending | pending | same voxel-aware estimator |
| Ag characteristic chord lengths (x/y/z) | pending | pending | distributions and resolution limits |
| Energy to first liquid | pending | pending | source integral and time bracket |
| First multi-voxel connected liquid | pending | pending | six-face, periodic xy/open z |
| Maximum liquid fraction / cluster volume | pending | pending | physical volume, not voxel count |
| Cluster residence / topology changes | pending | pending | overlap-tracked identity, sampled interval |
| Source-off solidification time | pending | pending | state boundary and energy ledger |
| Mesh/domain uncertainty | pending | pending | matched dt, refinement and RVE tests |

## Current limitations

The present calculation uses Case A (the Ti₃SiC₂ phase stores and conducts no heat), idealized Ag internal conduction, an imposed Gaussian/exponential source, and insulated outer faces. It can test geometric localization in that limit; it cannot yet predict an experimental weld threshold or establish a material ranking. Property-source audit, Case B comparison, interfacial thermal resistance, realistic external heat sinks, independently representative micro/nano morphologies and source calibration remain necessary.
