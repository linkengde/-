# Continuum applicability for the Ag70 pore geometry

## Geometry record

The read-only Stage 1 audit examined the Ag70 QSGS mask with SHA256 `03ef47b1dd534b3215cf03f42c8cf7e887dcdfcc185b33bab7453cd838c1faac`. The mask is 100×100×50 cells over approximately 29.9635×30.1857×15.8612 nm. The cell widths are about 0.300×0.302×0.317 nm. Its estimated widest bottom-to-top Ag path has a geometric bottleneck diameter near 1.55 nm along z, or roughly 4.9 z cells. This is a voxel-threshold topology estimate, not a directly resolved hydraulic aperture or atomically reconstructed channel.

## Scale assessment

The voxel width is on the order of an atomic spacing. A rough atomic-volume length for liquid Ag near melting, computed from the experimental density correlation, is about 0.27 nm. Thus the 1.55 nm bottleneck spans only about six such characteristic lengths. This comparison is an order-of-magnitude indicator; atomic-volume length is not a molecular mean free path or a liquid correlation length.

For a dense liquid, a gas-style Knudsen number is not the right continuum test. A credible pore-scale continuum description needs scale separation between the channel width and microscopic structural/relaxation lengths, and boundary conditions that represent the actual solid surface. Here neither an Ag/Ti₃SiC₂ slip length nor a dynamic wall-friction law is established, and the binary QSGS mask does not encode atomic surface termination, roughness, adsorbates, or local reconstruction. The ratio b/w therefore cannot be evaluated.

The current ~0.30 nm grid also leaves only about five cells across the z bottleneck. Refining the numerical mesh over the same binary geometry cannot create missing atomistic interface information or validate the pore geometry itself. The successful P0/P1 VOF tests establish numerical behavior on controlled continuum benchmarks; they do not establish quantitative Navier–Stokes validity in these nanothroats.

## Model choice gate

**Quantitative pore-scale Navier–Stokes in the existing Ag70 mask is not justified.** A coarse-grained Darcy/Brinkman route may be considered only after the permeability tensor and its temperature/liquid-fraction dependence are measured or independently computed with a validated pore model, and after interface momentum transfer is calibrated. A separate atomistic or mesoscopic study may be needed to derive the channel-scale transport closure, but that is a future calibration task, not a result of this audit.

No pressure gradient, body force, or flow rate is supplied by the present geometry record. Do not invent one to create visible flow. No Ag70 flow was simulated in this phase.
