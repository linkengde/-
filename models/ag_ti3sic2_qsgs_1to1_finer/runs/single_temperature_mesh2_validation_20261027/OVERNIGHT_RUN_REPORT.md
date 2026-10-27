# Ag–Ti₃SiC₂ single-temperature overnight run

Run date: 2026-10-27 UTC. Target composition: 70:30 Ag:Ti₃SiC₂ by volume. Run directory: `single_temperature_mesh2_validation_20261027`.

## Scope and safety

This work continues the isolated Case A project from the committed 0.2 ps 70:30 checkpoints on branch `codex/single-temperature-study`. It does not write to the original QSGS geometry or the separate double-temperature worktree. The source remains the 5 eV/(ps·Å²), 30 Å lateral Gaussian and 20 Å depth exponential. The 1× continuation is limited to absolute time 1.0 ps and a conservative 2300 K numerical cutoff. It saves after every short segment and terminates if fields become non-finite, liquid fraction leaves [0,1], or relative energy error exceeds 1e-8.

## Work completed before the staged run

- Restored the target-composition 1× and 2× 0.2 ps checkpoints. Both restore with consistent identity, temperature/liquid fields, tracking state and energy ledger.
- Ran the repository suite: 94 tests passed in 4.202 s.
- Added a checkpoint continuation driver, a source-off cooling driver, a geometry audit, and a checkpoint visualization script. All four scripts passed `py_compile`; the snapshot script generated an initial 1× 0.2 ps field image.
- Audited the binary QSGS mask: Ag 69.9996 vol%, one connected phase component per phase under periodic x/y and open z. The Ag–Ti₃SiC₂ interface area is 2.34018×10⁶ Å²; raster interface area density is 0.163125 Å⁻¹. Median directional lineal chords are Ag 14.98/15.09/15.86 Å and Ti₃SiC₂ 8.99/9.06/9.52 Å along x/y/z. These are raster descriptors, not continuum minimum neck widths; the voxel size is about 3 Å.
- Wrote the property-source audit. Exact primary citations are absent in the checkout; a direct Crossref query was blocked by the environment proxy (HTTP 403). All properties therefore remain candidates, not source-verified values.
- Drafted a matched micro/nano comparison plan; no micron-scale geometry or simulation was started.

## Staged heating and cooling results

The 1× run resumed at 0.2 ps and reached 1.0 ps in 904.39 s. Its maximum temperature was 906.2125 K (0.734 Tm), and its Ag liquid fraction, any-liquid voxel count and connected-liquid cluster count remained zero at every sample. First-liquid, first-connected and first multi-voxel-pool times were not observed through the 1 ps cutoff. The final checkpoint was read back successfully. Total injected energy was 1.5513589462365018×10⁻¹⁵ J; the final energy residual was 8.574918039752398×10⁻²⁸ J (relative 5.5274×10⁻¹³). There is no liquid checkpoint to cool, so a source-off “solidification time” is not defined for this 70:30 run.

The 2× calculation is continuing from its 0.2 ps checkpoint toward the same absolute 1.0 ps end time. It uses the unchanged heat source, geometry and candidate properties. No conclusion about 2× through 1.0 ps is made until its checkpoint completes. A cooling branch and full 3×3 transient are deferred because the 1× heating produced no liquid; a 3×3 source/energy preflight will be done only if enough budget remains after the 2× run.

| Item | Status | Result |
|---|---|---|
| First nonzero Ag liquid fraction | 1×, 0.2–1.0 ps | none observed; onset is outside the screened window or absent for this heating input |
| First connected liquid at threshold 0.5 | 1×, 0.2–1.0 ps | none observed |
| First ≥2-voxel six-face cluster | 1×, 0.2–1.0 ps | none observed |
| Source-off solidification | 1× | not applicable: no liquid initial state |
| Matched 2× mesh comparison | running to 1.0 ps | pending completion |
| New-composition 3×3 check | preflight only if time remains | no full transient launched yet |

## Interpretation limits

Case A gives the Ti₃SiC₂ skeleton zero heat capacity and zero conductivity, treats outer faces as insulated, and assumes the current idealized Ag network/source. A fully solid source-off state would mean latent-to-sensible conversion and internal redistribution in an insulated finite domain, not environmental cooling. The property data are not source-verified. A connected liquid cluster is not by itself a stable melt pool, and local melting is not evidence of metallurgical joining or a measured weld threshold.
