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

The 2× calculation was interrupted by cloud-environment unavailability after the progress log reached 0.38 ps. The last full-field checkpoint is 0.35 ps; after recovery no process was present. The 0.25, 0.30, and 0.35 ps checkpoints each passed a fresh `load_checkpoint()` readback and were pushed. The 0.36–0.38 ps entries are progress-log samples only and are not restart states. At matched 0.25, 0.30, and 0.35 ps, 1×/2× peak temperatures are 774.41/746.06 K, 794.07/763.41 K, and 809.70/777.19 K (2× lower by 28.35, 30.66, and 32.51 K); both grids remain fully solid. After block-averaging 2× to 1× at 0.35 ps, Ag-weighted MAE/RMSE are 0.709/1.954 K, while the hottest-voxel distance from source centre differs by about 12.16 Å and the peak is 32.51 K lower. These are early solid-state mesh-sensitivity results, not convergence and not melt-pool comparisons. The 2× run has not been resumed while the composition/source response is being reviewed. A cooling branch is not applicable while no liquid state exists. The full 3×3 transient remains deferred.

A read-only matched-time audit of the previously melting 41.6654 vol% geometry versus the Ag 69.9996 vol% geometry is in `AG70_NONMELTING_RESPONSE_AUDIT.md`, with reproducible diagnostics in `ag41_vs_ag70_static_diagnostics.json` and a matched 0.2 ps field image. It found about 75% more absorbed source power in the 70 vol% case, but about 68% more Case A Ag heat capacity; the total source-power/heat-capacity ratio differs by only about 4%. Mean Ag temperatures at 0.2 ps are nearly equal, while local peak and liquid fraction differ sharply. The audit points to source placement within distinct connected Ag networks and Case A's zero-conductivity/zero-capacity scaffold as key hypotheses, not a validated composition-only anti-weld result. No transient solve was performed for this audit.

| Item | Status | Result |
|---|---|---|
| First nonzero Ag liquid fraction | 1×, 0.2–1.0 ps | none observed; onset is outside the screened window or absent for this heating input |
| First connected liquid at threshold 0.5 | 1×, 0.2–1.0 ps | none observed |
| First ≥2-voxel six-face cluster | 1×, 0.2–1.0 ps | none observed |
| Source-off solidification | 1× | not applicable: no liquid initial state |
| Matched 2× mesh comparison | interrupted after log sample 0.38 ps; checkpointed through 0.35 ps | no liquid through logged 0.38 ps; full 1.0 ps comparison unavailable; do not infer beyond this interval |
| New-composition 3×3 check | preflight only if time remains | no full transient launched yet |

## Interpretation limits

Case A gives the Ti₃SiC₂ skeleton zero heat capacity and zero conductivity, treats outer faces as insulated, and assumes the current idealized Ag network/source. A fully solid source-off state would mean latent-to-sensible conversion and internal redistribution in an insulated finite domain, not environmental cooling. The property data are not source-verified. A connected liquid cluster is not by itself a stable melt pool, and local melting is not evidence of metallurgical joining or a measured weld threshold.
