# Ag–W V12 fine QSGS scaffold preview

## Purpose

This separate candidate tests the finer W-backbone scale requested from the supplied Ag–Ti₃SiC₂ reference image. It keeps the 200 × 200 × 200 Å material, 200 × 200 × 300 Å box, x/y periodic and z free, and approximately 1:1 Ag:W mass ratio. V11 and V5–V10 files are unchanged.

## Geometry

- QSGS seed: 20261125
- Core seed probability: 0.0015
- Base growth probability: 0.1
- Maximum generations: 8
- Target W voxel fraction: 0.3545
- Voxel grid: 80³ at 2.5 Å
- Atoms: 482,289 (178,353 W; 303,936 Ag)
- Ag:W mass ratio: 0.999897
- W strut local-diameter median: 5.00 Å; 90th percentile: 8.66 Å
- Atomic connectivity: W 98.35%; Ag 99.75%

The current V11 single-period 200 Å Gyroid has a much coarser W labyrinth (approximate median local diameter 25.7 Å).

## Status and limitation

**Geometry preview only; not an MD-ready structure.** The initial Ag/W lattice placement still has a minimum cross-phase distance of 0.08 Å. Static relaxation, 300 K stability, TTM mapping, and melt/quench validation have not been run. The image is for judging morphology and scaffold thickness only.

## Files

- `qsgs_fine_candidate/ag_w_qsgs.data.gz`: candidate atomistic geometry
- `qsgs_fine_candidate/qsgs_morphology.npz`: fine W/Ag voxel masks
- `qsgs_fine_candidate/geometry_report.json`: generated geometry statistics
- `qsgs_fine_candidate/ag_w_qsgs.extxyz.gz`: same candidate for ASE-compatible visualization
- `qsgs_fine_candidate/in.geometry_read_check`: LAMMPS data-file read check
- `render_v12_refined_qsgs_ovito.py`: OVITO rendering script
- `../../../scripts/build_ag_w_qsgs.py`: shared four-parameter QSGS builder; this candidate used seed 20261125, p_core 0.0015, p_grow 0.10, 8 generations, W fraction 0.3545
- `v12_fine_qsgs_morphology_preview.png`: full and cutaway OVITO views
