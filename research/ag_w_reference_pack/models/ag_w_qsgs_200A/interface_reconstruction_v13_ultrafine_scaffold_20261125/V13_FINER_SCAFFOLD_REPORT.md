# Ag–W V13 finer QSGS scaffold preview

## Geometry

V13 uses denser QSGS core seeding and keeps the Ag:W mass ratio near 1:1. It is separate from V12.

- Seed: 20261128
- Core seed probability: 0.0025
- Base growth probability: 0.1
- Maximum generations: 8
- W target voxel fraction: 0.355
- Grid: 80³, voxel spacing 2.5 Å
- Material: 200 × 200 × 200 Å; simulation box: 200 × 200 × 300 Å; boundaries p p f
- Atoms: 482,375 (178,382 W; 303,993 Ag)
- Ag:W mass ratio: 0.999922
- W largest atomic component: 98.19%; XY winding rank 2; spans full z
- Ag largest atomic component: 99.78%
- W local diameter quantiles (10/25/50/75/90%): 5.00, 5.00, 5.00, 7.07, 7.07 Å

Compared with V12, the W-mask 90th-percentile local diameter decreases from 8.66 Å to 7.07 Å (18.4%). The median remains 5 Å, the limit set by this 2.5 Å voxel grid.

## Status

**Geometry preview only; not an MD-ready structure.** The initial minimum Ag–W separation is 0.08 Å, so strong overlap remains. No static relaxation, 300 K MD, or TTM validation has been run.

## Files

- `qsgs_finer_massbalanced_candidate/ag_w_qsgs.data.gz`: LAMMPS atomistic candidate
- `qsgs_finer_massbalanced_candidate/ag_w_qsgs.extxyz.gz`: visualization copy
- `qsgs_finer_massbalanced_candidate/qsgs_morphology.npz`: W and Ag voxel masks
- `qsgs_finer_massbalanced_candidate/geometry_report.json`: geometry and connectivity report
- `render_v13_finer_qsgs_ovito.py`: OVITO rendering script
- `v13_finer_qsgs_morphology_preview.png`: full model and cutaway views

Rebuild command:

```bash
python research/ag_w_reference_pack/scripts/build_ag_w_qsgs.py --output <output> --seed 20261128 --p-core 0.0025 --p-grow 0.10 --generations 8 --target-w-fraction 0.3550 --grid 80
```
