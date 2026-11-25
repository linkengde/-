# Ag–W V14 finer QSGS scaffold candidate

## Result

V14 increases QSGS seed density and uses a finer 120³ voxel grid. The W scaffold is visibly finer than V13 while preserving a bicontinuous geometry and approximately 1:1 Ag:W mass ratio.

- QSGS seed: 20261129
- Core seed probability: 0.005
- Base growth probability: 0.07
- Maximum generations: 12 (target reached at generation 7)
- W target voxel fraction: 0.35437
- Grid: 120³; voxel spacing: 1.6667 Å
- Material: 200 × 200 × 200 Å; box: 200 × 200 × 300 Å; boundaries: p p f
- Atoms: 481,648 (178,106 W; 303,542 Ag)
- Ag/W mass ratio: 0.999985405 (relative error 0.00146%)
- W largest atomic connected component: 98.134%; XY winding rank 2; spans full z: True
- Ag largest atomic connected component: 99.908%; XY winding rank 2; spans full z: True
- W mask local diameter estimate q10/q25/q50/q75/q90/q95: 3.33, 3.33, 3.33, 3.33, 4.71, 4.71 Å
- Relative to V13, W mask q90 is reduced by 33.3% (V13 q90 ≈ 7.07 Å).

## Status and limitations

**Geometry-only candidate; not MD-ready.** The minimum initial Ag–W distance is 0.08 Å, so cross-phase overlap remains. No EFS minimization, 300 K MD, or TTM calculation has been performed.

The local-diameter numbers are voxel-mask estimates at 1.667 Å resolution, not physical atom-resolved strut diameters. The W atomic largest-component share is only modestly above the 98% geometry screen, so connectivity should be rechecked after any relaxation.

## Files

- `qsgs_v14_finer_candidate/ag_w_qsgs.data.gz`: LAMMPS atomic geometry candidate
- `qsgs_v14_finer_candidate/ag_w_qsgs.extxyz.gz`: visualization copy
- `qsgs_v14_finer_candidate/qsgs_morphology.npz`: W and Ag voxel masks and generation parameters
- `qsgs_v14_finer_candidate/geometry_report.json`: atom counts, mass ratio, topology, and QSGS settings
- `render_v14_finer_qsgs_ovito.py`: OVITO render script
- `v14_finer_qsgs_morphology_preview.png`: full model and cutaway preview

Rebuild:

```bash
python research/ag_w_reference_pack/scripts/build_ag_w_qsgs.py --output <output> --seed 20261129 --p-core 0.005 --p-grow 0.07 --generations 12 --target-w-fraction 0.35437 --grid 120
```
