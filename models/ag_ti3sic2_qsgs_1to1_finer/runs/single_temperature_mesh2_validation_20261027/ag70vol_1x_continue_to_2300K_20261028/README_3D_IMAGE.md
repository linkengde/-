# Ag70 1x 3D liquid-pool image

The images show the saved 1x Case A heating state at **18.286 ps**. The orange surface is the largest face-connected Ag liquid cluster, defined by liquid fraction `f_liquid >= 0.5`, extracted from the original QSGS voxel geometry.

At this checkpoint the peak temperature is 2285.37 K, the Ag-phase liquid fraction is 19.21%, and the main connected pool contains 67,099 voxels (1925.2 nm³). Its axis-aligned dimensions are approximately 23.37 × 23.24 × 12.05 nm (x × y × z). The full-sample phase clouds are sampled every third voxel for display; the highlighted pool uses the original grid.

- [Full sample and pool location](ag70_full_model_and_pool_shape_3d_18p286ps.png)
- [Pool shape close-up](ag70_pool_shape_3d_18p286ps.png)
- [Rendering script](render_current_3d.py)
- [Source checkpoint](checkpoints/caseA_ag70_1x_heat_18p286ps.npz)

Recreate the figures from this directory with:

```bash
MPLCONFIGDIR=/tmp/mpl-ag70 python render_current_3d.py
```

This is a thermal enthalpy visualization. The model has no fluid momentum, free-surface, or vaporization treatment, so the figure does not depict splashing or material ejection. The run-specific 2300 K stop guard is numerical protection for a model without vaporization; it is not a general physical standard.
