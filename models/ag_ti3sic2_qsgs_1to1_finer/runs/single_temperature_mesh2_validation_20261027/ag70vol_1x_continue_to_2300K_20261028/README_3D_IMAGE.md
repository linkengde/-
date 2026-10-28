# Ag70 1x 3D liquid-pool image

These images use the saved 1x Case A heating checkpoint at **18.636 ps**. Orange marks the largest face-connected Ag liquid cluster, defined by `f_liquid >= 0.5`. The 3D pool surface is transparent in the whole-sample view so the surrounding phase point clouds remain visible; the phase clouds are sampled every third voxel for display. The pool itself and the central phase slice use the original QSGS grid.

The phase-aware 3D view renders all four phase classes as solid, color-coded voxel blocks rather than point clouds. It uses one representative native voxel from each 4 × 4 × 4 region for a legible full-volume view, with a corner cutaway exposing the interior phases. The reported pool metrics below are computed from the full native grid, not the display sampling.

At this checkpoint the peak temperature is **2296.53 K**, the Ag-phase liquid fraction is **21.87%**, and the main connected pool contains **76,392 voxels (2191.8 nm³)**. Its axis-aligned dimensions are approximately **24.57 × 24.15 × 12.37 nm** (x × y × z). The horizontal slice at z = 11.26 nm distinguishes Ti3SiC2 scaffold, unmelted Ag, partially liquid Ag, and liquid Ag.

- [Whole sample with transparent pool](ag70_full_model_and_pool_shape_3d_final.png)
- [3D cutaway with all phases as voxel blocks](ag70_phase_map_3d_18p636ps.png)
- [Pool shape close-up](ag70_pool_shape_3d_final.png)
- [Phase slice through pool center](ag70_phase_slice_through_pool_final.png)
- [Rendering script](render_current_3d.py)
- [Phase-aware 3D rendering script](render_phase_map_3d.py)
- [Source checkpoint](checkpoints/caseA_ag70_1x_heat_18p636ps.npz)

Recreate the figures from this directory with:

```bash
MPLCONFIGDIR=/tmp/mpl-ag70 python render_current_3d.py
MPLCONFIGDIR=/tmp/mpl-ag70 python render_phase_map_3d.py
```

This is a thermal enthalpy visualization. The model has no fluid momentum, free-surface, or vaporization treatment, so it does not depict splashing or material ejection. The run-specific 2300 K stop guard is numerical protection for a model without vaporization; it is not a general physical standard.
