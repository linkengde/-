# Ag–Ti₃SiC₂ melt-pool geometry candidate 02

This folder isolates the user-selected candidate 02 from the geometry-candidate package. It is a reproducible **geometry-only** morphology design, not a calibrated phase-field prediction, atomic model, or validated melt-pool simulation input.

## Geometry

- Box: 150.1 × 300 × 300 Å (X × Y × Z).
- Grid: 76 × 150 × 150 voxels; spacing: 1.975 × 2 × 2 Å.
- Target composition: Ti₃SiC₂:Ag = 4:6 by mass. Using dense-phase densities of 4.53 and 10.4976 g/cm³ gives 60.7057 vol% Ti₃SiC₂ and 39.2943 vol% Ag. The voxelized result is 40.0000 wt% Ti₃SiC₂ to the reported precision.
- Generation: thresholded 3D correlated random field, seed 20261009, smoothing σ = (1.7, 2.7, 2.7) voxels. These are shape controls, not material parameters.
- Boundary: nonperiodic box with reflection at the field-generation boundary; porosity 0%.
- Connectivity: largest Ti₃SiC₂ component contains 99.9881% of that phase and spans X/Y/Z; largest Ag component contains 99.7308% and spans X/Y/Z (26-neighbour connectivity).
- Ti₃SiC₂ thickness proxy median: 8.855 Å. This is a geometric distance-transform measure, not a crystallographic thickness.

Density references recorded in the supplied material-input note: [Ti₃SiC₂](https://doi.org/10.1016/S1359-6462(99)00054-8) and [Ag](https://doi.org/10.1007/BF00917194).

## Files

- `candidate_02/geometry_3d_view.png`: 3D view, blue Ti₃SiC₂ and gray Ag.
- `candidate_02/geometry_preview.png`: 3D view and orthogonal slices.
- `candidate_02/morphology_reference_view_2d.png`: 2D slice view.
- `candidate_02/phase_masks_and_field.npz`: signed geometry field, two phase masks, voxel spacing, and box lengths.
- `candidate_02/geometry_stats.json`: machine-readable statistics and generation metadata.
- `selected_candidate.json`: selected-candidate record.
- `generate_selected_candidate.py`: standalone deterministic generator for candidate 02 only.
- `SHA256SUMS.txt`: checksums for the package files (excluding the checksum file itself).

## Regenerate

The script requires Python 3, NumPy, SciPy, and Matplotlib. In the prepared cloud environment:

```bash
XDG_CACHE_HOME=/workspace/.cache \
MPLCONFIGDIR=/workspace/.cache/matplotlib \
FONTCONFIG_FILE=/workspace/.setup/phasefield/fonts.conf \
/workspace/.venvs/phasefield/bin/python generate_selected_candidate.py
```

It writes the selected output under `candidate_02/` and refreshes `selected_candidate.json`.

## Limits

The masks do not define crystal lattices, atomic coordinates, crystal orientation, interface registry, or surface termination. The current geometry fills the box and has no vacuum region. It is not ready for LAMMPS or TTM-MD; a supported crystal structure, an independently validated Ag/Ti/Si/C interaction model, and Ag/Ti₃SiC₂ two-temperature parameters are still needed before a melt-pit simulation.
