# Ag–Ti₃SiC₂ morphology candidates

This directory contains three reproducible, geometry-only two-phase voxel candidates in a 150.1 × 300 × 300 Å box. Blue is Ti₃SiC₂ and light gray is Ag. The target is 40/60 wt% Ti₃SiC₂/Ag, converted to volume fraction using the supplied dense-phase densities.

The shapes use a seeded, correlated random field thresholded to the exact target voxel fraction. The correlation length was chosen to keep branches at a medium thickness. This is an exploratory geometric design, not a calibrated phase-field prediction, a physical spinodal decomposition calculation, or an atomistic model.

Run with the prepared environment:

```bash
XDG_CACHE_HOME=/workspace/.cache MPLCONFIGDIR=/workspace/.cache/matplotlib FONTCONFIG_FILE=/workspace/.setup/phasefield/fonts.conf \
  /workspace/.venvs/phasefield/bin/python generate_geometry.py
```

Each `candidate_NN` directory includes `morphology_reference_view_2d.png` (300 × 300 px, central YZ section), `geometry_preview.png` (3D view and three sections), `phase_masks_and_field.npz`, and `geometry_stats.json`. The masks define the phases; PNG colors are for visual review. `geometry_audit.md` records assumptions, volume conversion, connectivity, thickness proxies, and limitations.
