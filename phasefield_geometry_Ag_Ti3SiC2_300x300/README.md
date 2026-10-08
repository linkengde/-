# Branched Ag–Ti₃SiC₂ geometry candidates

This directory contains three reproducible, geometry-only two-phase voxel candidates at 150.1 × 300 × 300 Å. They target 40/60 wt% Ti₃SiC₂/Ag using the supplied room-temperature dense-phase densities. They are not calibrated phase-field predictions or atomistic models.

Run with the prepared environment:

```bash
XDG_CACHE_HOME=/workspace/.cache MPLCONFIGDIR=/workspace/.cache/matplotlib FONTCONFIG_FILE=/workspace/.setup/phasefield/fonts.conf \
  /workspace/.venvs/phasefield/bin/python generate_geometry.py
```

Read `geometry_audit.md` for assumptions, connectivity, thickness proxies, and limits. Candidate data are in each `candidate_NN/phase_masks_and_field.npz`; each preview is `candidate_NN/geometry_preview.png`.
