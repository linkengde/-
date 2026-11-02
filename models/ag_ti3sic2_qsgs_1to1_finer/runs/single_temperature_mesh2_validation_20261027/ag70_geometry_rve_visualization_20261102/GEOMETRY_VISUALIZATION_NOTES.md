# Ag70 current geometry visualization

This figure renders the existing `qsgs_morphology_ag70vol.npz`; it does not regenerate or edit the QSGS geometry.

- Grid: `100 × 100 × 50`.
- Domain: `29.9635 × 30.1857 × 15.8612 nm`.
- Voxel spacing: `0.29963 × 0.30186 × 0.31722 nm`.
- Volume fractions: Ag `69.9996%`; Ti₃SiC₂ `30.0004%`.
- Local inscribed-sphere diameter proxy on periodic-x/y EDT maxima: P10 `0.599 nm`, median `0.851 nm`, P90 `1.403 nm`, P95 `1.490 nm`, max `1.995 nm`.
- The EDT maxima are an approximate local-thickness statistic, not a uniquely segmented branch diameter or measured physical dimension. With approximately 0.30 nm voxels, a 0.60 nm estimate is only about two voxels and is resolution-limited.
- The 0.2 µm and 2 µm markers are reference scales from the user's supplied conceptual figure, not dimensions in this current nanoscale model.

The x/y periodicity is accounted for by tiling the mask 3×3 before the distance transform. The top and bottom two voxel planes are excluded from thickness statistics. The mesh in the 3-D panels is decimated only for display; all calculations use the original mask.
