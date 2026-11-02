# Ag70 scaffold zoom-in panels

Both inset windows are extracted from the same original `qsgs_morphology_ag70vol.npz` QSGS geometry.
They are not separate micro- and submicron structures.

- Geometry grid: `100 × 100 × 50`.
- Domain: `29.9635 × 30.1857 × 15.8612 nm`.
- Voxel spacing: `0.29963 × 0.30186 × 0.31722 nm`.
- Thin-region center voxel (x,y,z): `(9, 49, 39)`; selected local diameter proxy `0.5993 nm` (P10 target `0.5993 nm`).
- Thick/junction-region center voxel (x,y,z): `(9, 10, 13)`; selected local diameter proxy `1.4033 nm` (P90 target `1.4033 nm`).
- Each crop is `19 × 19 × 17` voxels; x/y crops wrap periodically.

The inset double arrows show twice the Euclidean distance-transform value at a local maximum inside the scaffold. This is a voxel-based proxy for local scaffold thickness, not a direct measurement of branch diameter or pore spacing. At approximately 0.30 nm per voxel, the lower-tail feature is only about two voxels wide and is resolution-limited. Mesh decimation and smoothing are applied only to the rendering; all reported geometry statistics come from the original voxel mask.
