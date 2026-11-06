# Same-time main and satellite Ag pools at 15.290000 ps

- Model: Ag70 vol%, native 1× QSGS geometry, single-temperature enthalpy method, Case A.
- Checkpoint: `caseA_ag70_1x_heat_15p29ps.npz`; this image uses one checkpoint only.
- Peak temperature: 2160.974625 K.
- Connected liquid components at `f_l >= 0.5` with face-neighbor connectivity: 2.
- Main pool: 20753 voxels; threshold-mask volume 595.442794 nm³; fraction-integrated liquid volume 593.811416 nm³; dimensions 15.581020 × 14.489159 × 8.565021 nm; centroid [15.277082430251047, 15.78551883242416, 12.709525050065052] nm.
- Largest separate satellite: 15 voxels; threshold-mask volume 0.430378 nm³; fraction-integrated liquid volume 0.426351 nm³; dimensions 1.198540 × 0.905572 × 0.634446 nm; centroid [9.408539, 15.777084529267377, 15.374741400000001] nm.
- Other separate components: 0.

The left panel locates both clusters in the same full-domain state. Its matrix background is sampled once per 4×4×4 native-voxel block, while the liquid clusters are shown at native resolution. The right panel shows the satellite neighborhood at native resolution. These rendering choices do not alter the saved thermal fields.

This thermal model does not include arc-root motion, evaporation, liquid flow, or splashing. The spatially fixed Gaussian source is not a moving arc model.
