# Ag70 Case A pool-size comparison

This composite places two real snapshots from the same native 1× Case A model side by side. They are at different simulation times and must not be interpreted as simultaneous pools.

- Early small pool: 4.52 ps; 65 connected native voxels at `f_l >= 0.5`; fraction-integrated liquid volume 1.838613 nm³; axis-aligned extent 2.397080 × 1.811145 × 0.951669 nm; peak temperature 1264.383 K.
- Late large main pool: 18.636 ps; 76,392 connected native voxels at `f_l >= 0.5`; fraction-integrated liquid volume 2183.181305 nm³; axis-aligned extent 24.570070 × 24.148599 × 12.371697 nm; peak temperature 2296.532 K.

The liquid-volume integral is `sum(f_l) × voxel_volume` over the thresholded connected component. The thermal model does not resolve an arc, evaporation, liquid flow, or splashing.
